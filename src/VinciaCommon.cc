// VinciaCommon.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the headers) for the Colour, Rambo, and
// VinciaCommon classes, and related auxiliary methods.

#include "Pythia8/VinciaCommon.h"

namespace Pythia8 {

//==========================================================================

// Rambo phase space generator.

//--------------------------------------------------------------------------

// Massless flat phase space generator. Generate a random (uniformly
// distributed) massless PS point with nOut particles and total sqrt(s) = eCM.

double Rambo::genPoint(double eCM, int nOut, vector<Vec4>& pOut) {

  // Set size of output vector
  pOut.resize(nOut);
  // Create momentum-sum four-vector
  Vec4 R;
  // Generate nParticles independent massless 4-momenta with isotropic angles
  for (int i = 0; i < nOut; ++i) {
    // Cos(theta), sin(theta), and phi
    double c   = 2.0*rndmPtr->flat() - 1.0;
    double s   = sqrt(1.0-pow2(c));
    double phi = 2.0*M_PI*rndmPtr->flat();
    // Norm
    double r12 = 0.0;
    while (r12 == 0.0) {
      double r1 = rndmPtr->flat();
      double r2 = rndmPtr->flat();
      r12 = r1*r2;
    }
    double En = -log(r12);
    pOut[i].e(En);
    pOut[i].pz(En*c);
    pOut[i].py(En*s*cos(phi));
    pOut[i].px(En*s*sin(phi));
    // Add to vector and add to sum
    R += pOut[i];
  }
  // Compute ECM and normalise to unity (with sign flip)
  double Rmass = R.mCalc();
  R /= -Rmass;
  // Transform momenta so add up to (eCM, 0, 0, 0)
  double a = 1.0/(1.0-R.e());
  double x = eCM/Rmass;
  for (int i = 0; i < nOut; ++i) {
    double bq = dot3(R, pOut[i]);
    pOut[i].px( x * (pOut[i].px()+R.px()*(pOut[i].e()+a*bq)) );
    pOut[i].py( x * (pOut[i].py()+R.py()*(pOut[i].e()+a*bq)) );
    pOut[i].pz( x * (pOut[i].pz()+R.pz()*(pOut[i].e()+a*bq)) );
    pOut[i].e(  x * (-R.e()*pOut[i].e()+bq) );
  }
  // The weight is always unity for the massless algorithm.
  return 1.0;

}

//--------------------------------------------------------------------------

// Massive flat phase space generator, generalised according to the
// original paper. The momenta are not distributed flat in phase
// space anymore, returns the weight of the phase space configutation.

double Rambo::genPoint(double eCM, vector<double> mIn, vector<Vec4>& pOut) {

  // Call the massless genPoint, initializing weight.
  int nOut = mIn.size();
  if (nOut <= 1 || eCM <= 0.) return 0.;
  double weight = genPoint(eCM, nOut, pOut);
  bool massesnonzero = false;

  // Set up the function determining the rescaling parameter xi.
  vector<double> energies;
  for (int i = 0; i < nOut; i++) {
    energies.push_back(pOut[i].e());
    if (pow2(mIn[i]/eCM) > TINY) massesnonzero = true;
  }

  // If none of the reduced masses is > TINY, return
  if (!massesnonzero) return weight;

  // Rescale all the momenta.
  TXiFunctor rhs = TXiFunctor(mIn, energies);
  double xi = zbrent(rhs, eCM, 0., 1., 1e-10);
  for (int iMom = 0; iMom < nOut; iMom++) {
    pOut[iMom].rescale3(xi);
    pOut[iMom].e( sqrt(pow2(mIn[iMom]) + pow2(xi)*pow2(pOut[iMom].e())) );
  }

  // Determine the quantities needed for the calculation of the weight.
  double sumP(0.), prodPdivE(1.), sumP2divE(0.);
  for (int iMom = 0; iMom < nOut; iMom++) {
    double pAbs2 = pOut[iMom].pAbs2();
    double pAbs  = sqrt(pAbs2);
    sumP      += pAbs;
    prodPdivE *= pAbs/pOut[iMom].e();
    sumP2divE += pAbs2/pOut[iMom].e();
  }

  // There's a typo in eq. 4.11 of the Rambo paper by Kleiss, Stirling
  // and Ellis, the Ecm below is not present there.
  weight *= pow(sumP/eCM,2*nOut-3)*prodPdivE*eCM/sumP2divE;
  return weight;
}

//==========================================================================

// The Colour class.

//--------------------------------------------------------------------------

// Initialize.

bool Colour::init() {

  // Sanity check.
  if (!isInitPtr) return false;

  // Set verbosity level.
  verbose = settingsPtr->mode("Vincia:verbose");

  // Set parameters. CR disabled in this version.
  inheritMode       = settingsPtr->mode("Vincia:CRinheritMode");

  // Sett initialization.
  isInit = true;
  return isInit;

}

//--------------------------------------------------------------------------

// Map a list of particles with ordinary Les Houches colour tags into
// a new set of Les Houches colour tags with the last index running
// from 1-9.
//
// Note on coloured resonances. These should have their colours
//  defined by whatever system they were produced by, and should not
//  be recoloured when decaying. This means that only colour lines for
//  which both the colour AND the corresponding anticolour line are
//  found inside the same system should be redefined (assuming the
//  coloured resonance itself does not appear explicitly in its
//  corresponding decay partonSystem.

bool Colour::colourise(int iSys, Event& event) {

  // Sanity checks.
  if (!isInit) {
    printOut("Colour::colourise","ERROR! Colour not initialised");
    return false;
  }
  else if (partonSystemsPtr->sizeAll(iSys) <= 1) return false;

  // Construct colour map and assign new tags.
  map<int, int> colourMap;
  int startTag = event.lastColTag();
  int nTries = 0;
  bool accept = false;

  // Do not recolour decaying resonances (as they will already have
  // been coloured when they were produced).
  int colRes(0), acolRes(0);
  if (partonSystemsPtr->hasInRes(iSys)) {
    int iRes = partonSystemsPtr->getInRes(iSys);
    colRes  = event[iRes].col();
    acolRes = event[iRes].acol();
  }

  // Reject assignments that would produce (subleading) singlet gluons.
  while (!accept && ++nTries < 10) {
    colourMap.clear();

    // Example: if startTag is 220-229, nextTagBase is 230.
    int nextTagBase = 10*int((startTag/10)+1);
    for (int i=0; i<partonSystemsPtr->sizeAll(iSys); ++i) {
      int i1 = partonSystemsPtr->getAll(iSys,i);
      if (i1 <= 0) continue;
      Particle* partonPtr = &event[i1];
      if (partonPtr->colType() == 0) continue;
      int col, acol;

      // Cross initial-state colours.
      if (i < partonSystemsPtr->sizeAll(iSys)
        - partonSystemsPtr->sizeOut(iSys)) {
        acol = partonPtr->col();
        col  = partonPtr->acol();
        if (col  == acolRes)  col = 0;
        if (acol ==  colRes) acol = 0;
      }
      else {
        col  = partonPtr->col();
        acol = partonPtr->acol();
        if (col == colRes)   col  = 0;
        if (acol == acolRes) acol = 0;
      }
      if (col == 0 && acol == 0) continue;
      int colIndx = colourMap[col];
      int acolIndx = colourMap[acol];
      if (col != 0) {

        // First time a tag is encountered, mark it negative -> one end found.
        if (colIndx == 0) {
          // Ensure gluons have different col and acol indices.
          while (colIndx == 0 || colIndx == acolIndx) {
            colIndx = nextTagBase + int(rndmPtr->flat()*9) + 1;
          }
          colourMap[col]  = -colIndx;
        }
        // Second time mark it positive -> both ends found
        else colourMap[col] = abs(colourMap[col]);
      }

      if (acol != 0) {
        // First time a tag is encountered, mark it negative -> one end found
        if (acolIndx == 0) {
          // Ensure gluons have different col and acol indices
          while (acolIndx == 0 || colIndx == acolIndx) {
            acolIndx = nextTagBase + int(rndmPtr->flat()*9) + 1;
          }
          colourMap[acol] = -acolIndx;
        }
        // Second time mark it positive -> both ends found
        else colourMap[acol] = abs(colourMap[acol]);
      }
      // Update nextTagBase
      nextTagBase += 10;
    }

    // Check if these assignments would produce any singlet gluons
    accept = true;
    for (int i=0; i<partonSystemsPtr->sizeAll(iSys); ++i) {
      int i1 = partonSystemsPtr->getAll(iSys,i);
      Particle* partonPtr = &event[i1];
      if (partonPtr->colType() != 2) continue;
      int colIndexNew  = colourMap[partonPtr->col()] % 10;
      int acolIndexNew = colourMap[partonPtr->acol()] % 10;
      if (colIndexNew == acolIndexNew) {
        accept=false;
        break;
      }
    }
  }

  // Check for failure to find acceptable conf.
  if (!accept) {
    if (verbose >= 3) event.list();
    printOut("Colour::colourise","Warning! failed to avoid singlet gluon(s).");
  }

  // Update event.
  for (int i = 0; i < partonSystemsPtr->sizeAll(iSys); ++i) {
    int ip = partonSystemsPtr->getAll(iSys,i);
    Particle* partonPtr = &event[ip];
    if (partonPtr->colType() == 0) continue;
    if ( colourMap[partonPtr->col()] > 0 )
      partonPtr->col(colourMap[partonPtr->col()]);
    if ( colourMap[partonPtr->acol()] > 0 )
      partonPtr->acol(colourMap[partonPtr->acol()]);

    // Update max used colour tag.
    int lastTag = event.lastColTag();
    int colMax  = max(abs(partonPtr->col()),abs(partonPtr->acol()));
    while (colMax > lastTag) lastTag = event.nextColTag();
  }

  // Return successful.
  return true;
}

//--------------------------------------------------------------------------

// Order a list of partons in colour sequence.

vector<int> Colour::colourSort(vector<Particle*> partons) {

  // Output vector (starts empty).
  vector<int> iSorted;

  // Shorthand for final-state parton multiplicities.
  int nPartons=partons.size();
  if (nPartons <= 1) return iSorted;

  // Find string endpoints and colour types of particles.
  vector<int> iTrip, iAnti, iOct, iOtherIn, iOtherOut;

  // Definition of colType (classified by multiplet up to total charge
  // p+q = 4).
  //  +- 1 : triplet (e.g., Red)        [1,0] / [0,1] {quark, antidiquark}
  //     2 : octet (e.g., R-Gbar)       [1,1] {
  //           gluon, incoherent qqbar, partially coherent gluon-gluon
  //           eg R-(Bbar-B)-Gbar -> R-Gbar (no junction)
  //           or (R-B)-(Bbar-Gbar) -> Gbar-R (junction-antijunction)}
  //  +- 3 : sextet (e.g., 2 x Red)     [2,0] / [0,2] {2 incoherent quarks}
  //  +- 4 : fifteen (e.g., R-R-Gbar)   [2,1] / [1,2] {incoherent qg}
  //  +- 5 : decuplet (e.g., 3 x Red)   [3,0] / [0,3] {
  //           3 incoherent quarks / partially coherent gluon-gluon
  //           eg R-Gbar-R-Bbar -> R-R-R}
  //     6 : vigintiseptet (e.g., 2 x R-Gbar)   [2,2] {2 incoherent gluons}
  //  +- 7 : fifteen' (e.g., 4 x Red)   [4,0] / [0,4] {4 incoherent quarks}
  //  +- 8 : vigintiquartet (e.g., R-R-R-Gbar)  [3,1] / [1,3]

  map<int, int> iOfAcol;
  for (int i=partons.size()-1; i>=0; --i) {
    int sign    = (partons[i]->isFinal() ? 1 : -1);
    int colType = particleDataPtr->colType(partons[i]->id());

    // Store indices of anticolour partners.
    if (sign == 1 && ( colType == -1 || colType == 2))
      iOfAcol[partons[i]->acol()] = i;
    else if (sign == -1 && ( colType == 1 || colType == 2 ))
      iOfAcol[partons[i]->col()] = i;

    // Construct list of triplets (= starting points).
    if (colType * sign == 1) iTrip.push_back(i);

    // List of antitriplets.
    else if (colType * sign == -1) iAnti.push_back(i);
    // List of octets.
    else if (colType == 2) iOct.push_back(i);

    // Higher representations.
    else if (abs(colType) >= 3) {
      cout << "colourSort(): ERROR! handling of coloured particles in "
           << "representations higher than triplet or octet is not implemented"
           << endl;
    }

    // Colourless particles.
    else if (sign == -1) iOtherIn.push_back(i);
    else iOtherOut.push_back(i);
  }

  // Now sort particles.
  int  i1 = -1;
  bool beginNewChain = true;

  // Keep looping until we have sorted all particles.
  while (iSorted.size() < partons.size()) {

    // Start new piece (also add colourless particles at front and end).
    if (beginNewChain) {

      // Insert any incoming colourless particles at front of iSorted.
      if (iOtherIn.size() > 0) {
        iSorted.push_back(iOtherIn.back());
        iOtherIn.pop_back();

      // Triplet starting point (charge += 1).
      } else if (iTrip.size() > 0) {
        beginNewChain = false;
        iSorted.push_back(iTrip.back());
        iTrip.pop_back();

      // Octet starting point if no triplets/sextets available.
      } else if (iOct.size() > 0) {
        beginNewChain = false;
        iSorted.push_back(iOct.back());
        iOct.pop_back();

      } else if (iOtherOut.size() > 0) {
        iSorted.push_back(iOtherOut.back());
        iOtherOut.pop_back();
      }

      // Index of current starting parton.
      i1 = iSorted.back();

    // Step to next parton in this chain.
    } else {
      bool isFinal = partons[iSorted.back()]->isFinal();
      int col = (isFinal) ? partons[iSorted.back()]->col()
        : partons[iSorted.back()]->acol();
      int iNext = iOfAcol[col];

      // Sanity check.
      if (iNext < 0) {
        cout << "colourSort(): ERROR! cannot step to < 0" << endl;
        beginNewChain = true;

      // Catch close of gluon ring.
      } else if (iNext == i1) {
        beginNewChain = true;

      // Step to next parton; end if not gluon (antiquark or IS quark).
      } else {
        // Add to sorted list.
        iSorted.push_back(iNext);
        // If endpoint reached, begin new chain.
        if (particleDataPtr->colType(partons[iNext]->id()) != 2)
          beginNewChain = true;
        // Octet: continue chain and erase this octet from list.
        else {
          beginNewChain = false;
          // Erase this endpoint from list.
          for (int i=0; i<(int)iOct.size(); ++i) {
            if (iOct[i] == iNext) {
              iOct.erase(iOct.begin()+i);
              break;
            }
          }
        }
      }
    } // End step to next parton.
  }

  // Normal return.
  return iSorted;

}

//--------------------------------------------------------------------------

// Make colour maps and construct list of parton pairs that form QCD dipoles.

void Colour::makeColourMaps(const int iSysIn, const Event& event,
  map<int,int>& indexOfAcol, map<int,int>& indexOfCol,
  vector< pair<int,int> >& antLC, const bool findFF, const bool findIX) {

  // Loop over all parton systems.
  int iSysBeg = (iSysIn >= 0) ? iSysIn : 0;
  int iSysEnd = (iSysIn >= 0) ? iSysIn + 1: partonSystemsPtr->sizeSys();
  for (int iSys = iSysBeg; iSys < iSysEnd; ++iSys) {

    // Loop over a single parton system.
    int sizeSystem = partonSystemsPtr->sizeAll(iSys);
    for (int i = 0; i < sizeSystem; ++i) {
      int i1 = partonSystemsPtr->getAll( iSys, i);
      if ( i1 <= 0 ) continue;

      // Save to colour maps.
      int col  = event[i1].col();
      int acol = event[i1].acol();

      // Switch colours for initial partons.
      if (!event[i1].isFinal()) {
        col  = acol;
        acol = event[i1].col();
      }

      // Save colours (taking negative-index sextets into account).
      if (col > 0) indexOfCol[col] = i1;
      else if (col < 0) indexOfAcol[-col] = i1;
      if (acol > 0) indexOfAcol[acol] = i1;
      else if (acol < 0) indexOfCol[-acol] = i1;

      // Look for partner on colour side.
      if (col > 0 && indexOfAcol.count(col) == 1) {
        int i2 = indexOfAcol[col];
        if ( event[i1].isFinal() && event[i2].isFinal() ) {
          if (findFF) antLC.push_back( make_pair(i1,i2) );
        } else if (findIX) antLC.push_back( make_pair(i1,i2) );
      }

      // Look for partner on anticolour side.
      if (acol > 0 && indexOfCol.count(acol) == 1) {
        int i2 = indexOfCol[acol];
        // Coloured parton first -> i2, i1 instead of i1, i2)
        if (event[i1].isFinal() && event[i2].isFinal()) {
          if (findFF) antLC.push_back( make_pair(i2, i1) );
        } else if (findIX) antLC.push_back( make_pair(i2,i1) );
      }

      // Allow for sextets: negative acol -> extra positive col.
      if (acol < 0 && indexOfAcol.count(-acol) == 1) {
        int i2 = indexOfAcol[-acol];
        if (event[i1].isFinal() && event[i2].isFinal()) {
          if (findFF) antLC.push_back( make_pair(i1,i2) );
        } else if (findIX) antLC.push_back( make_pair(i1,i2) );
      }
      if (col < 0 && indexOfCol.count(-col) == 1) {
        int i2 = indexOfAcol[-acol];
        if (event[i1].isFinal() && event[i2].isFinal()) {
          if (findFF) antLC.push_back( make_pair(i1,i2) );
        } else if (findIX) antLC.push_back( make_pair(i1,i2) );
      }
    }
  }
  return;

}

//--------------------------------------------------------------------------

// Determine which of two antennae inherits the old colour tag after a
// branching. Default is that the largest invariant has the largest
// probability to inherit, with a few alternatives also implemented.

bool Colour::inherit01(double s01, double s12) {
  // Initialization check.
  if (!isInit) {
    printOut("VinciaColour::inherit01",
      "ERROR! Colour not initialised");
    if (isInitPtr && rndmPtr->flat() < 0.5) return false;
    else return true;
  }

  // Mode 0: Purely random.
  if (inheritMode == 0) {
    if (rndmPtr->flat() < 0.5) return true;
    else return false;
  }

  // Safety checks: small, or approximately equal s01, s12.
  double a12 = abs(s01);
  double a23 = abs(s12);

  // Inverted mode (smallest invariant inherits - should only be used
  // for extreme variation checks).
  if (inheritMode < 0) {
    a12 = abs(s12);
    a23 = abs(s01);
    inheritMode = abs(inheritMode);
  }

  // Winner-takes-all mode.
  if (inheritMode == 2) {
    if (a12 > a23) return true;
    else return false;
  }
  double p12 = 0.5;
  if ( max(a12,a23) > TINY ) {
    if ( a12 < TINY ) p12 = 0.;
    else if ( a23 < TINY ) p12 = 1.;
    else {
      double r = a23/a12;
      if (r < TINY) p12 = 1. - r;
      else if (r > 1./TINY) p12 = 1./r;
      else p12 = 1./(1. + r);
    }
  }
  if (rndmPtr->flat() < p12) return true;
  else return false;

}

//==========================================================================

// The Resolution class.

//--------------------------------------------------------------------------

// Initialize.

bool Resolution::init() {

  // Check that pointers initialized.
  if (!isInitPtr) {
    printOut("Resolution::init","Cannot initialize, pointers not set.");
    return false;
  }

  // Set members.
  verbose          = settingsPtr->mode("Vincia:verbose");
  nFlavZeroMassSav = settingsPtr->mode("Vincia:nFlavZeroMass");
  isInit           = true;
  return isInit;

}

//--------------------------------------------------------------------------

// Sector resolution functions.

double Resolution::q2sector2to3(const Particle* a, const Particle* b,
  const Particle* j, bool) {

  // Construct basic 4-products.
  double saj  = 2*a->p()*j->p();
  double sjb  = 2*b->p()*j->p();
  double sab  = 2*a->p()*b->p();
  // Gluon emission.
  if (j->id() == 21) {
    // FF emission.
    if (a->isFinal() && b->isFinal()) {
      return saj * sjb / (saj + sjb + sab) ;

    // IF emission.
    } else if (b->isFinal()) {
      return saj * sjb / (saj + sab) ;

    // IF emission (a <-> b).
    } else if (a->isFinal()) {
      return saj * sjb / (sjb + sab) ;

    // II emission.
    } else {
      return saj * sjb / sab;
    }

  // FF Gluon splitting.
  } else if (a->isFinal() && b->isFinal()) {
    // Assume b is recoiler.
    double m2j  = pow2(j->m());
    double m2qq = saj + 2*m2j;
    // Find colour-connected invariant.
    double sX   = 0;
    if (j->col() != 0 && j->col() == b->acol())
      sX = sjb + m2j;
    else {
      sX = sab + m2j;
    }
    // Normalisation.
    double sAnt   = saj + sjb + sab + 2 * m2j;

    // Return the g->qq sector variable defined in arXiv:1109.3608.
    return m2qq * sqrt(sX/sAnt);
  } else {
    cout << "Sector criterion not implemented for II/IF splittings/conversions"
         << endl;
  }
  return -1;

}

//--------------------------------------------------------------------------

// Sector resolution function for 3->4 branchings (currently only used
// for gluon splitting, with m2qq as the measure).

double Resolution::q2sector3to4(const Particle*, const Particle* ,
  const Particle* j1, const Particle* j2) {
  Vec4   pqq  = j1->p() + j2->p();
  double m2qq = pqq.m2Calc();
  return m2qq;
}

//--------------------------------------------------------------------------

// Sector resolution function for 2->4 branchings (double emission).
// Assume j1 and j2 are colour connected, with a and b hard recoilers.

double Resolution::q2sector2to4(const Particle* a, const Particle* b,
  const Particle* j1, const Particle* j2) {
  return min(q2sector2to3(a,j2,j1),q2sector2to3(j1,b,j2));
}

//--------------------------------------------------------------------------

// Sector resolution function for 3->5 branchings
// (emission + splitting).

double Resolution::q2sector3to5(Particle* a, Particle* b,
  Particle* j1, Particle* j2, Particle* j3) {

  // j3 is gluon.
  Particle* gluPtr;
  Particle* qPtr;
  Particle* qBarPtr;
  if (j1->id() == 21) {
    gluPtr  = j1;
    qPtr    = (j2->id() > 0 ? j2 : j3);
    qBarPtr = (j2->id() < 0 ? j2 : j3);
  } else if (j2->id() == 21) {
    gluPtr  = j2;
    qPtr    = (j1->id() > 0 ? j1 : j3);
    qBarPtr = (j1->id() < 0 ? j1 : j3);
  } else if (j3->id() == 21) {
    gluPtr  = j3;
    qPtr    = (j2->id() > 0 ? j2 : j1);
    qBarPtr = (j2->id() < 0 ? j2 : j1);
  } else {
    cout << " q2sector3to5: unable to identify branching type" << endl;
    return 1.e19;
  }
  Vec4   pqq  = qPtr->p() + qBarPtr->p();
  double m2qq = pqq.m2Calc();
  Particle* colPtr = a;
  if (a->col()   != gluPtr->acol()) colPtr  = j1;
  if (j1->col()  != gluPtr->acol()) colPtr  = j2;
  if (j2->col()  != gluPtr->acol()) colPtr  = j3;
  if (j3->col()  != gluPtr->acol()) colPtr  = b;
  Particle* acolPtr = b;
  if (b->acol()  != gluPtr->col())  acolPtr = j3;
  if (j3->acol() != gluPtr->col())  acolPtr = j2;
  if (j2->acol() != gluPtr->col())  acolPtr = j1;
  if (j1->acol() != gluPtr->col())  acolPtr = a;
  double q2emit = q2sector2to3(colPtr,acolPtr,gluPtr);
  return min(q2emit,m2qq);

}

//--------------------------------------------------------------------------

// Sector accept function. Optionally prevent g->qq clusterings if
// that would reduce the number of fermion lines below some minimum
// (cheap way to indicate that Z->qq/ll always has at least one
// fermion pair).

double Resolution::findSector(vector<int>& iSec, vector<Particle> state,
  int nFmin) {

  map<int, int> colMap, acolMap;
  map<int, vector<int> > flavMap;
  vector<int> helX;
  int nFerm(0), nIn(0);
  for (int i=0; i<(int)state.size(); ++i) {
    colMap[state[i].isFinal() ? state[i].col() : state[i].acol()] = i;
    acolMap[state[i].isFinal() ? state[i].acol() : state[i].col()] = i;
    flavMap[state[i].isFinal() ? state[i].id() :
      -state[i].id()].push_back(i);
    helX.push_back(state[i].isFinal() ? state[i].pol() : -state[i].pol());
    if (state[i].isQuark() || state[i].isLepton()) ++nFerm;
    if (!state[i].isFinal()) ++nIn;
  }
  nFerm /= 2;

  // Initialise.
  double q2min = 1e19;
  iSec.clear();

  // Do all possible 2->3 clusterings. Note, needs modification for
  // 3->4 and 2->4 showers.
  for (int ir = 0; ir < (int)state.size(); ++ir) {
    // There are no sectors for emission into the initial state.
    if (!state[ir].isFinal()) continue;
    // If a gluon, compute LC pT from colour partners.
    if (state[ir].isGluon()) {
      int ia = colMap[state[ir].acol()];
      int ib = acolMap[state[ir].col()];
      double q2this = q2sector2to3(&state[ia],&state[ib],&state[ir]);
      if (q2this < q2min) {
        q2min = q2this;
        iSec.clear();
        iSec.push_back(ia);
        iSec.push_back(ib);
        iSec.push_back(ir);
      }
    }

    // If fermion, check if we are allowed to look at this clustering.
    if (state[ir].isQuark()) {
      // Skip if this clustering would reduce below number requested.
      if (nFerm <= nFmin) continue;
      // Cluster quark: recoiler is the anticolour partner.
      // Cluster antiquark: recoiler is colour partner.
      int ib = (state[ir].id() > 0) ? acolMap[state[ir].col()]
        : colMap[state[ir].acol()];

      // Loop over all same-flavour (anti)quarks (must in principle
      // also require opposite helicities).
      vector<int> aList;
      aList = flavMap[-state[ir].id()];
      for (int ifa = 0; ifa < (int)aList.size(); ++ifa) {
        int ia = aList[ifa];
        if (helX[ir] != helX[ia] || helX[ia] == 9 || helX[ir] == 9) {
          double q2this = q2sector2to3(&state[ia],&state[ib],&state[ir]);
          if (q2this < q2min) {
            q2min = q2this;
            iSec.clear();
            iSec.push_back(ia);
            iSec.push_back(ib);
            iSec.push_back(ir);
          }
        }
      }
    }
    // Other sectors (QED, EW, ...) to be implemented here.
  } // End loop over ir.

  // Return smallest scale found.
  return q2min;

}

//==========================================================================

// The VinciaCommon class.

//--------------------------------------------------------------------------

// Initialize the class.

bool VinciaCommon::init() {

  // Check initPtr.
  if (!isInitPtr) {
    if (verbose >= 1)
      printOut("VinciaCommon::init","Error! pointers not initialized");
    return false;
  }

  // Verbosity level and checks.
  verbose   = settingsPtr->mode("Vincia:verbose");
  epTolErr  = settingsPtr->parm("Check:epTolErr");
  epTolWarn = settingsPtr->parm("Check:epTolWarn");
  mTolErr   = settingsPtr->parm("Check:mTolErr");
  mTolWarn  = settingsPtr->parm("Check:mTolWarn");

  // Counters
  nUnkownPDG    = 0;
  nIncorrectCol = 0;
  nNAN          = 0;
  nVertex       = 0;
  nChargeCons   = 0;
  nMotDau       = 0;
  nUnmatchedMass.resize(2);
  nEPcons.resize(2);
  for (int i=0; i<2; i++) {
    nUnmatchedMass[i] = 0;
    nEPcons[i]        = 0;
  }

  // Quark masses
  mt                 = particleDataPtr->m0(6);
  if (mt < TINY) mt  = 171.0;
  mb                 = min(mt,particleDataPtr->m0(5));
  if (mb < TINY) mb  = min(mt,4.8);
  mc                 = min(mb,particleDataPtr->m0(4));
  if (mc < TINY) mc  = min(mb,1.5);
  ms                 = min(mc,particleDataPtr->m0(3));
  if (ms < TINY) ms  = min(mc,0.1);

  // Number of flavours to treat as massless in clustering and
  // kinematics maps.
  nFlavZeroMass = settingsPtr->mode("Vincia:nFlavZeroMass");

  // Default alphaS, with and without CMW.
  double alphaSvalue = settingsPtr->parmDefault("Vincia:alphaSvalue");
  int    alphaSorder = settingsPtr->modeDefault("Vincia:alphaSorder");
  int    alphaSnfmax = settingsPtr->modeDefault("Vincia:alphaSnfmax");
  bool   useCMW      = settingsPtr->flagDefault("Vincia:useCMW");
  alphaStrongDef.init(    alphaSvalue, alphaSorder, alphaSnfmax, false);
  alphaStrongDefCMW.init( alphaSvalue, alphaSorder, alphaSnfmax, true);

  // Strong coupling for use in merging.
  alphaSvalue  = settingsPtr->parm("Vincia:alphaSvalue");
  alphaSorder  = settingsPtr->mode("Vincia:alphaSorder");
  alphaSnfmax  = settingsPtr->mode("Vincia:alphaSnfmax");
  useCMW       = settingsPtr->flag("Vincia:useCMW");
  alphaS.init(alphaSvalue, alphaSorder, alphaSnfmax, useCMW);

  // User alphaS, with and without CMW.
  alphaSvalue  = settingsPtr->parm("Vincia:alphaSvalue");
  alphaSorder  = settingsPtr->mode("Vincia:alphaSorder");
  alphaSnfmax  = settingsPtr->mode("Vincia:alphaSnfmax");
  useCMW       = settingsPtr->flag("Vincia:useCMW");
  alphaStrong.init(    alphaSvalue, alphaSorder, alphaSnfmax, false);
  alphaStrongCMW.init( alphaSvalue, alphaSorder, alphaSnfmax, true);

  // Freeze and minimum scales.
  mu2freeze    = pow2(settingsPtr->parm("Vincia:alphaSmuFreeze"));
  alphaSmax    = settingsPtr->parm("Vincia:alphaSmax");

  // Find the overall minimum scale. Take into account the freezeout
  // scale, Lambda pole, and alphaSmax.
  double muMin = max(sqrt(mu2freeze),1.05*alphaS.Lambda3());
  double muMinASmax;
  if (alphaStrong.alphaS(mu2min) < alphaSmax) {
    muMinASmax = muMin;
  } else if (settingsPtr->mode("Vincia:alphaSorder") == 0) {
    muMinASmax = muMin;
  } else {
    muMinASmax = muMin;
    while (true) {
      if (alphaS.alphaS(pow2(muMinASmax)) < alphaSmax) break;
      muMinASmax += 0.001;
    }
  }
  mu2min = pow2( max(muMinASmax, muMin) );

  // EM coupling for use in merging. Dummy, as no EW clusterings.
  alphaEM.init(1, settingsPtr);

  // Return.
  isInit = true;
  return true;

}

//--------------------------------------------------------------------------

// Function to find the lowest meson mass for a parton pair treating
// gluons as down quarks. Used to determine hadronisation boundaries
// consistent with assumption of string length > 1 meson.

double VinciaCommon::mHadMin(const int id1in, const int id2in) {
  int id1 = abs(id1in);
  if (id1 == 21 || id1 <= 2) id1 = 1;
  int id2 = abs(id2in);
  if (id2 == 21 || id1 <= 2) id2 = 1;
  int idMes = max(id1,id2)*100 + min(id1,id2)*10 + 1;

  // Special for ssbar, use eta rather than eta'.
  if (idMes == 331) idMes = 221;
  return particleDataPtr->m0(idMes);

}

//--------------------------------------------------------------------------

// Function to check the event after each branching, mostly copied
// from Pythia8.

bool VinciaCommon::showerChecks(Event& event, bool ISR) {

  // Only for verbose >= 4.
  if (verbose < 4) return true;

  // Count incoming partons (beam daughters) with negative momentum and charge.
  Vec4 pSum;
  double chargeSum = 0.0;
  bool beamRemnantAdded = false;
  for (int i = 0; i < event.size(); ++i){
    if(!beamRemnantAdded){
      if (!event[i].isFinal()) {
        if ( (event[i].mother1() == 1) || (event[i].mother1() == 2) ) {
          pSum      -= event[i].p();
          chargeSum -= event[i].charge();
        }
      }
      if(abs(event[i].status())==63){
        beamRemnantAdded=true;
        pSum = -(event[1].p() + event[2].p());
        chargeSum = -(event[1].charge()+event[2].charge());
      }
    }
  }
  double eLab = abs(pSum.e());

  // Loop over particles in the event.
  for (int i = 0; i < event.size(); ++i) {

    // Look for any unrecognized particle codes.
    int id = event[i].id();
    if (id == 0 || !particleDataPtr->isParticle(id)) {
      nUnkownPDG++;
      if (nUnkownPDG == 1){
        cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
             << ": unknown particle code"
             << ", i = " << i << ", id = " << id << endl;
        return false;
      }
    }

    // Check that colour assignments are the expected ones.
    else {
      int colType = event[i].colType();
      int col     = event[i].col();
      int acol    = event[i].acol();
      if (    (colType ==  0 && (col  > 0 || acol  > 0))
        || (colType ==  1 && (col <= 0 || acol  > 0))
        || (colType == -1 && (col  > 0 || acol <= 0))
        || (colType ==  2 && (col <= 0 || acol <= 0)) ) {
        nIncorrectCol++;
        if (nIncorrectCol == 1){
          cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
               << ": incorrect colours"
               << ", i = " << i << ", id = " << id << " cols = " << col
               << " " << acol << endl;
          return false;
        }
      }
    }

    // Look for particles with mismatched or not-a-number energy/momentum/mass.
    if (abs(event[i].px()) >= 0.0 && abs(event[i].py()) >= 0.0
        && abs(event[i].pz()) >= 0.0 && abs(event[i].e())  >= 0.0
        && abs(event[i].m())  >= 0.0 ) {
      double errMass  = abs(event[i].mCalc() - event[i].m()) /
        max( 1.0, event[i].e());

      if (errMass > mTolErr) {
        nUnmatchedMass[0]++;
        if (nUnmatchedMass[0] == 1){
          cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
               << ": unmatched particle energy/momentum/mass"
               << ", i = " << i << ", id = " << id << endl;
          return false;
        }
      } else if (errMass > mTolWarn) {
        nUnmatchedMass[1]++;
        if (nUnmatchedMass[1] == 1){
          cout << "WARNING in Vincia::ShowerChecks"
               << (ISR ? "(ISR)" : "(FSR)")
               << ": not quite matched particle energy/momentum/mass"
               << ", i = " << i << ", id = " << id << endl;
        }
      }
    } else {
      nNAN++;
      if (nNAN == 1){
        cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
             << ": not-a-number energy/momentum/mass"
             << ", i = " << i << ", id = " << id << endl;
        return false;
      }
    }

    // Look for particles with not-a-number vertex/lifetime.
    if (abs(event[i].xProd()) >= 0.0 && abs(event[i].yProd()) >= 0.0
        && abs(event[i].zProd()) >= 0.0 && abs(event[i].tProd()) >= 0.0
        && abs(event[i].tau())   >= 0.0) {
    } else {
      nVertex++;
      if (nVertex == 1){
        cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
             << ": not-a-number vertex/lifetime"
             << ", i = " << i << ", id = " << id << endl;
        return false;
      }
    }

    // Add final-state four-momentum and charge.
    if (event[i].isFinal()) {
      pSum      += event[i].p();
      chargeSum += event[i].charge();
    }
  } // End of particle loop.

  // Check energy-momentum/charge conservation.
  double epDev = abs( pSum.e()) + abs(pSum.px()) + abs(pSum.py())
    + abs(pSum.pz() );
  if (epDev > epTolErr * eLab) {
    nEPcons[0]++;
    if (nEPcons[0] == 1) {
      cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
           << ": energy-momentum not conserved" << endl;
      cout << " epDev = " << epDev << " Ein = " << eLab
           << " pSum = " << pSum << endl;
      return false;
    }
  } else if (epDev > epTolWarn * eLab) {
    nEPcons[1]++;
    if (nEPcons[1] == 1)
      cout << "WARNING in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
           << ": energy-momentum not quite conserved" << endl;
  }
  if (abs(chargeSum) > 0.1) {
    nChargeCons++;
    if (nChargeCons == 1){
      cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
           << ": charge not conserved" << endl;
      return false;
    }
  }

  // Check that mother and daughter information match for each particle.
  vector<int> noMot, noDau;
  vector< pair<int,int> > noMotDau;

  // Loop through the event and check that there are beam particles.
  bool hasBeams = false;
  for (int i = 0; i < event.size(); ++i) {
    int status = event[i].status();
    if (abs(status) == 12) hasBeams = true;

    // Check that mother and daughter lists not empty where not expected to.
    vector<int> mList = event[i].motherList();
    vector<int> dList = event[i].daughterList();
    if (mList.size() == 0 && abs(status) != 11 && abs(status) != 12)
      noMot.push_back(i);
    if (dList.size() == 0 && status < 0 && status != -11)
      noDau.push_back(i);

    // Check that the particle appears in the daughters list of each mother.
    for (int j = 0; j < int(mList.size()); ++j) {
      if (event[mList[j]].daughter1() <= i
          && event[mList[j]].daughter2() >= i) continue;
      vector<int> dmList = event[mList[j]].daughterList();
      bool foundMatch = false;
      for (int k = 0; k < int(dmList.size()); ++k)
        if (dmList[k] == i) {
          foundMatch = true;
          break;
        }
      if (!hasBeams && mList.size() == 1 && mList[0] == 0) foundMatch = true;
      if (!foundMatch) {
        bool oldPair = false;
        for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == mList[j] && noMotDau[k].second == i) {
            oldPair = true;
            break;
          }
        if (!oldPair) noMotDau.push_back( make_pair( mList[j], i) );
      }
    }

    // Check that the particle appears in the mothers list of each daughter.
    for (int j = 0; j < int(dList.size()); ++j) {
      if (event[dList[j]].statusAbs() > 80
          && event[dList[j]].statusAbs() < 90
          && event[dList[j]].mother1() <= i
          && event[dList[j]].mother2() >= i ) continue;
      vector<int> mdList = event[dList[j]].motherList();
      bool foundMatch = false;
      for (int k = 0; k < int(mdList.size()); ++k)
        if (mdList[k] == i) {
          foundMatch = true;
          break;
        }
      if (!foundMatch) {
        bool oldPair = false;
        for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == i && noMotDau[k].second == dList[j]) {
            oldPair = true;
            break;
          }
        if (!oldPair) noMotDau.push_back( make_pair( i, dList[j]) );
      }
    } // End loop through the event.
  }

  // Warn if any errors were found.
  if (noMot.size() > 0 || noDau.size() > 0 || noMotDau.size() > 0) {
    nMotDau++;
    if (nMotDau == 1) {
      cout << "ERROR in Vincia::ShowerChecks" << (ISR ? "(ISR)" : "(FSR)")
           << ": mismatch in daughter and mother lists" << endl;
      // Print some more info.
      if (noMot.size() > 0) {
        cout << " missing mothers for particles ";
        for (int i = 0; i < int(noMot.size()); ++i) cout << noMot[i] << " ";
        cout << endl;
      }
      if (noDau.size() > 0) {
        cout << " missing daughters for particles ";
        for (int i = 0; i < int(noDau.size()); ++i) cout << noDau[i] << " ";
        cout << endl;
      }
      if (noMotDau.size() > 0) {
        cout << " inconsistent history for (mother,daughter) pairs ";
        for (int i = 0; i < int(noMotDau.size()); ++i)
          cout << "(" << noMotDau[i].first << ","
               << noMotDau[i].second << ") ";
        cout << endl;
      }
      return false;
    }
  }

  //Made it to here: no major problems.
  return true;

}

//--------------------------------------------------------------------------

// Get the shower starting scale.

double VinciaCommon::getShowerStartingScale(int iSys,
  PartonSystems* partonSystemsPtr, const Event& event, double sbbSav) {

  // Depending on user choice shower starts at q2maxFudge *
  // factorization scale of phase space maximum.
  int    qMaxMatch    = settingsPtr->mode("Vincia:QmaxMatch");
  double q2maxFudge   = pow2(settingsPtr->parm("Vincia:QmaxFudge"));
  bool   hasFSpartons = false;
  int nOut = partonSystemsPtr->sizeOut(iSys);
  vector<int> iFS;
  for (int iOut = 0; iOut < nOut; ++iOut) {
    int i = partonSystemsPtr->getOut(iSys,iOut);
    int idAbs = event[i].idAbs();
    if ((idAbs >= 1 && idAbs <= 5) || idAbs == 21) hasFSpartons = true;
    iFS.push_back(i);
  }
  if (qMaxMatch == 1 || (qMaxMatch == 0 && hasFSpartons) ) {
    double Q2facSav = sbbSav;
    double Q2facFix  = settingsPtr->parm("SigmaProcess:factorFixScale");
    double Q2facMult = settingsPtr->parm("SigmaProcess:factorMultFac");

    // Ask Pythia about 2 -> 1 scale.
    if (nOut == 1) {
      int factorScale1 = settingsPtr->mode("SigmaProcess:factorScale1");
      Q2facSav = ( factorScale1 == 1 ? Q2facMult*event[iFS[0]].m2Calc() :
        Q2facFix );

    // Ask Pythia about 2 -> 2 scale.
    } else if (iFS.size() == 2) {
      int factorScale2 = settingsPtr->mode("SigmaProcess:factorScale2");
      double mT21 = event[iFS[0]].mT2(), mT22 = event[iFS[1]].mT2();
      double sHat = m2(event[iFS[0]],event[iFS[1]]);
      double tHat = m2(event[3],event[iFS[0]]);
      if (factorScale2 == 1) Q2facSav = min(mT21,mT22);
      else if (factorScale2 == 2) Q2facSav = sqrt(mT21*mT22);
      else if (factorScale2 == 3) Q2facSav = 0.5*(mT21+mT22);
      else if (factorScale2 == 4) Q2facSav = sHat;
      else if (factorScale2 == 5) Q2facSav = Q2facFix;
      else if (factorScale2 == 6) Q2facSav = abs(-tHat);
      if (factorScale2 != 5) Q2facSav *= Q2facMult;

    // Ask Pythia about 2 -> 3 scale.
    } else if (iFS.size() == 3) {
      int factorScale3 = settingsPtr->mode("SigmaProcess:factorScale3");
      double mT21 = event[iFS[0]].mT2(), mT22 = event[iFS[1]].mT2(),
        mT23 = event[iFS[2]].mT2();
      double mT2min1 = min(mT21,min(mT22,mT23));
      double mT2min2 = max(max(min(mT21,mT22),min(mT21,mT23)),min(mT22,mT23));
      double sHat    = m2(event[iFS[0]],event[iFS[1]],event[iFS[2]]);
      if (factorScale3 == 1) Q2facSav = mT2min1;
      else if (factorScale3 == 2) Q2facSav = sqrt(mT2min1*mT2min2);
      else if (factorScale3 == 3) Q2facSav = pow(mT21*mT22*mT23,1.0/3.0);
      else if (factorScale3 == 4) Q2facSav = (mT21+mT22+mT23)/3.0;
      else if (factorScale3 == 5) Q2facSav = sHat;
      else if (factorScale3 == 6) Q2facSav = Q2facFix;
      if (factorScale3 != 6) Q2facSav *= Q2facMult;

    // Unknown, leave as is, all emissions allowed now.
    } else {;}
    return q2maxFudge*Q2facSav;
  }
  return sbbSav;

}

//--------------------------------------------------------------------------

// FF clustering map(s) for massless partons. Inputs are as follows:
//   kMapType = map number ( 1 : ARIADNE, 2 : PS, 3 : Kosower)
//   pIn      = Vec4 list (arbitrary length, but greater than 3)
//   a,r,b    = indices of 3 particles in pIn to be clustered
// Outputs are as follows:
//   pClu     = pIn but with the a and b momenta replaced by the clustered
//              aHat and bHat and with r erased.
// For example:
//   pIn(...,a,...,r-1,r,r+1,...,b,...) ->
//   pOut(...,aHat,...,r-1,r+1,...,bHat,...)
// with {a,r,b} -> {aHat,bHat} using kinematics map kMapType.

bool VinciaCommon::map3to2FFmassless(vector<Vec4>& pClu, vector<Vec4> pIn,
  int kMapType, int a, int r, int b) {

  // Intialize and sanity check.
  pClu=pIn;
  if (max(max(a,r),b) > int(pIn.size()) || min(min(a,r),b) < 0) {
    if (verbose >= 3)
      printOut("VinciaCommon::map3to2FFmassless",
        "Error! Unable to cluster (a,r,b) = "+
        num2str(a)+num2str(r)+num2str(b)+" p.size ="
        +num2str(int(pIn.size())));
    return false;
  }

  // Verbose output.
  if (verbose >= superdebug) {
    printOut("VinciaCommon:map3to2FFmassless", "called with ");
    cout << "pi = " << pIn[a];
    cout << "pj = " << pIn[r];
    cout << "pk = " << pIn[b];
  }

  // Compute total invariant mass squared.
  Vec4 pSum    = pIn[a] + pIn[r] + pIn[b];
  double m2Ant = pSum.m2Calc();
  if (m2Ant < 1e-20) {
    printOut("VinciaCommon::map3to2FFmassless",
             "Massless or spacelike system. Cannot find rest frame");
    return false;
  }

  // ARIADNE and PS maps (recoded from old F77 maps for v.1.2.01)
  // (maps -1 and -2 are special: force A or B to take all recoil).
  if (kMapType == 1 || kMapType == 2 || kMapType == -1 || kMapType == -2) {

    // Make copies of PA, PB, and compute sum of LAB momenta and CM energy.
    Vec4 paDum   = pIn[a];
    Vec4 pbDum   = pIn[b];
    double eCM = sqrt(m2Ant);
    paDum.bstback(pSum);
    pbDum.bstback(pSum);

    // Rotate so a goes into upper half of (x,z) plane.
    double phiA = paDum.phi();
    paDum.rot(0.,-phiA);
    pbDum.rot(0.,-phiA);

    // Rotate so a goes onto z axis.
    double theta = paDum.theta();
    pbDum.rot(-theta,0.);

    // Rotate so (r,b) go into (x,z) plane.
    double phiB = pbDum.phi();

    // Compute psi(a,ahat) from A, B energies and theta(AB).
    double thetaAB = pbDum.theta();
    double psi = 0.0;

    // ARIADNE angle (smoothly varying antenna recoil).
    if (kMapType == 1) {
      psi = pbDum.e()*pbDum.e()/(paDum.e()*paDum.e() + pbDum.e()*pbDum.e())
        * (M_PI - thetaAB);

    // PS angle (direction a fixed if sar > srb, and vice versa). Use
    // org momenta, since new ones may not all have been rotated.
    } else if (kMapType == 2) {
      Vec4 pAR = pIn[a] + pIn[r];
      Vec4 pRB = pIn[r] + pIn[b];
      double sar = pAR.m2Calc();
      double srb = pRB.m2Calc();
      if (sar > srb) psi = 0.0;
      else psi = (M_PI - thetaAB);

    // Force A to be the emitter. B recoils longitudinally.
    } else if (kMapType == -1) {
      psi = M_PI - thetaAB;

    // Force B to be the emitter. A recoils longitudinally.
    } else if (kMapType == -2) {
      psi = 0.0;
    }
    // Now we know everything:
    // CM -> LAB : -PSI, PHIB, THETA, PHIA, BOOST

    // Set up initial massless AHAT, BHAT with AHAT along z axis.
    pClu[a] = Vec4(0.,0.,eCM/2.,eCM/2.);
    pClu[b] = Vec4(0.,0.,-eCM/2.,eCM/2.);

    // 1) Take into account antenna recoil, and rotate back in phiB.
    pClu[a].rot(-psi,phiB);
    pClu[b].rot(-psi,phiB);

    // 2) Rotate back in theta and phiA.
    pClu[a].rot(theta,phiA);
    pClu[b].rot(theta,phiA);

    // 3) Boost back to LAB.
    pClu[a].bst(pSum);
    pClu[b].bst(pSum);

  // kMapType = 3, 4 (and catchall for undefined kMapType
  // values). Implementation of the massless Kosower antenna map(s).
  } else {

    // Compute invariants.
    double s01  = 2*pIn[a]*pIn[r];
    double s12  = 2*pIn[r]*pIn[b];
    double s02  = 2*pIn[a]*pIn[b];

    // Check whether the arguments need to be reversed for kMapType == 4.
    if (kMapType == 4 && ( ! (s01 < s12) ) ) {
      if (verbose >= superdebug) {
        printOut("VinciaCommon::map3to2FFmassless",
          "choose parton i as the recoiler");
      }
      // Call the function with reverse arguments, then return.
      return map3to2FFmassless(pClu, pIn, kMapType, b, r, a);
    }
    double sAnt = s01 + s12 + s02;

    // Compute coefficients.
    //  kMapType  = 3: GGG choice
    //  kMapType >= 4: r=1
    double rMap = kMapType == 3 ? s12/(s01 + s12) : 1;
    double rho  = sqrt(1.0+(4*rMap*(1.0-rMap)*s01*s12)/sAnt/s02);
    double x    = 0.5/(s01+s02)*((1+rho)*(s01+s02)+(1+rho-2*rMap)*s12);
    double z    = 0.5/(s12+s02)*((1-rho)*sAnt-2*rMap*s01);

    // Compute reclustered momenta.
    pClu[a]     =     x*pIn[a] +     rMap*pIn[r] +     z*pIn[b];
    pClu[b]     = (1-x)*pIn[a] + (1-rMap)*pIn[r] + (1-z)*pIn[b];
  }

  // A dimensionless quantitiy to compare with TINY.
  if (pClu[a].m2Calc()/m2Ant >= TINY || pClu[b].m2Calc()/m2Ant >= TINY) {
    if (verbose >= 3)
      printOut("VinciaCommon::map3to2FFmassless",
               "on-shell check failed. m2I/sIK ="
               + num2str(pClu[a].m2Calc()/m2Ant)+" m2K/m2Ant ="
               + num2str(pClu[b].m2Calc()/m2Ant)+" m2Ant = "+num2str(m2Ant));
    return false;
  }

  // Erase the clustered momentum and return.
  pClu.erase(pClu.begin()+r);
  return true;

}

//--------------------------------------------------------------------------

// Implementations of FF clustering maps for massive partons. See
// VinciaCommon::map3to2FFmassless for details but with the additional
// input of masses mI and mK for the first and second parent
// particles.

bool VinciaCommon::map3to2FFmassive(vector<Vec4>& pClu, vector<Vec4> pIn,
  int kMapType, int a, int r, int b, double mI, double mK) {

  // If both parent masses are negligible and all the momenta are
  // measure off-shellness normalised to average energy of the partons
  // to be clustered, avoids small denominator for collinear massless
  // p_a, p_r, p_b.
  double eAv = 1.0/3.0*( pIn[a].e() + pIn[r].e() + pIn[b].e() );
  if (mI/eAv < TINY && mK/eAv < TINY && pIn[a].mCalc()/eAv < TINY
    && pIn[r].mCalc()/eAv < TINY && pIn[b].mCalc()/eAv < TINY ) {
    return map3to2FFmassless(pClu, pIn, kMapType, a, r, b);
  }

  // Intialize and sanity check.
  pClu = pIn;
  if (max(max(a,r),b) > int(pIn.size()) || min(min(a,r),b) < 0) return false;

  // Verbose output.
  if (verbose >= superdebug) {
    printOut("VinciaCommon:map3to2FFmassive", "called with ");
    cout << "p0 = " << pIn[a];
    cout << "p1 = " << pIn[r];
    cout << "p2 = " << pIn[b];
  }

  // ARIADNE map not defined for massive partons; use Kosower map instead.
  if (kMapType == 1) kMapType = 3;

  // Longitudinal map; use Kosower map with r = 1.
  if (kMapType == 2) kMapType = 4;

  // Forced-longitudinal maps.
  if (kMapType < 0) {
    printOut("VinciaCommon::map3to2FFmassive", "longitudinal clustering maps "
             "not yet implemented for massive partons.");
    return false;
  }

  // Compute invariants.
  double m0   = pIn[a].mCalc();
  double m1   = pIn[r].mCalc();
  double m2   = pIn[b].mCalc();
  double s01  = 2*pIn[a]*pIn[r];
  double s12  = 2*pIn[r]*pIn[b];
  double s02  = 2*pIn[a]*pIn[b];

  // Check whether the arguments need to be reversed for mapType == 4.
  if (kMapType == 4 && (! (s01 < s12) ) ) {
    return map3to2FFmassive(pClu, pIn, kMapType, b, r, a, mK, mI);
  }
  Vec4   pAnt    = pIn[a] + pIn[r] + pIn[b];
  double m2Ant   = pAnt.m2Calc();
  double mAnt    = sqrt(m2Ant);

  // Rewrite the determination in terms of dimensionless variables.
  // Note normalisation choice here is mAnt, rather than sAnt.
  double mu0 = m0/mAnt;
  double mu1 = m1/mAnt;
  double mu2 = m2/mAnt;
  double y01 = s01/m2Ant;
  double y12 = s12/m2Ant;
  double y02 = s02/m2Ant;
  double y01min = 2*mu0*mu1;
  double y12min = 2*mu1*mu2;
  double y02min = 2*mu0*mu2;
  double muI = mI/mAnt;
  double muK = mK/mAnt;
  double yIK = 1. - pow2(muI) - pow2(muK);
  double yIKmin = 2*muI*muK;
  double sig2 = 1 + pow2(muI) - pow2(muK);
  double gdetdimless = gramDet(y01,y12,y02,mu0,mu1,mu2);
  double gdetdimless01 = (y02*y12-2.*pow2(mu2)*y01)/4.;
  double gdetdimless12 = (y02*y01-2.*pow2(mu0)*y12)/4.;
  double rMap = 1.;
  if ( kMapType == 3) {
    rMap =
      (
        sig2
        + sqrt( pow2(yIK) - pow2(yIKmin) )
        *( y12 - y12min - ( y01 - y01min ) )
        /( y12 - y12min + ( y01 - y01min ) )
       )/2.;

  // Fallback: map with massless r -> 1.
  } else  {
    double mu01squa = pow2(mu0) + pow2(mu1) + y01;
    double lambda = 1 + pow2(mu01squa) + pow4(mu2) - 2*mu01squa - 2*pow2(mu2)
      - 2*mu01squa*pow2(mu2);
    rMap = (sig2 + sqrt(pow2(yIK) - pow2(yIKmin))
            *(1 - pow2(mu0) - pow2(mu1) + pow2(mu2) - y01)/sqrt(lambda))/2.;
  }

  // Compute reclustered momenta.
  double bigsqr = sqrt(
    16.*( rMap*(1.-rMap) - (1.-rMap)*pow2(muI) - rMap*pow2(muK) )*gdetdimless
    + ( pow2(y02) - pow2(y02min) )*( pow2(yIK) - pow2(yIKmin) ));
  double x = (
    sig2*( pow2(y02) - pow2(y02min) + 4.*gdetdimless01)
    + 8.*rMap*( gdetdimless - gdetdimless01 )
    + bigsqr*( 1. - pow2(mu0) - pow2(mu1) + pow2(mu2) - y01)
              )/(2.*( 4.*gdetdimless + pow2(y02) - pow2(y02min) ));
  double z = (
    sig2*( pow2(y02) - pow2(y02min) + 4.*gdetdimless12)
    + 8.*rMap*( gdetdimless - gdetdimless12 )
    - bigsqr*( 1. + pow2(mu0) - pow2(mu1) - pow2(mu2) - y12)
              )/(2.*( 4.*gdetdimless + pow2(y02) - pow2(y02min)));
  pClu[a] =     x*pIn[a] +     rMap*pIn[r] +     z*pIn[b];
  pClu[b] = (1-x)*pIn[a] + (1-rMap)*pIn[r] + (1-z)*pIn[b];

  // Check if on-shell.
  double offshellnessI = abs(pClu[a].m2Calc() - pow2(mI))/m2Ant;
  double offshellnessK = abs(pClu[b].m2Calc() - pow2(mK))/m2Ant;
  if (offshellnessI > TINY || offshellnessK > TINY) {
    if (verbose >= 3) {
      printOut("VinciaCommon::map3to2FFmassive","on-shell check failed");
    }
    return false;
  }

  // Erase the clustered parton and return.
  pClu.erase(pClu.begin()+r);
  return true;

}

//--------------------------------------------------------------------------

// Implementations of IF clustering maps for massive partons.

bool VinciaCommon::map3to2IFmassive(vector<Vec4>& pClu, vector<Vec4>& pIn,
  double saj, double sjk, double sak) {
  double sAK = saj + sak - sjk;
  Vec4 pA = pIn[0];
  pA.rescale4(sAK/(sAK + sjk));
  Vec4 pK = pA - pIn[0] + pIn[1] + pIn[2];
  pClu.push_back(pA);
  pClu.push_back(pK);
  return true;

}

//--------------------------------------------------------------------------

// Implementations of II clustering maps for massive partons.

bool VinciaCommon::map3to2IImassive(vector<Vec4>& pClu, vector<Vec4>& pIn,
  vector<Vec4>& pRec, double saj, double sjb, double sab, bool doBoost) {

  // Scale factors and momenta.
  double sAB = sab - saj - sjb;
  double rescaleFacA = 1./sqrt(sab/sAB * (sab - saj)/(sab - sjb));
  double rescaleFacB = 1./sqrt(sab/sAB * (sab - sjb)/(sab - saj));
  Vec4 pA = pIn[0];
  pA.rescale4(rescaleFacA);
  pClu.push_back(pA);
  Vec4 pB = pIn[2];
  pB.rescale4(rescaleFacB);
  pClu.push_back(pB);
  Vec4 pInSum = pIn[0] + pIn[2] - pIn[1];
  Vec4 pCluSum = pClu[0] + pClu[1];

  // Perform boost - if doBoost, we boost back to the lab frame.
  if (doBoost) {
    for (int i=0; i<(int)pRec.size(); i++) {
      pRec[i].bstback(pInSum);
      pRec[i].bst(pCluSum);
    }

  // Otherwise stay in the current frame. Adjust clustered momenta.
  } else {
    for (int i=0; i<(int)pClu.size(); i++) {
      pClu[i].bstback(pCluSum);
      pClu[i].bst(pInSum);
    }
  }
  return true;

}

//--------------------------------------------------------------------------

// 2->3 branching kinematics: output arguments first, then input
// arguments.  The output is p[i,j,k] whil the inputs are p[I,K],
// kMapType, inviariants[sIK,sij,sjk], phi, and
// masses[mi,mj,mk]. Note, sab defined as 2*pa.pb.

bool VinciaCommon::map2to3FFmassive(vector<Vec4>& pThree,
  const vector<Vec4>& pTwo, int kMapType, const vector<double>& invariants,
  double phi, vector<double> masses) {

  // Hand off to massless map if all masses << sIK.
  if (masses.size() < 3 ||
    (masses[0] <= TINY && masses[1] <= TINY && masses[2] <= TINY))
    return map2to3FFmassless(pThree,pTwo,kMapType,invariants,phi);

  // Antenna invariant mass and sIK = 2*pI.pK.
  double m2Ant = m2(pTwo[0],pTwo[1]);
  double mAnt  = sqrt(m2Ant);
  double sAnt  = invariants[0];
  if (sAnt <= 0.0) return false;

  // Masses and invariants
  double mass0 = max(0.,masses[0]);
  double mass1 = max(0.,masses[1]);
  double mass2 = max(0.,masses[2]);
  // Check for totally closed phase space. Should normally have
  // happened before generation of invariants but put last-resort
  // check here since not caught by Gram determinant.
  if (mAnt < mass0 + mass1 + mass2) {
    cout <<" (VinciaCommon::map2to3FFmassive:) "
         <<"ERROR! unphysical phase space\n";
    return false;
  }
  double s01 = max(0.,invariants[1]);
  double s12 = max(0.,invariants[2]);
  double s02 = m2Ant - s01 - s12 - pow2(mass0) - pow2(mass1) - pow2(mass2);
  if (s02 <= 0.) return false;

  // Check whether we are inside massive phase space.
  double gDet = gramDet(s01, s12, s02, mass0, mass1, mass2);
  if (gDet <= 0.) {
    if (verbose >= 9)
      cout << "   map2to3FFmassive: failed massive phase space check" << endl;
    return false;
  }

  // Verbose output.
  if (verbose >= 7) {
    cout << " (VinciaCommon::map2to3FFmassive:) m  = " << num2str(mAnt)
         << " sqrtsIK = " << num2str(sqrt(sAnt)) << "   sqrts(ij,jk,ik) ="
         << num2str(sqrt(s01)) << " " << num2str(sqrt(s12)) << " "
         << num2str(sqrt(s02)) << "   m(i,j,k) =" << num2str(mass0) << " "
         << num2str(mass1) << " " << num2str(mass2) << " D = " << gDet << endl;
    RotBstMatrix M;
    Vec4 p1cm = pTwo[0];
    Vec4 p2cm = pTwo[1];
    M.toCMframe(p1cm,p2cm);
    p1cm.rotbst(M);
    p2cm.rotbst(M);
    Vec4 tot = p1cm+p2cm;
    cout << " (VinciaCommon::map2to3FFmassive:) starting dipole in CM\n"
         << " p1cm = " << p1cm << " p2cm = " << p2cm
         << " total = " << tot << endl;
  }

  // Set up kinematics in rest frame.
  double E0 = (pow2(mass0) + s01/2 + s02/2)/mAnt;
  double E1 = (pow2(mass1) + s12/2 + s01/2)/mAnt;
  double E2 = (pow2(mass2) + s02/2 + s12/2)/mAnt;

  // Make sure energies > masses (should normally be ensured by
  // combination of open phase space and positive Gram determinant).
  if (E0 < mass0 || E1 < mass1 || E2 < mass2) {
    cout <<" (VinciaCommon::map2to3FFmassive:) ERROR! "
         <<"Unphysical energy value(s)\n";
    return false;
  }
  double ap0 = sqrt( pow2(E0) - pow2(mass0) );
  double ap1 = sqrt( pow2(E1) - pow2(mass1) );
  double ap2 = sqrt( pow2(E2) - pow2(mass2) );
  double cos01 = (E0*E1 - s01/2)/(ap0*ap1);
  double cos02 = (E0*E2 - s02/2)/(ap0*ap2);

  // Protection: num. precision loss for small (ultracollinear) invariants.
  if ( 1-abs(cos01) < 1e-15 ) cos01 = cos01 > 0 ? 1. : -1.;
  if ( 1-abs(cos02) < 1e-15 ) cos02 = cos02 > 0 ? 1. : -1.;

  // Use positive square root for sine.
  double sin01 = (abs(cos01) < 1) ? sqrt(abs(1.0 - pow2(cos01))) : 0.0;
  double sin02 = (abs(cos02) < 1) ? sqrt(abs(1.0 - pow2(cos02))) : 0.0;

  // Set momenta in CMz frame (frame with 1 oriented along positive z
  // axis and event in (x,z) plane).
  Vec4 p1(0.0,0.0,ap0,E0);
  Vec4 p2(-ap1*sin01,0.0,ap1*cos01,E1);
  Vec4 p3(ap2*sin02,0.0,ap2*cos02,E2);

  // Verbose output.
  if (verbose >= superdebug) {
    Vec4 tot = p1+p2+p3;
    cout  << " (map2to3FFmassive:) configuration in CM* (def: 1 along z)\n";
    cout  << " k1* =  " << p1 << " k2* =  " << p2 << " k3* =  " << p3
          << " total = " << tot << endl;
  }

  // Choose global rotation around axis perpendicular to event plane.
  double psi;

  // kMapType = -2(-1): force A(B) to be purely longitudinal recoiler.
  if (kMapType == -2) psi = 0.0;
  else if (kMapType == -1) psi = M_PI - acos(cos02);
  // Else general antenna-like recoils.
  else {
    double m2I = max(0.0,m2(pTwo[0]));
    double m2K = max(0.0,m2(pTwo[1]));
    double sig2 = m2Ant + m2I - m2K;
    double sAntMin = 2*sqrt(m2I*m2K);
    double s01min = max(0.0,2*mass0*mass1);
    double s12min = max(0.0,2*mass1*mass2);
    double s02min = max(0.0,2*mass0*mass2);

    // The r and R parameters in arXiv:1108.6172.
    double rAntMap = ( sig2 + sqrt( pow2(sAnt) - pow2(sAntMin) )
      * ( s12-s12min - (s01-s01min) )
      / ( s01-s01min + s12-s12min ) ) / (2*m2Ant);
    double bigRantMap2 = 16*gDet * ( m2Ant*rAntMap * (1.-rAntMap)
      - (1.-rAntMap)*m2I - rAntMap*m2K )
      + ( pow2(s02) - pow2(s02min) )
      * ( pow2(sAnt) - pow2(sAntMin) );
    if(bigRantMap2 < 0.){
      stringstream ss;
      ss<<"On line "<<__LINE__;
      infoPtr->errorMsg("Warning in "+__METHOD_NAME__
        +": kinematics map is broken.",ss.str());
      return false;
    }
    double bigRantMap = sqrt( bigRantMap2 );
    double p1dotpI = (sig2*(pow2(s02) - pow2(s02min))*
      (m2Ant + pow2(mass0) - pow2(mass1) - pow2(mass2) - s12)
      +8*rAntMap*(m2Ant + pow2(mass0) - pow2(mass1) - pow2(mass2) - s12)*gDet
      -bigRantMap*(pow2(s02) - pow2(s02min) + s01*s02-2*s12*pow2(mass0)))
      /(4*(4*gDet + m2Ant*(pow2(s02) - pow2(s02min))));

    // Norm of the three-momentum and the energy of the first parent
    // particle (could also be obtained by explicitly boosting
    // incoming particles to CM).
    double apInum2 = pow2(m2Ant) + pow2(m2I) + pow2(m2K) - 2*m2Ant*m2I
      - 2*m2Ant*m2K - 2*m2I*m2K;
    if (apInum2 < 0.) {
      stringstream ss;
      ss<<"On line "<<__LINE__;
      infoPtr->errorMsg("Warning in "+__METHOD_NAME__
        +": kinematics map is broken.",ss.str());
      return false;
    }
    double apI = sqrt(apInum2)/(2*mAnt);
    double EI = sqrt( pow2(apI) + m2I );

    // Protect against rounding errors before taking acos.
    double cospsi = ((E0*EI) - p1dotpI)/(ap0*apI);
    if (cospsi >= 1.0) {
      psi = 0.;
    } else if (cospsi <= -1.0) {
      psi = M_PI;
    } else if(isnan(cospsi)){
      psi= 0.;
      stringstream ss;
      ss << "ap0 = " << ap0;
      ss << " apI = " << apI;
      ss << " E0 = " << E0;
      ss << " mass0 = " << mass0;
      ss << " mAnt = " << mAnt;
      ss << " sum1 = " << pow2(sAnt) + pow2(m2I) + pow2(m2K);
      ss << " sum2 = " << 2*sAnt*m2I + 2*sAnt*m2K + 2*m2I*m2K;
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": cos(psi) = nan.",ss.str());
      return false;
    }
    else{
      psi = acos( cospsi );
    }
  }

  // Perform global rotations.
  p1.rot(psi,phi);
  p2.rot(psi,phi);
  p3.rot(psi,phi);

  // Verbose output.
  if (verbose >= 7) {
    Vec4 tot = p1+p2+p3;
    printOut("VinciaCommon::map2to3FFmassive:", "phi = " + num2str(phi,6)
             + " cospsi = " + num2str(cos(psi),6) +" psi = " + num2str(psi,6));
    printOut("VinciaCommon::map2to3FFmassive:",
             "mapType = " + num2str(kMapType));
    printOut("VinciaCommon::map2to3FFmassive:","final momenta in CM");
    cout << " k1cm = " << p1 << " k2cm = " << p2 << " k3cm = " << p3
         << " total = " << tot << endl;
  }

  // Rotate and boost to lab frame.
  RotBstMatrix M;
  M.fromCMframe(pTwo[0],pTwo[1]);
  if (verbose >= superdebug) {
    Vec4 tot = pTwo[0]+pTwo[1];
    cout << " (VinciaCommon::map2to3FFmassive) boosting to LAB frame "
         << "defined by\n";
    cout << " p1 =   " << pTwo[0] << " p2 =   " << pTwo[1]
         << " total = " << tot << endl;
  }
  p1.rotbst(M);
  p2.rotbst(M);
  p3.rotbst(M);
  if (verbose >= 7) {
    Vec4 tot = p1+p2+p3;
    cout << " (VinciaCommon::map2to3FFmassive:) final momenta in LAB\n";
    cout << " k1 =   " << p1 << " k2 =   " << p2 << " k3 =   " << p3
         << " total = " << tot << endl;
  }

  // Save momenta.
  pThree.resize(0);
  pThree.push_back(p1);
  pThree.push_back(p2);
  pThree.push_back(p3);

  Vec4 total = pTwo[0] + pTwo[1];
  total -= (p1+p2+p3);
  if (abs(total.e()) > SMALL || abs(total.px()) > SMALL
    || abs(total.py()) > SMALL || abs(total.pz()) >  SMALL ){
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+
      ": Failed momentum conservation test. Aborting.");
    infoPtr->setAbortPartonLevel(true);
    return false;
  }
  if (isnan(total)) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": (E,p) = nan.");
    return false;
  }

  // Return.
  return true;

}

//--------------------------------------------------------------------------

// 2->3 branching kinematics, massless with output arguments first,
// then input arguments. The output is p3, while the inputs are
// kMapType, invariants(sIK, s01, s12), and phi.

bool VinciaCommon::map2to3FFmassless(vector<Vec4>& pThree,
  const vector<Vec4>& pTwo, int kMapType, const vector<double>& invariants,
  double phi) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");

  // Antenna invariant mass.
  double m2Ant = m2(pTwo[0],pTwo[1]);
  double mAnt  = sqrt(m2Ant);

  // Compute invariants (massless relation).
  double s01 = invariants[1];
  double s12 = invariants[2];
  double s02 = m2Ant - s01 - s12;

  // Can check alternative hadronization vetos here.
  if (verbose >= 7) {
    cout << " (VinciaCommon::map2to3FFmassless:) m  = " << num2str(mAnt)
         << "   m12 =" << num2str(sqrt(s01))
         << "   m23 =" << num2str(sqrt(s12))
         << "   m13 =" << num2str(sqrt(s02)) << endl;
    RotBstMatrix M;
    Vec4 p1cm = pTwo[0];
    Vec4 p2cm = pTwo[1];
    M.toCMframe(p1cm,p2cm);
    p1cm.rotbst(M);
    p2cm.rotbst(M);
    Vec4 tot = p1cm+p2cm;
    cout << " (VinciaCommon::map2to3FFmassless:) starting dipole in CM\n"
         << " p1cm = " << p1cm << " p2cm = " << p2cm
         << " total = " << tot<<endl;
  }

  // Set up kinematics in rest frame.
  double E0 = 1/mAnt*(s01/2 + s02/2);
  double E1 = 1/mAnt*(s01/2 + s12/2);
  double E2 = 1/mAnt*(s02/2 + s12/2);
  double ap0 = E0;
  double ap1 = E1;
  double ap2 = E2;
  double cos01 = (E0*E1 - s01/2)/(ap0*ap1);
  double cos02 = (E0*E2 - s02/2)/(ap0*ap2);

  // Protection: num. precision loss for small (ultracollinear) invariants.
  if ( 1-abs(cos01) < 1e-15 ) cos01 = cos01 > 0 ? 1. : -1.;
  if ( 1-abs(cos02) < 1e-15 ) cos02 = cos02 > 0 ? 1. : -1.;
  double sin01 = (abs(cos01) < 1) ? sqrt(abs(1.0 - pow2(cos01))) : 0.0;
  double sin02 = (abs(cos02) < 1) ? sqrt(abs(1.0 - pow2(cos02))) : 0.0;

  // Set momenta in CMz frame (with 1 oriented along positive z axis
  // and event in (x,z) plane).
  Vec4 p1(0.0,0.0,ap0,E0);
  Vec4 p2(-ap1*sin01,0.0,ap1*cos01,E1);
  Vec4 p3(ap2*sin02,0.0,ap2*cos02,E2);

  // Verbose output.
  if (verbose >= superdebug) {
    Vec4 tot = p1+p2+p3;
    cout << " (map2to3FFmassless:) configuration in CM* (def: 1 along z)\n"
         << " k1* =  " << p1 << " k2* =  " << p2 << " k3* =  " << p3
         << " total = " << tot << endl;
  }

  // Choose global rotation around axis perpendicular to event plane.
  double psi;

  // Force A to be longitudinal recoiler.
  if (kMapType == -2) {
    psi = 0.0;

  // Force B to be longitudinal recoiler.
  } else if (kMapType == -1) {
    psi = M_PI - acos(max(-1.,min(1.,cos02)));

  // ARIADNE map.
  } else if (kMapType == 1) {
    psi = E2*E2/(E0*E0+E2*E2)*(M_PI-acos(cos02));

  // psi PYTHIA-like. "Recoiler" remains along z-axis.
  } else if (kMapType == 2) {
    psi = 0.;
    if (s01 < s12 || (s01 == s12 && rndmPtr->flat() > 0.5) )
      psi = M_PI-acos(cos02);

  // Kosower's map. Similar to ARIADNE.
  } else {
    double rMap(1);
    if (kMapType == 3) rMap = s12/(s01+s12);
    double rho=sqrt(1.0+4.0*rMap*(1.0-rMap)*s01*s12/s02/m2Ant);
    double s00=-( (1.0-rho)*m2Ant*s02 + 2.0*rMap*s01*s12 ) / 2.0 /
      (m2Ant - s01);
    psi=acos(1.0+2.0*s00/(m2Ant-s12));
  }

  // Perform global rotations.
  p1.rot(psi,phi);
  p2.rot(psi,phi);
  p3.rot(psi,phi);

  // Verbose output.
  if (verbose >= 7) {
    Vec4 tot = p1+p2+p3;
    printOut(__METHOD_NAME__, "phi = " + num2str(phi,6) + "psi = " +
             num2str(psi,6));
    printOut(__METHOD_NAME__, "final momenta in CM");
    cout << " k1cm = " << p1 << " k2cm = " << p2 << " k3cm = " << p3
         << " total = " << tot << endl;
  }

  // Rotate and boost to lab frame.
  RotBstMatrix M;
  M.fromCMframe(pTwo[0],pTwo[1]);
  Vec4 total = pTwo[0] + pTwo[1];
  if (verbose >= superdebug) {
    cout  << " (VinciaCommon::map2to3FFmassless:) boosting to LAB frame "
          << "defined by\n" << " p1 =   " << pTwo[0] << " p2 =   " << pTwo[1]
          << " total = " << total << endl;
  }
  p1.rotbst(M);
  p2.rotbst(M);
  p3.rotbst(M);
  if (verbose >= 7) {
    Vec4 tot = p1 + p2 + p3 ;
    printOut(__METHOD_NAME__,"final momenta in LAB");
    cout <<" k1 =   "<<p1<<" k2 =   "<<p2<<" k3 =   "<<p3
         <<" total = "<<tot<<endl;
  }

  // Save momenta.
  pThree.resize(0);
  pThree.push_back(p1);
  pThree.push_back(p2);
  pThree.push_back(p3);

  // Check momentum conservation.
  Vec4 diff = total - (p1+p2+p3);
  if(abs(diff.e())  / abs(total.e()) > SMALL ||
     abs(diff.px()) / abs(total.e()) > SMALL ||
     abs(diff.py()) / abs(total.e()) > SMALL ||
     abs(diff.pz()) / abs(total.e()) > SMALL) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": (E,p) not conserved.","Aborting.");
    cout << setprecision(10) << " difference = " << total.px() << " "
         << total.py() << " " << total.pz() << " " << total.e() << endl;
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // Return.
  return true;

}

//--------------------------------------------------------------------------

// Implementations of RF clustering maps for massive partons.

bool VinciaCommon::map2to3RFmassive(vector<Vec4>& pThree, vector<Vec4> pTwo,
  vector<double> invariants,double phi,
  vector<double> masses) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");

  // Get momenta and boost to lab frame.
  if(pTwo.size() != 2){
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Wrong number of momenta provided.");
    return false;
  }

  // Momentum of recoiler(s), final state parton, and (modified) resonance.
  Vec4 pAKBefore = pTwo.at(0);
  Vec4 pKBefore  = pTwo.at(1);
  Vec4 pABefore  = pKBefore + pAKBefore;
  Vec4 pACoM     = pABefore;

  // Boost to CoM frame of (modified) resonance.
  pKBefore.bstback(pABefore);
  pAKBefore.bstback(pABefore);
  pACoM.bstback(pABefore);

  // Get the polar and phi angle in CoM frame of K.
  double thetaK = pKBefore.theta();
  double phiK = pKBefore.phi();


  // Construct recoiled momenta in (modified) resonance CoM
  // frame. Recover masses and unscaled invariants.
  double saj = invariants[1];
  double sjk = invariants[2];
  double sak = invariants[3];
  double mA  = masses[0];
  double mj  = masses[1];
  double mk  = masses[2];
  double mAK = masses[3];
  double invDiff = mA*mA + mj*mj + mk*mk -saj-sak+sjk;
  invDiff -= mAK*mAK;

  // Get energies.
  double EjAfter = saj/(2.0*mA);
  double EkAfter = sak/(2.0*mA);
  if (EkAfter < mk)  return false;
  else if (EjAfter < mj) return false;
  else if (invDiff > SMALL) return false;

  // Get cosTheta.
  double cosTheta = getCosTheta(EjAfter,EkAfter, mj,mk, sjk);
  if (abs(cosTheta) > 1.0) return false;
  double sinTheta = sqrt(1.0 - cosTheta*cosTheta);
  double pk = sqrt(EkAfter*EkAfter-mk*mk);
  double pj = sqrt(EjAfter*EjAfter-mj*mj);

  // Construct three momenta, assuming k was along z.
  Vec4 pkAfter(0.,0.,pk, EkAfter);
  Vec4 pjAfter(pj*sinTheta,0.,pj*cosTheta, EjAfter);
  Vec4 pajkAfter = pACoM - pkAfter - pjAfter;

  // Give some transverse recoil to k.
  double thetaEff = -(M_PI-pajkAfter.theta());
  pkAfter.rot(thetaEff,0.);
  pjAfter.rot(thetaEff,0.);
  pajkAfter.rot(thetaEff,0.);

  // Rotate by arbitrary phi.
  pkAfter.rot(0.,phi);
  pjAfter.rot(0.,phi);
  pajkAfter.rot(0.,phi);

  // Rotate to recover original orientation of frame.
  pkAfter.rot(thetaK,phiK);
  pjAfter.rot(thetaK,phiK);
  pajkAfter.rot(thetaK,phiK);

  // Boost to lab frame.
  pkAfter.bst(pABefore);
  pjAfter.bst(pABefore);
  pajkAfter.bst(pABefore);
  pThree.clear();
  pThree.push_back(pajkAfter);
  pThree.push_back(pjAfter);
  pThree.push_back(pkAfter);

  // Return.
  return true;

}

//--------------------------------------------------------------------------

// Implementations of resonance kineatic maps for massive partons. Inputs
// are as follows:
//   pBefore    = momenta of resonance and  all downstream recoilers
//                before emission.
//   posF       = position in pBefore of the momentum of the F end of the
//                antenna.
//   invariants = yaj and yjk scaled invariants.
//   phi        = azimuthal angle of gluon emission.
//   mui        = masses of a, j, k.
// The output is as follows:
//  pAfter = momenta of resonance, emission and all downstream recoilers
//           after emission.
//           [0]   = pa - will be unchanged
//           [1]   = pj
//           [2]   = pk
//           [i>3] = recoilers

bool VinciaCommon::map2toNRFmassive(vector<Vec4>& pAfter, vector<Vec4> pBefore,
  unsigned int posR, unsigned int posF, vector<double> invariants,double phi,
  vector<double> masses) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");
  pAfter.clear();

  // Momentum of "R", "F" end of antenna, and sum of downstream recoilers.
  Vec4 pR = pBefore.at(posR);
  Vec4 pF = pBefore.at(posF);
  Vec4 pSum(0.,0.,0.,0.);
  vector<Vec4> pRec;
  for(unsigned int imom = 0; imom < pBefore.size(); imom++){
    if (imom==posF || imom==posR) {
      continue;
    } else {
      pSum+= pBefore.at(imom);
      pRec.push_back(pBefore.at(imom));
    }
  }
  vector<Vec4> pTwo;
  vector<Vec4> pThree;
  pTwo.push_back(pSum);
  pTwo.push_back(pF);

  // Recoil AK system.
  if (!map2to3RFmassive(pThree,pTwo,invariants,phi,masses)) {
    return false;
  } else if (pThree.size() != 3) {
    return false;
  }

  // Add pa, pj, and k. Check mass.
  pAfter.push_back(pR);
  pAfter.push_back(pThree.at(1));
  pAfter.push_back(pThree.at(2));
  Vec4 pSumAfter = pThree.at(0);
  if (abs(pSumAfter.mCalc() - pSum.mCalc()) > SMALL) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Failed to conserve mass of system.");
    return false;
  }

  // If only a single recoiler, it just takes the remaining momentum.
  if (pRec.size() == 1) {
    pAfter.push_back(pSumAfter);

  // Boost the downstream recoilers appropriately
  } else {
    for(unsigned int imom = 0; imom < pRec.size(); imom++) {
      double mRecBef = pRec.at(imom).mCalc();
      pRec.at(imom).bstback(pSum,pSum.mCalc());
      pRec.at(imom).bst(pSumAfter,pSum.mCalc());
      double mRecAfter = pRec.at(imom).mCalc();

      // Check mass.
      if (abs(mRecAfter- mRecBef) > SMALL) {
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed to conserve mass of recoilers.");
        return false;
      }
      pAfter.push_back(pRec.at(imom));
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// 2 -> 3 kinematics map for initial-initial antennae, for general mj.

bool VinciaCommon::map2to3II(vector<Vec4>& pNew, vector<Vec4>& pRec,
  vector<Vec4>& pOld, double sAB, double saj, double sjb, double sab,
  double phi, double mj2) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");

  // Hand off to massless map if mj2 = 0.
  if (mj2 == 0.0)
    return map2to3IImassless(pNew, pRec, pOld, sAB, saj, sjb, sab, phi);

  // Do massive mapping.
  pNew.clear();
  pNew.resize(3);

  // Force incoming momenta on shell (massless) with mass squared = sAB.
  pOld[0].py(0.);
  pOld[0].px(0.);
  pOld[1].py(0.);
  pOld[1].px(0.);
  double sCM = m2(pOld[0] + pOld[1]);
  double fac = sqrt(sAB/sCM);
  double e0 = pOld[0].e();
  double e1 = pOld[1].e();
  if (abs(1. - fac) > TINY) {
    if (verbose >= 3 && abs(1.-fac) > 1.01)
      printOut("VinClu::map2to3II", "Warning: scaling AB so m2(AB) = sAB");
    e0 *= fac;
    e1 *= fac;
  }
  int sign = pOld[0].pz() > 0 ? 1 : -1;
  pOld[0].e(e0);
  pOld[0].pz(sign * e0);
  pOld[1].e(e1);
  pOld[1].pz(-sign * e1);

  // Initialise new momenta.
  pNew[0] = pOld[0];
  pNew[2] = pOld[1];

  // Check if we're inside massive phase space.
  double G = saj*sjb*sab - mj2*sab*sab;
  if (G < 0. || sab < 0.) return false;
  if ((sab <= sjb) || (sab <= saj)) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Incompatible invariants.");
    return false;
  }

  // Incoming momenta.
  double rescaleFacA = sqrt(sab/sAB * (sab - saj)/(sab - sjb));
  double rescaleFacB = sqrt(sab/sAB * (sab - sjb)/(sab - saj));
  pNew[0].rescale4(rescaleFacA);
  pNew[2].rescale4(rescaleFacB);

  // Emission.
  double preFacA = sjb*sqrt((sab - saj)/(sab - sjb)/sab/sAB);
  double preFacB = saj*sqrt((sab - sjb)/(sab - saj)/sab/sAB);
  double preFacT = sqrt(saj*sjb/sab - mj2);
  Vec4 pTrans(cos(phi), sin(phi), 0.0, 0.0);
  pNew[1] = preFacA*pOld[0] + preFacB*pOld[1] + preFacT*pTrans;

  if (verbose >= superdebug) {
    printOut("VinClu::map2to3II","Invariants are");
    cout << scientific << "    sAB = " << sAB << " saj = " << saj
         << " sjb = " << sjb << " sab = " << sab << endl
         << " Given momenta are" << endl;
    for (int i = 0; i < 2; i++) cout << "    " << pOld[i];
    cout << " New momenta are" << endl;
    for (int i = 0; i < 3; i++) cout << "    " << pNew[i];
  }

  // Check the invariants, allow 0.0001% difference.
  double check = 1e-6;
  double sajNew = 2*pNew[0]*pNew[1];
  double sjbNew = 2*pNew[1]*pNew[2];
  double sabNew = 2*pNew[0]*pNew[2];
  if (abs(sabNew - sab)/sab > check) {
    if (verbose >= 5) {
      printOut("VinClu::map2to3II","ERROR! Invariants differ!");
      cout << scientific << " sab (" << sab << ") fracdiff = ydiff = "
           << abs(sabNew-sab)/sab << endl << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
    }
    return false;
  } else if (abs(sajNew - saj)/sab > check) {
    if (verbose >= 5) {
      printOut("VinClu::map2to3II","ERROR! Invariants differ!");
      cout << scientific << " saj (" << saj << ") fracdiff = "
           << abs(sajNew-saj)/saj << " ydiff = "
           << abs(sajNew-saj)/sab << endl << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
    }
    return false;
  } else if (abs(sjbNew - sjb)/sab > check) {
    if (verbose >= 5) {
      printOut("VinClu::map2to3II","ERROR! Invariants differ!");
      cout << scientific << " sjb (" << sjb << ") fracdiff = "
           << abs(sjbNew-sjb)/sjb << " ydiff = " << abs(sjbNew-sjb)/sab
           << endl << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
    }
    return false;
  }

  // Change the final state recoiler. The recoiler is currently sum
  // of initial guys => E,0,0,pz. Boost in center-of-mass frame AB
  // E,0,0,0.
  Vec4 pSum = pOld[0] + pOld[1];
  Vec4 pRecSumBefore(0.,0.,0.,0.);
  Vec4 pRecSumAfter(0.,0.,0.,0.);
  for (int i=0; i<(int)pRec.size(); i++){
    pRecSumBefore+=pRec[i];
    pRec[i].bstback(pSum);
  }

  // Now boost from E,0,0,0 to E',px',py',pz'.
  Vec4 pPrime = pNew[0] + pNew[2] - pNew[1];
  for (int i=0; i<(int)pRec.size(); i++) {
    pRec[i].bst(pPrime, pPrime.mCalc());
    pRecSumAfter+=pRec[i];
  }
  if(verbose >= superdebug) {
    Vec4 total= pOld[0]+pOld[1];
    cout << " Total In before" << total <<  endl
         << " Total Out before" << pRecSumBefore <<  endl;
    total = pNew[0] + pNew[2] - pNew[1];
    cout << " Total In After" << total <<  endl
         << " Total Out After" << pRecSumAfter <<  endl
         << " Total diff After" << total-pRecSumAfter <<  endl;
  }
  return true;

}

//--------------------------------------------------------------------------

// 2 -> 3 kinematics map for initial-initial antennae, for mj= 0.

bool VinciaCommon::map2to3IImassless(vector<Vec4>& pNew, vector<Vec4>& pRec,
  vector<Vec4>& pOld, double sAB, double saj, double sjb, double sab,
  double phi) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");
  pNew.clear();
  pNew.resize(3);

  // Force incoming momenta on shell (massless) with mass squared = sAB.
  pOld[0].py(0.);
  pOld[0].px(0.);
  pOld[1].py(0.);
  pOld[1].px(0.);
  double sCM = m2(pOld[0] + pOld[1]);
  double fac = sqrt(sAB/sCM);
  double e0 = pOld[0].e();
  double e1 = pOld[1].e();
  if (abs(1. - fac) > TINY) {
    if (verbose >= 3 && abs(1. - fac) > 1.01)
      printOut("VinciaCommon::map2to3IImassless",
               "Warning: scaling AB so m2(AB) = sAB");
    e0 *= fac;
    e1 *= fac;
  }
  int sign = pOld[0].pz() > 0 ? 1 : -1;
  pOld[0].e(e0);
  pOld[0].pz(sign * e0);
  pOld[1].e(e1);
  pOld[1].pz(-sign * e1);

  // Initialise new momenta.
  pNew[0] = pOld[0];
  pNew[2] = pOld[1];

  // Incoming momenta.
  double rescaleFacA = sqrt(sab/(sAB+saj) * (1. + sjb/sAB));
  double rescaleFacB = sqrt(sab/(sAB+sjb) * (1. + saj/sAB));
  pNew[0].rescale4(rescaleFacA);
  pNew[2].rescale4(rescaleFacB);

  // Emission.
  double preFacA = sjb*sqrt((sAB+sjb)/(sAB+saj)/sab/sAB);
  double preFacB = saj*sqrt((sAB+saj)/(sAB+sjb)/sab/sAB);
  double preFacT = sqrt(saj*sjb/sab);
  Vec4 pTrans(cos(phi), sin(phi), 0.0, 0.0);
  pNew[1] = preFacA*pOld[0] + preFacB*pOld[1] + preFacT*pTrans;

  // Debugging info.
  if (verbose >= superdebug) {
    stringstream ss;
    ss << "Invariants are: "
       << scientific << "    sAB = " << sAB << " saj = " << saj
       << " sjb = " << sjb << " sab = " << sab;
    printOut(__METHOD_NAME__, ss.str());
    printOut(__METHOD_NAME__, "Given momenta are");
    for (int i = 0; i < 2; i++) cout << "    " << pOld[i];
    printOut(__METHOD_NAME__, "New momenta are");
    for (int i = 0; i < 3; i++) cout << "    " << pNew[i];
  }

  // Check the invariants, allow 0.0001% difference.
  double check = 1e-6;
  double sajNew = 2*pNew[0]*pNew[1];
  double sjbNew = 2*pNew[1]*pNew[2];
  double sabNew = 2*pNew[0]*pNew[2];
  if (abs(sabNew - sab)/sab > check) {
    if (verbose >= 5) {
      printOut("VinciaCommon::map2to3IImassless","ERROR! Invariants differ!");
      cout << scientific << " sab (" << sab << ") fracdiff = ydiff = "
           << abs(sabNew - sab)/sab << endl
           << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
    }
    return false;
  } else if (abs(sajNew - saj)/sab > check) {
    if (verbose >= 5) {
      printOut("VinciaCommon::map2to3IImassless","ERROR! Invariants differ!");
      cout << scientific << " saj (" << saj << ") fracdiff = "
           << abs(sajNew-saj)/saj << " ydiff = "<< abs(sajNew - saj)/sab
           << endl << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
    }
    return false;
  } else if ( abs(sjbNew-sjb)/sab > check ) {
    if (verbose >= 5) {
      printOut("VinciaCommon::map2to3IImassless","ERROR! Invariants differ!");
      cout << scientific << " sjb (" << sjb << ") fracdiff = "
           << abs(sjbNew-sjb)/sjb << " ydiff = "<< abs(sjbNew - sjb)/sab
           << endl << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
    }
    return false;
  }

  // Change the final state recoiler. The recoiler is currently sum
  // of initial guys => E,0,0,pz. Boost in center-of-mass frame AB
  // E,0,0,0.
  Vec4 pSum = pOld[0] + pOld[1];
  for (int i=0; i<(int)pRec.size(); i++) pRec[i].bstback(pSum);

  // Now boost from E,0,0,0 to E',px',py',pz' and return.
  Vec4 pPrime = pNew[0] + pNew[2] - pNew[1];
  for (int i=0; i<(int)pRec.size(); i++) pRec[i].bst(pPrime);
  return true;

}

//--------------------------------------------------------------------

// 2->3 kinematics map for local recoils, for general mj,mk. Assumes
// partons from proton explicitly massless

bool VinciaCommon::map2to3IFlocal(vector<Vec4>& pNew, vector<Vec4>& pOld,
  double sAK, double saj, double sjk, double sak, double phi,
  double mK2, double mj2, double mk2) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");
  pNew.clear();
  pNew.resize(3);
  if (verbose >= superdebug) {
    printOut("VinClu::map2to3IFlocal","Invariants are");
    cout << "    sAK = " << sAK << " saj = " << saj
         << " sjk = " << sjk << " sak = " << sak << endl
         << "    mK = " << sqrt(mK2) << " mj = " << sqrt(mj2)
         << " mk = " << sqrt(mk2) << endl
         << " Given momenta are" << endl;
    for (int i=0; i<2; i++) cout << "    " << pOld[i];
  }

  // Check invariants.
  double check = 1.e-3;
  double inv1Norm = (saj + sak)/(sAK + sjk);
  double inv2Norm = 1.0  + (mj2 + mk2 - mK2)/(sAK + sjk);
  double diff = abs(inv1Norm-inv2Norm);
  if(diff > check) {
    if (verbose >= 2) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Inconsistent invariants.","Aborting.");
      cout <<" yaj + yak = " << inv1Norm
           << " 1 + muj2 + muk2 - muK2 = "<< inv2Norm
           << " Diff = " << diff << endl;
    }
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // Check if we're inside massive phase space.
  double G = saj*sjk*sak - mj2*sak*sak - mk2*saj*saj;
  if (G < 0. || sak < 0.) return false;

  // Set up some variables for boosting pT.
  Vec4 pSum = pOld[0] + pOld[1];
  Vec4 pOldBst = pOld[0];
  pOldBst.bstback(pSum);
  double thetaRot = pOldBst.theta();
  double phiRot = pOldBst.phi();
  Vec4 pTrans(cos(phi), sin(phi), 0.0, 0.0);

  // Rotate and boost.
  pTrans.rot(thetaRot, phiRot);
  pTrans.bst(pSum);

  // Check if pT was properly boosted, allow 0.1% difference.
  if (pTrans*pOld[0] > pOld[0][0]*1e-3 || pTrans*pOld[1] > pOld[1][0]*1e-3) {
    if (verbose >= normal) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": The transverse momentum is not transverse after boosting");
    }
    return false;
  }
  double sig = sak + saj;
  double cj = (sig*(sak + mj2 - mk2) + mK2*(sak - saj) - sAK*sak)/(sAK*sig);
  double ck = (sig*(saj - mj2 + mk2) + mK2*(saj - sak) - sAK*saj)/(sAK*sig);
  double dj = saj/sig;
  double dk = sak/sig;

  // Construct new momenta p_a, p_j, and p_k.
  pNew[0] = pOld[0];
  pNew[0].rescale4(sig/sAK);
  pNew[1] = cj*pOld[0] + dj*pOld[1] + (sqrt(G)/sig)*pTrans;
  pNew[2] = ck*pOld[0] + dk*pOld[1] - (sqrt(G)/sig)*pTrans;

  // Check the invariants, allow 0.1% difference (due to boost).
  double sakNew = pNew[0]*pNew[2]*2;
  double sajNew = pNew[0]*pNew[1]*2;
  double sjkNew = pNew[1]*pNew[2]*2;
  if (abs(sakNew - sak)/sak > check) {
    if (verbose >= 2) {
      printOut("VinClu::map2to3IFlocal","ERROR! sak is inconsistent!");
      cout << scientific << " sak (" << sak << ") diff = "
           << abs(sakNew-sak)/sak << endl
           << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
      cout <<"Masses:    mK2 = " << mK2 << " mj2 = " << mj2
           << " mk2 = " << mk2 << endl;
    }
    return false;
  } else if (abs(sajNew - saj)/saj > check) {
    if (verbose >= 2 ) {
      printOut("VinClu::map2to3IFlocal","ERROR! saj is inconsistent!");
      cout << scientific << " saj (" << saj << ") diff = ";
      cout << abs(sajNew-saj)/saj << endl;
      cout << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
      cout <<"Masses:    mK2 = " << mK2 << " mj2 = " << mj2
           << " mk2 = " << mk2 << endl;
    }
    return false;
  } else if ( abs(sjkNew-sjk)/sjk > check ) {
    if (verbose >= 2 ){
      printOut("VinClu::map2to3IFlocal","ERROR! sjk is inconsistent!");
      cout << scientific << " sjk (" << sjk << ") diff = ";
      cout << abs(sjkNew-sjk)/sjk << endl;
      cout << " Old momenta are" << endl;
      for (int i=0; i<2; i++) cout << "    " << pOld[i];
      cout << " New momenta are" << endl;
      for (int i=0; i<3; i++) cout << "    " << pNew[i];
      cout <<"Masses:    mK2 = " << mK2 << " mj2 = " << mj2
           << " mk2 = " << mk2 << endl;
    }
    infoPtr->setAbortPartonLevel(true);
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------

// 2->3 kinematics map for global recoils, for general mj,mk.  Assumes
// partons from proton explicitly massless.

bool VinciaCommon::map2to3IFglobal(vector<Vec4>& pNew,
  vector<Vec4>& pRec, vector<Vec4>& pOld, Vec4 &pB,
  double sAK, double saj, double sjk, double sak, double phi,
  double mK2, double mj2, double mk2) {

  if (verbose >= superdebug) printOut(__METHOD_NAME__, "begin --------------");
  pNew.clear();
  pNew.resize(3);

  // Set up some variables for boosting pT.
  Vec4 pSum = pOld[0] + pOld[1];
  Vec4 pAcm = pOld[0];
  pAcm.bstback(pSum);
  double thetaRot = pAcm.theta();
  double phiRot = pAcm.phi();
  Vec4 pTrans(cos(phi), sin(phi), 0.0, 0.0);

  // Rotate and boost.
  pTrans.rot(thetaRot, phiRot);
  pTrans.bst(pSum);

  // Check if pT was properly boosted, allow 0.1% difference.
  if (pTrans*pOld[0] > pOld[0].e()*1e-3 || pTrans*pOld[1] > pOld[1].e()*1e-3) {
    if (verbose >= superdebug) {
      printOut("VinciaCommon:map2to3IFglobal",
               "The transverse momentum is not transverse after boosting");
    }
    return false;
  }

  // Solution vector, start from massless values.
  vector<double> v; v.resize(5);
  double sig = sAK - saj;
  v[0] = sak/sig;
  v[1] = (saj*sjk)/(sAK*sig);
  v[2] = sjk/sig;
  v[3] = (saj*sak)/(sAK*sig);
  v[4] = sjk*saj*sak/pow2(sig);

  // Root finding with Newton-Raphson in 5D.
  int nCount = 0;
  while (true) {
    nCount ++;

    // Construct function.
    vector<double> f(5, 0);
    f[0] = v[0]*v[1]*sAK + pow2(v[1])*mK2 - v[4];
    f[1] = v[2]*v[3]*sAK + pow2(v[3])*mK2 - v[4] - mj2;
    f[2] = (v[0] - v[2] - 1)*(v[1] - v[3] + 1)*sAK
      + pow2(v[1] - v[3] + 1)*mK2 - mk2;
    f[3] = (v[0]*v[3] + v[1]*v[2])*sAK + 2*v[1]*v[3]*mK2 - 2*v[4] - saj;
    f[4] = (v[2]*(v[1] - v[3] + 1) + v[3]*(v[0] - v[2] - 1))*sAK
      + 2*v[3]*(v[1] - v[3] + 1)*mK2 - sjk;

    // Construct Jacobian.
    vector<vector<double> > A(5, vector<double>(5, 0));
    A[0][0] = v[1]*sAK;
    A[0][1] = v[0]*sAK + 2*v[1]*mK2;
    A[0][2] = 0;
    A[0][3] = 0;
    A[0][4] = -1;
    A[1][0] = 0;
    A[1][1] = 0;
    A[1][2] = v[3]*sAK;
    A[1][3] = v[2]*sAK + 2*v[3]*mK2;
    A[1][4] = -1;
    A[2][0] = (v[1] - v[3] + 1)*sAK;
    A[2][1] = (v[0] - v[2] - 1)*sAK + 2*(v[1] - v[3] + 1)*mK2;
    A[2][2] = -(v[1] - v[3] + 1)*sAK;
    A[2][3] = -( (v[0] - v[2] - 1)*sAK + 2*(v[1] - v[3] + 1)*mK2 );
    A[2][4] = 0;
    A[3][0] = v[3]*sAK;
    A[3][1] = v[2]*sAK + 2*v[3]*mK2;
    A[3][2] = v[1]*sAK;
    A[3][3] = v[0]*sAK + 2*v[1]*mK2;
    A[3][4] = -2;
    A[4][0] = v[3]*sAK;
    A[4][1] = v[2]*sAK + 2*v[3]*mK2;
    A[4][2] = (v[1] - 2*v[3] + 1)*sAK;
    A[4][3] = (v[0] - 2*v[2] - 1)*sAK + 2*(v[1] - 2*v[3] + 1)*mK2;
    A[4][4] = 0;

    // Invert Jacobian and append identity.
    int n = 5;
    for (int i = 0; i < n; i++) {
      A[i].resize(2*n);
      A[i][n+i] = 1;
    }

    for (int i = 0; i < n; i++) {
      // Find column max.
      double eleMax = abs(A[i][i]);
      int rowMax = i;
      for (int k=i+1; k<n; k++) {
        if (abs(A[k][i]) > eleMax) {
          eleMax = A[k][i];
          rowMax = k;
        }
      }

      // Swap maximum row with current row.
      A[rowMax].swap(A[i]);

      // Reduce current column.
      for (int k = i+1; k < n; k++) {
        double c = -A[k][i]/A[i][i];
        for (int j = i; j < 2*n; j++) {
          if (i == j) {
            A[k][j] = 0;
          } else {
            A[k][j] += c * A[i][j];
          }
        }
      }
    }

    // Solve equation Ax = b for upper triangular matrix A.
    for (int i = n-1; i >= 0; i--) {
      for (int k = n; k < 2*n; k++) {
        A[i][k] /= A[i][i];
      }
      for (int j = i-1; j >= 0; j--) {
        for (int k = n; k < 2*n; k++) {
          A[j][k] -= A[i][k] * A[j][i];
        }
      }
    }

    // Remove remaining identity.
    for (int i = 0; i < n; i++) {
      A[i].erase(A[i].begin(), A[i].begin()+n);
    }

    // Find next iteration.
    vector<double> vNew(5, 0);
    for (int i = 0; i < n; i++) {
      vNew[i] = v[i];
      for (int j=0; j<n; j++) {
        vNew[i] -= A[i][j]*f[j];
      }
    }

    // Check if done.
    double eps = 0;
    for (int i=0; i<n; i++) {
      if (abs(vNew[i] - v[i])/abs(v[i]) > eps) {
        eps = abs(vNew[i] - v[i])/abs(v[i]);
      }
    }
    v = vNew;

    // Break if iterations/precision reached or randomly vary the variables.
    if (nCount == 1000) {return false;}
    if ((nCount%100) == 0) {
      for (int i=0; i<n; i++) {
        v[i] *= (2*rndmPtr->flat() - 1);
      }
    }
    if (eps < 1E-6) {break;}
  }

  // Construct post-branching momenta.
  pNew[0] = v[0]*pOld[0] + v[1]*pOld[1] - sqrt(v[4])*pTrans;
  pNew[1] = v[2]*pOld[0] + v[3]*pOld[1] - sqrt(v[4])*pTrans;
  pNew[2] = (v[0] - v[2] - 1)*pOld[0] + (v[1] - v[3] + 1)*pOld[1];

  // Set up the boost.
  Vec4 pa = pNew[0];
  Vec4 pA = pOld[0];
  double qaB = pa*pB;
  double qAB = pA*pB;
  double qAa = pA*pa;

  // Perform boost.
  for (int i=0; i<3; i++) {
    Vec4 p = pNew[i];
    pNew[i] += pB*((pa*p)/qaB) - pa*((pB*p)/qaB) + pA*((pB*p)/qAB)
      - pB*((pA*p)/qAB) + pB*(qAa*(pB*p)/(qAB*qaB));

    // Force the initial state to be on the beam axis.
    if (i==0) {
      double ea = pNew[i].e();
      double sign = (pNew[0].pz() > 0) ? 1 : -1;
      pNew[0] = Vec4(0, 0, sign*ea, ea);
    }
  }

  // Perform boost on the rest of the system and return.
  for (int i=0; i<(int)pRec.size(); i++) {
    Vec4 p = pRec[i];
    pRec[i] += pB*((pa*p)/qaB) - pa*((pB*p)/qaB) + pA*((pB*p)/qAB)
      - pB*((pA*p)/qAB) + pB*(qAa*(pB*p)/(qAB*qaB));
  }
  double sajNew = 2*pNew[0]*pNew[1];
  double sjkNew = 2*pNew[1]*pNew[2];
  double sakNew = 2*pNew[0]*pNew[2];
  if (verbose >= 5) {
    if (abs(sajNew - saj)/saj > 1E-3)
      printOut("VinciaCommon:map2to3IFglobal", "saj not quite correct");
    if (abs(sjkNew - sjk)/sjk > 1E-3)
      printOut("VinciaCommon:map2to3IFglobal", "sjk not quite correct");
    if (abs(sakNew - sak)/sak > 1E-3)
      printOut("VinciaCommon:map2to3IFglobal", "sak not quite correct");
  }
  return true;
}

//--------------------------------------------------------------------------

// pA + pK -> pa + pj + pk

bool VinciaCommon::map2to3RFmassive(vector<Vec4>& pNew, vector<Vec4>& pRec,
  vector<Vec4> pOld, double saj, double sjk, double phi,
  double mA2 = 0.0, double mj2 = 0.0, double mK2 = 0.0) {

  // Get momenta and boost to lab frame.
  if (pOld.size() != 2) {return false;}
  Vec4 pABefore = pOld[0];
  Vec4 pKBefore = pOld[1];
  Vec4 pAKBefore = pABefore - pKBefore;
  Vec4 pACoM = pABefore;
  double sAK = 2.0*pABefore*pKBefore;
  double sak = sAK - saj + sjk;

  // Check if inside phase space boundaries.
  if (sak < 0) {return false;}
  if (saj*sjk*sak - mA2*sjk*sjk - mj2*sak*sak - mK2*saj*saj < 0) return false;

  // Boost to CoM frame.
  pKBefore.bstback(pABefore);
  pAKBefore.bstback(pABefore);
  pACoM.bstback(pABefore);

  // Get energies, cos(theta), and sin(theta).
  double EjAfter    = saj/(2.0*sqrt(mA2));
  double pVecjAfter = sqrt(pow2(EjAfter) - mj2);
  double EkAfter    = sak/(2.0*sqrt(mA2));
  double pVeckAfter = sqrt(pow2(EkAfter) - mK2);
  double cosTheta   = (2.*EjAfter*EkAfter - sjk)/2./pVecjAfter/pVeckAfter;
  if (abs(cosTheta) > 1) {return false;}
  double sinTheta   = sqrt(1.0 - cosTheta*cosTheta);

  // Construct three momenta.
  Vec4 pkAfter(0., 0., pVeckAfter, EkAfter);
  Vec4 pjAfter(pVecjAfter*sinTheta*cos(phi), pVecjAfter*sinTheta*sin(phi),
               pVecjAfter*cosTheta, EjAfter);
  Vec4 pajkAfter = pACoM - pkAfter - pjAfter;

  // Boost to lab frame.
  pkAfter.bst(pABefore);
  pjAfter.bst(pABefore);
  pajkAfter.bst(pABefore, sqrt(mA2));
  pNew.clear();
  pNew.push_back(pABefore);
  pNew.push_back(pjAfter);
  pNew.push_back(pkAfter);

  // In case of a single recoiler, it just takes the remaining momentum.
  if (pRec.size() == 1) {
    pRec[0] = pajkAfter;

  // Perform Lorentz boost for multiple recoilers, fails for a
  // single massless recoiler
  } else {
    for(int i=0; i<(int)pRec.size(); i++) {
      pRec[i].bstback(pAKBefore);
      pRec[i].bst(pajkAfter);
    }
  }
  return true;

}

//--------------------------------------------------------------------------

// Check if 2-particle system is on-shell and rescale if not.

bool VinciaCommon::onShellCM(Vec4& p1, Vec4& p2, double m1, double m2,
  double tol) {

  double s1     = pow2(m1);
  double s2     = pow2(m2);
  double s01    = Vec4(p1+p2).m2Calc();
  double s1Calc = p1.m2Calc();
  double s2Calc = p2.m2Calc();
  if (abs(s1Calc-s1)/s01 > tol || abs(s2Calc-s2)/s01 > tol) {
    if (verbose >= 3)
      printOut("VinClu::onShellCM","forcing particles on mass shell");
    RotBstMatrix M;
    M.fromCMframe(p1,p2);

    // Define massive on-shell momenta.
    double E0 = (s01 + s1 - s2)/(2*sqrt(s01));
    double E1 = (s01 - s1 + s2)/(2*sqrt(s01));
    double pz = pow2(E0)-s1;
    Vec4 p1new = Vec4(0.0,0.0,-pz,E0);
    Vec4 p2new = Vec4(0.0,0.0,pz,E1);
    p1new.rotbst(M);
    p2new.rotbst(M);
    double s1Test = p1new.m2Calc();
    double s2Test = p2new.m2Calc();
    if (verbose >= 3) {
      cout << " p1   : " << p1 << " p1new: " << p1new
           << " p2   : " << p1 << " p2new: " << p1new;
    }

    // If this got them closer to mass shell, replace momenta.
    if (abs(s1Test-s1)/s01 <= abs(s1Calc-s1)/s01
      && abs(s2Test-s2)/s01 <= abs(s2Calc-s2)/s01) {
      p1 = p1new;
      p2 = p2new;
    }
    return false;
  }
  else return true;

}

//--------------------------------------------------------------------------

// Map partons partonSystems[iSys] to equivalent massless ones.
// Return true if succeeded. Note, a method using only Particles or
// Vec4 as input could in principle be split off from this, if needed,
// but has not been required so far.

bool VinciaCommon::mapToMassless(int iSys, Event& event,
  PartonSystems* partonSystemsPtr, bool makeNewCopies) {

  // Start by making new copies, if requested to do so.
  if (makeNewCopies) {
    int iOld, iNew;

    // Copy incoming partons, interpret the copying operation as the
    // two incoming partons recoiling off each other. Assign status
    // code -42, incoming copy of recoiler (as mother).
    if (partonSystemsPtr->hasInAB(iSys)) {
      iOld = partonSystemsPtr->getInA(iSys);
      iNew = event.copy(iOld,-42);
      partonSystemsPtr->replace(iSys, iOld, iNew);
      iOld = partonSystemsPtr->getInB(iSys);
      iNew = event.copy(iOld,-42);
      partonSystemsPtr->replace(iSys, iOld, iNew);
    }
    // Note, a decaying resonance is not copied to preserve structure
    // of production and decay. Copy outgoing partons (use status code
    // 52).
    for (int i=0; i<partonSystemsPtr->sizeOut(iSys); ++i) {
      iOld = partonSystemsPtr->getOut(iSys,i);
      iNew = event.copy(iOld, 52);
    } // End loop to make new copies.
  } // End if new copies requested.

  // Initial-state partons, always assumed massless in VINCIA.
  if (partonSystemsPtr->hasInAB(iSys) ) {
    int iA = partonSystemsPtr->getInA(iSys);
    int iB = partonSystemsPtr->getInB(iSys);
    if (event[iA].mCalc() != 0.0 || event[iB].mCalc() != 0.0) {

      // Below we assume iA is the one with pz > 0; swap if opposite case.
      if (event[iA].pz() < 0 || event[iB].pz() > 0) {
        iA = partonSystemsPtr->getInB(iSys);
        iB = partonSystemsPtr->getInA(iSys);
      }

      // Transverse components assumed zero: check.
      if (event[iA].pT() > 1.e-6 || event[iB].pT() > 1.e-6) {
        cout<<"Vincia::VinciaCommon::mapToMassless) Error: incoming partons "
          "have nonvanishing transverse momenta for system iSys = "
            <<iSys<<"; giving up"<<endl;
        return false;
      }

      // Verbose output.
      if (verbose > superdebug) {
        stringstream ss;
        ss<<"Warning: forcing initial"
          "-state partons to be massless for system "<<iSys;
        printOut(__METHOD_NAME__,ss.str());
      }

      // Define explicitly massless momenta (same as in Pythia::PartonLevel).
      double pPos = event[iA].pPos() + event[iB].pPos();
      double pNeg = event[iA].pNeg() + event[iB].pNeg();
      event[iA].pz( 0.5 * pPos);
      event[iA].e ( 0.5 * pPos);
      event[iA].m ( 0.);
      event[iB].pz(-0.5 * pNeg);
      event[iB].e ( 0.5 * pNeg);
      event[iB].m ( 0.);
    }
  } // End make initial-state partons massless.

  // Final-state partons.
  if (partonSystemsPtr->sizeOut(iSys) >= 2) {
    vector<Vec4> momenta;
    vector<double> massOrg;
    bool makeMassless = false;
    Vec4 pSysOrg;
    for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
      momenta.push_back(event[partonSystemsPtr->getOut(iSys,i)].p());
      massOrg.push_back(event[partonSystemsPtr->getOut(iSys,i)].m());
      if (massOrg[i] > 0. && event[partonSystemsPtr->getOut(iSys,i)].idAbs()
        <= nFlavZeroMass) makeMassless = true;
      pSysOrg += momenta[i];
    }
    // Return if nothing needs to be done.
    if (!makeMassless) return true;

    // Create copy of original momenta (in case of failure).
    vector<Vec4> momentaOrg = momenta;

    // Boost to CM if original system not at rest.
    double sCM = m2(pSysOrg);
    bool isInCM = ( pow2(pSysOrg.pAbs())/sCM < 1e-10 );
    if (!isInCM)
      for (int i=0; i<(int)momenta.size(); ++i) momenta[i].bstback(pSysOrg);

    // Define vector for computing CM energy of modified system.
    Vec4 pSysNew;

    // Identify particles to be made massless (by ID code) and rescale
    // their momenta along direction of motion.
    for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
      int ip = partonSystemsPtr->getOut(iSys,i);
      if (event[ip].idAbs() <= nFlavZeroMass && event[ip].m() != 0.) {
        double facInv = momenta[i].pAbs()/momenta[i].e();
        // Sanity check.
        if (facInv <= 0.) {
          if (verbose >= 2)
            printOut("VinciaCommon:mapToMassless",
              "Remap failed. Particle is spacelike or at rest.");
          // Restore masses in case any were already changed.
          for (int j=0; j < partonSystemsPtr->sizeOut(iSys); ++j)
            event[partonSystemsPtr->getOut(iSys,j)].m(massOrg[j]);
          // Failed.
          return false;
        }
        momenta[i].rescale3(1./facInv);
        event[ip].m(0.);
        // Check new 4-vector.
        double mNew = momenta[i].mCalc();
        if (pow2(mNew/momenta[i].e()) > TINY) {
          printOut("VinciaCommon:mapToMassless","Warning: rounding problem.");
          if (verbose >= 7) {
            cout<<scientific << "(p,e) = "<<momenta[i].pAbs() << "  "
                << momenta[i].e() << " facInv = " << facInv
                << " 1/facInv = " << 1./facInv << " mNew = " << mNew << endl;
          }
        } // End check new 4-vector.
      } // End if massless flavour with mass > 0.

      // Add to new CM momentum.
      pSysNew += momenta[i];
    } // End loop over FS particles.

    // New system generally has smaller invariant mass and some
    // motion. Determine if additional scalings or boosts are needed.
    Vec4 delta = pSysOrg - pSysNew;
    if (delta.e()/sqrt(sCM) < TINY && delta.pAbs()/sqrt(sCM) < TINY) {
      // Update event record (masses already updated above).
      for (int i = 0; i < (int)momenta.size(); ++i)
        event[partonSystemsPtr->getOut(iSys,i)].p(momenta[i]);
      if (verbose >= superdebug)
        printOut("VinciaCommon:mapToMassless","No further rescalings needed.");

      //Check momentum conservation.
      if (!checkCoM(iSys,event,partonSystemsPtr)) {
        infoPtr->setAbortPartonLevel(true);
        infoPtr->errorMsg("Error in "+__METHOD_NAME__+
          ": Failed (E,p) conservation check.","Aborting.");
        return false;
      }
      return true;
    }

    // If the new system has a different CM energy, rescale all
    // energies and momenta to restore the same CM energy as before.
    double sCMnew   = m2(pSysNew);
    double scaleFac = sqrt(sCM/sCMnew);
    if (verbose >= 7 && pow2(scaleFac-1.0) > TINY)
      printOut("VinciaCommon:mapToMassless",
               "Rescaling 4-vectors to restore eCM");
    Vec4 pSysNewB;
    for (int i = 0; i < (int)momenta.size(); ++i) {
      momenta[i].rescale4(scaleFac);
      pSysNewB += momenta[i];
    }
    double sCMnewB = m2(pSysNewB);
    if (verbose >= superdebug)
      cout << "old CM energy = " << sqrt(sCM) << " intermediate CM energy = "
           << sqrt(sCMnew) << " new CM energy = " << sqrt(sCMnewB) << endl;
    // Then boost to CM frame (preserves CM energy).
    if (verbose >= 7)
      printOut("VinciaCommon:mapToMassless","Boosting back to CM frame");
    for (int i=0; i<(int)momenta.size(); ++i) {
      // Boost to new CM frame
      momenta[i].bstback(pSysNewB);
      // If required, also boost back to frame of input system
      if (!isInCM) momenta[i].bst(pSysOrg);

      // Update event record (masses already updated above).
      event[partonSystemsPtr->getOut(iSys,i)].p(momenta[i]);

    } // End do boosts.

    // Verbose output: final configuration.
    if (verbose >= 7) {
      cout << "Final configuration:" << endl;
      for (int i = 0; i < (int)momenta.size(); ++i)
        cout << "  " << i << " " << momenta[i];
    }

  } // End make final-state momenta massless.

  // Check momentum conservation.
  if(!checkCoM(iSys, event,partonSystemsPtr)){
    infoPtr->setAbortPartonLevel(true);
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Failed (E,p) conservation check.","Aborting.");
    return false;
  }

  // We made it.
  return true;

}

//--------------------------------------------------------------------------

// More lightweight function to check conservation of momentum.

bool VinciaCommon::checkCoM(int iSys, Event& event,
  PartonSystems* partonSystemsPtr){
  Vec4 total(0.,0.,0.,0.);
  if (!partonSystemsPtr->hasInRes(iSys)){
    if (partonSystemsPtr->getInA(iSys) > 0)
      total+= event[partonSystemsPtr->getInA(iSys)].p();
    if (partonSystemsPtr->getInB(iSys) > 0)
      total+= event[partonSystemsPtr->getInB(iSys)].p();
  } else total+= event[partonSystemsPtr->getInRes(iSys)].p();
  double sysMass = total.mCalc();

  // Loop over members of current system.
  for( int iPart=0; iPart<partonSystemsPtr->sizeOut(iSys); iPart++){
    int iOut = partonSystemsPtr->getOut(iSys,iPart);

    // Sum total FS momentum.
    if (event[iOut].isFinal()) total -= event[iOut].p();
    else {
      stringstream ss;
      ss << "iSys = " << iSys << " iOut = " << iOut;
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": isFinal()=false for outgoing parton.",ss.str());
      partonSystemsPtr->list();
      event.list();
      return false;
    }
  }
  total/=sysMass;
  if(abs(total.e()) > SMALL || abs(total.px()) > SMALL
     || abs(total.py()) > SMALL || abs(total.pz()) >  SMALL) {
    event.list();
    cout << "total = " << setprecision(10) << total.e() << " " << total.px()
         << " " << total.py() << " " << total.pz() << endl;
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +" Failed (E,p) conservation check.");
    return false;
  } else if(isnan(total)){
    event.list();
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +" Failed (E,p) isnan check.");
    return false;
  } else return true;

}

//==========================================================================

// TXiFunctor helper class.

//--------------------------------------------------------------------------

// Constructor.

TXiFunctor::TXiFunctor(vector<double> mIn, vector<double> energiesIn) {
  if (mIn.size() != energiesIn.size()) {
    m = vector<double>(0); energies = vector<double>(0);
  } else {m = mIn; energies = energiesIn;}
}

double TXiFunctor::operator() (double xi) {
  double retval = 0.;
  for (vector<double>::size_type i = 0; i < m.size(); i++)
    retval += sqrt( pow2(m[i]) + pow2(xi)*pow2(energies[i]));
  return retval;
}

//==========================================================================

// VINCIA Auxiliary helper functions.

//--------------------------------------------------------------------------

// External auxiliaries, extra four-products.

double m(const Vec4& v) {
  double s = m2(v); return (s >= 0) ? sqrt(s) : -sqrt(-s);}

double m2(const Vec4& v) {
  return pow2(v.e()) - pow2(v.px()) - pow2(v.py()) - pow2(v.pz());}

double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3) {
  return pow2(v1.e() + v2.e() + v3.e())
    - pow2(v1.px() + v2.px() + v3.px())
    - pow2(v1.py() + v2.py() + v3.py())
    - pow2(v1.pz() + v2.pz() + v3.pz());}

double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Vec4& v4) {
  return pow2(v1.e() + v2.e() + v3.e() + v4.e())
    - pow2(v1.px() + v2.px() + v3.px() + v4.px())
    - pow2(v1.py() + v2.py() + v3.py() + v4.py())
    - pow2(v1.pz() + v2.pz() + v3.pz() + v4.pz());}

double m2(const Particle& p1, const Particle& p2, const Particle& p3) {
  return m2(p1.p(), p2.p(), p3.p());}

double dot4(const Particle& p1, const Particle& p2) {return p1.p()*p2.p();}

double getCosTheta(double E1, double E2, double m1, double m2, double s12){
  return  (2.0*E1*E2 - s12)/(2.0*sqrt(E1*E1 - m1*m1)*sqrt(E2*E2 - m2*m2));}

//--------------------------------------------------------------------------

// External auxiliaries, string manipulation.

string num2str(int i, int width) {
  ostringstream tmp;
  if (width <= 1) tmp << i;
  else if (abs(i) < pow(10.0, width - 1) || ( i > 0 && i < pow(10.0, width)))
    tmp << fixed << setw(width) << i;
  else {
    string ab = "k";
    double r = i;
    if      (abs(i) < 1e5)       {r/=1e3;}
    else if (abs(i) < 1e8)  {r/=1e6;  ab = "M";}
    else if (abs(i) < 1e11) {r/=1e9;  ab = "G";}
    else if (abs(i) < 1e14) {r/=1e12; ab = "T";}
    tmp << fixed << setw(width - 1)
        << (r > 10 ? setprecision(width-4) : setprecision(width-3)) << r << ab;
  }
  return tmp.str();
}

string num2str(double r, int width) {
  ostringstream tmp;
  if (width <= 0) tmp << r;
  else if (r == 0.0 || (abs(r) > 0.1 && abs(r) < pow(10., max(width-3,1)))
           || width <= 8) tmp << fixed << setw(max(width,3))
                              << setprecision(min(3, max(1, width - 2))) << r;
  else tmp << scientific << setprecision(max(2, width - 7))
           << setw(max(9, width)) << r;
  return tmp.str();
}

string bool2str(bool b, int width) {
  string tmp = b ? "on" : "off";
  int nPad = width - tmp.length();
  for (int i = 1; i <= nPad; ++i) tmp = " " + tmp;
  return tmp;
}

void printOut(string place, string message) {
  cout.setf(ios::internal);
  cout << " (" << (place + ") ") << message << "\n";
}

//--------------------------------------------------------------------------

// Gram determinant, invariants used in the argument = 2*pi*pj.

double gramDet( double s01tilde, double s12tilde, double s02tilde,
  double m0, double m1, double m2) {
  return ((s01tilde*s12tilde*s02tilde - pow2(s01tilde)*pow2(m2)
           - pow2(s02tilde)*pow2(m1) - pow2(s12tilde)*pow2(m0))/4
          + pow2(m0)*pow2(m1)*pow2(m2));
}

double gramDet(Vec4 p0, Vec4 p1, Vec4 p2) {
  return gramDet(2*p0*p1, 2*p1*p2, 2*p0*p2, p0.mCalc(), p1.mCalc(),
    p2.mCalc());
}

//--------------------------------------------------------------------------

// Math support auxiliaries.

// Dilogarithm.
double Li2(const double x, const double kmax, const double xerr) {
  if (x < 0.0) return 0.5*Li2(x*x) - Li2(-x);
  if (x <= 0.5) {
    double sum(x), term(x);
    for (int k = 2; k < kmax; k++) {
      double rk = (k-1.0)/k;
      term *= x*rk*rk;
      sum += term;
      if (abs(term/sum) < xerr) return sum;
    }
    cout << "Maximum number of iterations exceeded in Li2" << endl;
    return sum;
  }
  if (x < 1.0)  return M_PI*M_PI/6.0 - Li2(1.0 - x) - log(x)*log(1.0 - x);
  if (x == 1.0) return M_PI*M_PI/6.0;
  if (x <= 1.01) {
    const double eps(x - 1.0), lne(log(eps)),
      c0(M_PI*M_PI/6.0),         c1(  1.0 - lne),
      c2(-(1.0 - 2.0*lne)/4.0),  c3( (1.0 - 3.0*lne)/9.0),
      c4(-(1.0 - 4.0*lne)/16.0), c5( (1.0 - 5.0*lne)/25.0),
      c6(-(1.0 - 6.0*lne)/36.0), c7( (1.0 - 7.0*lne)/49.0),
      c8(-(1.0 - 8.0*lne)/64.0);
    return c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*(c5 + eps*(
                     c6 + eps*(c7 + eps*c8)))))));
  }
  double logx = log(x);
  if (x<=2.0) return M_PI*M_PI/6.0 + Li2(1.0 - 1.0/x) -
                logx*(log(1.0 - 1.0/x) + 0.5*logx);
  return M_PI*M_PI/3.0 - Li2(1.0/x) - 0.5*logx*logx;
}


// Standard factorial.
double factorial(const int n) {
  double fac = 1;
  for (int i = 2; i <= n; i++) fac *= i;
  return fac;}

// Binomial coefficient.
int binomial(const int n, const int m) {
  if (m < 0 || m > n) return 0;
  else if (m == n || m == 0) return 1;
  else if (m == 1 || m == n - 1) return n;
  else return factorial(n)/factorial(m)/factorial(n - m) + 0.01;
}

// Lambert W function using the rational fit from Darko Veberic's
// paper, arXiv:1209.0735v2.  Should give 5 digits of precision for
// positive arguments x not too large (fit region was 0.3, 2e, but
// still has 5-digit accuracy at zero).  Precision quickly drops for
// negative values, but he has extra functions that can be implemented
// if those are needed, and for very large values the asymptotic
// log(x), log(log(x)) form could be used if precise solutions for
// large values are needed. For now just write a warning if we are
// ever asked for a value far outside region of validity.
double LambertW(const double x) {
  if (x == 0.) return 0.;
  if (x < -0.2) {
    cout << "Warning in "<<__METHOD_NAME__
         << ": Accuracy less than three decimal places for x < -0.2";
  } else if (x > 10.) {
    cout << "Warning in "<<__METHOD_NAME__
         <<": Accuracy less than three decimal places for x > 10.";
  }
  return x*(1. + x*(2.445053 + x*(1.343664 + x*(0.14844 + 0.000804*x))))
    /(1. + x*(3.444708 + x*(3.292489 + x*(0.916460 + x*(0.053068)))));
}

// Version of zbrent using function pointers, solve fun(x) - r = 0 for x.
double zbrent(TFunctor& fun, double r, double x1, double x2, double tol) {
  int iter;
  double a(x1), b(x2), c(x2), d(x2-x1), e(x2-x1), min1, min2;
  double fa(fun(a) - r), fb(fun(b) - r), fc, p, q, r1, s, tol1, xm;
  double REALTINY = min(TINY, 1e-12);
  tol = max(tol, REALTINY);

  // Check if there is a single zero in range.
  if (fa*fb > 0) return 0.0;

  // Start search.
  fc = fb;
  for (iter = 1; iter < max(1000, int(1.0/sqrt(tol))); iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a; fc = fa; e = d = b-a;
    }
    if (abs(fc) < abs(fb)) {
      a = b; b = c; c= a; fa = fb; fb = fc; fc = fa;
    }
    tol1 = 2.0*REALTINY*abs(b) + 0.5*tol;
    xm = 0.5*(c-b);
    if (abs(xm) <= tol1 || fb == 0.0) return b;
    if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
      s = fb/fa;
      if (a == c) {p = 2.0*xm*s; q = 1.0-s;}
      else {
        q = fa/fc; r1 = fb/fc;
        p = s*(2.0*xm*q*(q - r1) - (b - a)*(r1 - 1.0));
        q = (q - 1.0)*(r1 - 1.0)*(s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = abs(p);
      min1 = 3.0*xm*q - abs(tol1*q);
      min2 = abs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {e = d; d= p/q;}
      else {d = xm; e = d;}
    } else {d = xm; e = d;}
    a = b;
    fa = fb;
    b += abs(d) > tol1 ? d : (xm > 1) ? tol1 : - tol1;
    fb = fun(b) - r;
  }
  cerr << "(brent:) -> Maximum number of iterations exceeded" << endl;
  return 0.0;

}

//==========================================================================

} // end namespace Pythia8
