// VinciaFSR.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the VinciaFSR class
// and auxiliary classes.

#include "Pythia8/VinciaFSR.h"

namespace Pythia8 {

//==========================================================================

// The Brancher class, base class containing a generic set of "parent
// partons" as well as virtual methods for generating trial
// branchings.

//--------------------------------------------------------------------------

// Reset.

void Brancher::reset(int iSysIn, Event& event, vector<int> iIn) {

  // Save info on parents and resize vectors.
  iSav         = iIn;
  hasTrialSav  = false;
  systemSav    = iSysIn;
  Vec4 pSum;
  int nMassive = 0;
  idSav.resize(iIn.size());
  hSav.resize(iIn.size());
  colTypeSav.resize(iIn.size());
  colSav.resize(iIn.size());
  acolSav.resize(iIn.size());
  mSav.resize(iIn.size());
  for (unsigned int i = 0; i < iIn.size(); ++i) {
    idSav[i]      = event[iIn[i]].id();
    hSav[i]       = event[iIn[i]].pol();
    colTypeSav[i] = event[iIn[i]].colType();
    colSav[i]     = event[iIn[i]].col();
    acolSav[i]    = event[iIn[i]].acol();
    mSav[i]       = event[iIn[i]].m();
    if (mSav[i] != 0.0) nMassive += 1;
    // Compute and store antenna invariant mass.
    pSum += event[iIn[i]].p();
  }
  m2AntSav     = pSum.m2Calc();
  mAntSav      = (m2AntSav >= 0) ? sqrt(m2AntSav) : -sqrt(-m2AntSav);
  // Massless parents: sIK = m2IK and kallenFac = 1.0.
  sAntSav      = m2AntSav;
  kallenFacSav = 1.0;
  // Mass corrections to sAnt and kallenFac.
  if (nMassive != 0) {
    // sIK = m2IK - m2I - m2K.
    for (unsigned int i = 0; i < iIn.size(); ++i) sAntSav -= pow2(mSav[i]);
    // Phase-space correction non-unity if both parents massive.
    // Note, so far only defined for 2-parton systems.
    if (nMassive == 2 && iIn.size() == 2)
      kallenFacSav = sAntSav/sqrt(pow2(sAntSav) - 4*pow2(mSav[0]*mSav[1]));
  }
  // Call method to initialise variables in derived classes.
  init();

}

//--------------------------------------------------------------------------

// Compute pT scale of trial branching.
double Brancher::getpTscale() {

  if (invariantsSav.size() == 3) {
    double sIK = invariantsSav[0];
    double y12 = invariantsSav[1] / sIK;
    double y23 = invariantsSav[2] / sIK;
    return sIK * y12 * y23;
  } else return 0.;

}

//--------------------------------------------------------------------------

// Return Xj.
double Brancher::getXj() {

  if (invariantsSav.size() == 3) {
    double sIK = invariantsSav[0];
    double y12 = invariantsSav[1] / sIK;
    double y23 = invariantsSav[2] / sIK;
    return y12 + y23;
  } else return 1.0;

}

//--------------------------------------------------------------------------

// Simple print utility, showing the contents of the Brancher.

void Brancher::list(string header) const {

  // Check if we are asked to output a header.
  if (header != "none")
    cout << " --------  " << std::left << setw(30) << header
         << "  ----------------"
         << "--------------------------------------------- \n"
         << "  sys type     mothers         colTypes              ID codes "
         << "            hels             m    qNewSav \n"
         << fixed << std::right << setprecision(3);
  cout << setw(5) << system() << " ";
  string type = "FF";
  if (iSav.size() == 3) type = "FFF";
  else if (iSav.size() >= 4) type="FS";
  cout << setw(4) << type << " "
       << setw(5) << i0() << " " << setw(5) << i1() << " "
       << setw(5) << (i2() > 0 ? num2str(i2(), 5) : " ") << "   "
       << setw(3) << colType0() << " " << setw(3) << colType1() << " "
       << setw(3) << (i2() > 0 ? num2str(colType2(), 3) : " ") << " "
       << setw(9) << id0() << setw(9) << id1()
       << setw(9) << (i2() > 0 ? num2str(id2(), 9) : " ") << "   "
       << setw(2) << h0() << " " << setw(2) <<  h1() << " "
       << setw(2) << (i2() > 0 ? num2str(h2(), 2) : " ") << " "
       << num2str(mAnt(), 10);
  if (hasTrial()) {
    if (q2NewSav > 0.) cout << " " << num2str(sqrt(q2NewSav), 10);
    else cout << " " << num2str(0.0, 10);
  }
  else cout << " " << setw(10) << "-";
  cout << endl;

}

//--------------------------------------------------------------------------

// Set post-branching IDs and masses. Base class is for gluon emission.

void Brancher::setidPost() {
  idPostSav.clear();
  idPostSav.push_back(id0());
  idPostSav.push_back(21);
  idPostSav.push_back(id1());
}

vector<double> Brancher::setmPostVec() {
  mPostSav.clear();
  mPostSav.push_back(mSav[0]); // mi
  mPostSav.push_back(0.0);     // mj
  mPostSav.push_back(mSav[1]); // mk
  return mPostSav;
}

void Brancher::setStatPost() {
  statPostSav.resize(iSav.size() + 1, 51);}

void Brancher::setMaps(int){
  mothers2daughters.clear(); daughters2mothers.clear();}

//--------------------------------------------------------------------------

// Return index of new particle (slightly arbitrary choice for splittings).

int Brancher::iNew(){

  if (i0() > 0) {
    if(mothers2daughters.find(i0()) != mothers2daughters.end())
      return mothers2daughters[i0()].second;
    else return 0;
  }
  else return 0;

}


//==========================================================================

// Class BrancherEmitFF, branch elemental for 2->3 gluon emissions.

//--------------------------------------------------------------------------

// Method to initialise members specific to BrancherEmitFF.

void BrancherEmitFF::init() {

  branchType = 1;
  if (colType0() == 2 && colType1() == 2) iAntSav = iGGemitFF;
  else if (colType1() == 2) iAntSav = iQGemitFF;
  else if (colType0() == 2) iAntSav = iGQemitFF;
  else iAntSav = iQQemitFF;

}

//--------------------------------------------------------------------------

// Generate a new Q2 value, soft-eikonal 2/yij/yjk implementation.

double BrancherEmitFF::genQ2(int evTypeIn, double q2BegIn, Rndm* rndmPtr,
  const EvolutionWindow* evWindowIn, double colFacIn,
  vector<double> headroomIn, vector<double> enhanceIn,
  int) {

  // Initialise output value and save input parameters.
  q2NewSav       = 0.;
  evTypeSav      = evTypeIn;
  evWindowSav    = evWindowIn;
  colFacSav      = colFacIn;
  q2BegSav       = q2BegIn;
  headroomSav    = (headroomIn.size() >=1) ?  headroomIn[0] : 1.0 ;
  enhanceSav     = (enhanceIn.size() >=1) ?  enhanceIn[0] : 1.0 ;
  // Set flag that this call will produce a saved trial.
  hasTrialSav    = true;
  // Overall normalisation factor (independent of evType).
  double normFac = colFacSav * kallenFac() * headroomSav * enhanceSav;
  // pT evolution.
  if (evTypeSav == 1) {
    // Constant trial alphaS.
    if (evWindowSav->runMode <= 0) {
      double coeff  = normFac * evWindowSav->alphaSmax / 4. / M_PI;
      double lnR    = log(rndmPtr->flat());
      double xT1    = q2BegSav / sAnt();
      double ln2xT2 = pow2(log(xT1)) - lnR / coeff;
      double xT2    = exp(-sqrt(ln2xT2));
      q2NewSav      = xT2 * sAnt();

    // Running trial alphaS
    } else {
      double q2Min  = pow2(evWindowSav->qMin);
      double d      = sqrt(1 - 4 * q2Min/sAnt());
      double deltaY = log( (1.+d)/(1.-d) );
      double coeff  = normFac * deltaY / 2. / M_PI;
      double b0     = evWindowSav->b0;
      double powR   = pow(rndmPtr->flat(),b0/coeff);
      double facLam = evWindowSav->lambda2 / evWindowSav->kMu2;
      q2NewSav = facLam * pow(q2BegSav/facLam,powR);
    }
  }
  if (q2NewSav > q2BegIn) {
    string errorMsg = "Error in "+__METHOD_NAME__
      +": Generated q2New > q2BegIn."+" Returning 0.";
    cout<<errorMsg<<endl;
    q2NewSav = 0.;
  }
  return q2NewSav;

}

//--------------------------------------------------------------------------

// Generate invariants.

bool BrancherEmitFF::genInvariants(vector<double>& invariants,
  Rndm* rndmPtr, int){

  // Clear output vector, check if we have a sensible q2New scale.
  invariants.clear();
  if (q2NewSav <= 0.) return false;

  // pT evolution.
  if (evTypeSav == 1) {
    double xT  = q2NewSav / sAnt();
    if (xT > 0.25) return false;
    // Rapidity boundaries for constant trial alphaS.
    double yMinTrial, yMaxTrial;
    if (evWindowSav->runMode <= 0) {
      yMinTrial = log(xT)/2.;
      yMaxTrial = -yMinTrial;

    // Rapidity boundaries for running trial alphaS
    } else {
      double q2MinNow = pow2(evWindowSav->qMin);
      double xTMinNow = q2MinNow/sAnt();
      double dTrial = sqrt(1. - 4*xTMinNow);
      yMaxTrial = 0.5 * log( (1.+dTrial)/(1.-dTrial ) );
      yMinTrial = -yMaxTrial;
    }
    double yTrial    = yMinTrial + rndmPtr->flat()*(yMaxTrial - yMinTrial);
    double dPhys     = sqrt(1. - 4*xT);
    double yMaxPhys  = 0.5 * log((1. + dPhys)/(1. - dPhys));
    double yMinPhys  = -yMaxPhys;
    if (yTrial < yMinPhys || yTrial > yMaxPhys) return false;
    double r = exp(yTrial);
    double rootXT = sqrt(xT);
    double yij = rootXT / r;
    double yjk = rootXT * r;
    double yik = 1.0  - yij - yjk;
    if (yij < 0. || yjk < 0. || yik < 0.) {
      cout << " Problem in genInvariants yij = " << yij << " yjk = "
           << yjk << endl;
      return false;
    }
    double sij = sAnt() * yij;
    double sjk = sAnt() * yjk;
    double sik = sAnt() * yik;

    // Element 0 is the antenna sAnt.
    invariants.push_back(sAnt());
    invariants.push_back(sij);
    invariants.push_back(sjk);
    invariants.push_back(sik);

    // Save the generated invariants (for future query and use by pAccept).
    invariantsSav = invariants;
    setmPostVec();

    // Veto if the point outside the available phase space.
    double det = gramDet(sij,sjk,sik,mPostSav[0],mPostSav[1],mPostSav[2]);
    if (det > 0.) return true;
    else return false;
  }
  else return false;

}

//--------------------------------------------------------------------------

// Compute antPhys / antTrial for gluon emissions, given antPhys.

double BrancherEmitFF::pAccept(const double antPhys, int) {

  // Evaluate trial function taking headroom factor into account.
  if (invariantsSav.size() <= 2) return 0.;
  double antTrial = headroomSav
    *2.*invariantsSav[0]/invariantsSav[1]/invariantsSav[2];
  // pT evolution.
  if (evTypeSav == 1) {
    // Constant trial alphaS.
    if (evWindowSav->runMode <= 0)
      antTrial *= colFacSav * evWindowSav->alphaSmax;
    // Running trial alphaS: pars = colFac, b0, facLam = Lambda^2/kMuR^2.
    else {
      double b0         = evWindowSav->b0;
      double facLam     = evWindowSav->lambda2/evWindowSav->kMu2;
      double lnx        = log(q2NewSav/facLam);
      double alphaTrial = 1./b0/lnx;
      antTrial *= colFacSav * alphaTrial;
    }
  }
  return antPhys/antTrial;

}

//--------------------------------------------------------------------------

// Return the maximum Q2.

double BrancherEmitFF::getQ2Max(int evType) {

  if      (evType == 1) return sAntSav/4.;
  else if (evType == 2) return sAntSav/9.;
  else if (evType == 3) return sAntSav/2.;
  else return 0.;

}

//--------------------------------------------------------------------------

// Method to make mothers2daughters and daughters2mothers pairs.

void BrancherEmitFF::setMaps(int sizeOld) {

  mothers2daughters.clear();
  daughters2mothers.clear();

  //For updating the children of existing parents.
  mothers2daughters[i0()] = make_pair(sizeOld, sizeOld + 1);
  mothers2daughters[i1()] = make_pair(sizeOld + 1, sizeOld + 2);

  //For adding mothers of new children.
  daughters2mothers[sizeOld]   = make_pair(i0(), 0);
  daughters2mothers[sizeOld+1] = make_pair(i0(), i1());
  daughters2mothers[sizeOld+2] = make_pair(i1(), 0);

}

//--------------------------------------------------------------------------

// Generic getter method. Assumes setter methods called earlier.

bool BrancherEmitFF::getNewParticles(Event& event, vector<Vec4> momIn,
  vector<int> hIn, vector<Particle> &pNew, Rndm* rndmPtr, Colour* colourPtr) {

  // Initialize.
  unsigned int nPost = iSav.size() + 1;
  pNew.clear();
  pNew.resize(nPost);
  setidPost();
  setStatPost();
  double scaleNew = sqrt(q2NewSav);
  setMaps(event.size());

  // Check everything set.
  if (momIn.size() != nPost || hIn.size() != nPost ||
    mPostSav.size() != nPost || idPostSav.size() != nPost ||
    statPostSav.size() != nPost || invariantsSav.size() < 3) return false;

  // Who inherits the colour?
  double sij  = invariantsSav[1];
  double sjk  = invariantsSav[2];
  bool inh01  = colourPtr->inherit01(sij,sjk);
  int lastTag = event.lastColTag();
  vector<int> col(nPost, 0);
  vector<int> acol(nPost, 0);
  acol[0] = event[i0()].acol();
  col[0]  = event[i0()].col();
  acol[2] = event[i1()].acol();
  col[2]  = event[i1()].col();

  // Generate a new colour tag.
  int colNew = lastTag + 1 + rndmPtr->flat()*10;
  // 0 keeps colour.
  if (inh01) {
    while (colNew%10 == col[2]%10 || colNew%10 == 0)
      colNew = lastTag + 1 + rndmPtr->flat()*10;
    acol[1]=col[0];
    col[1]=colNew;
    acol[2]=colNew;
  // 2 keeps colour.
  } else {
    while (colNew%10 == acol[0]%10 || colNew%10 == 0)
      colNew = lastTag + 1 + rndmPtr->flat()*10;
    col[0]=colNew;
    acol[1]=colNew;
    col[1]=acol[2];
  }

  // Now populate particle vector.
  for (unsigned int ipart = 0; ipart < nPost; ++ipart) {
    pNew[ipart].status(statPostSav[ipart]);
    pNew[ipart].id(idPostSav[ipart]);
    pNew[ipart].pol(hIn[ipart]);
    pNew[ipart].p(momIn[ipart]);
    pNew[ipart].m(mPostSav[ipart]);
    pNew[ipart].setEvtPtr(&event);
    pNew[ipart].scale(scaleNew);
    pNew[ipart].daughters(0,0);
    pNew[ipart].col(col[ipart]);
    pNew[ipart].acol(acol[ipart]);
  }
  colTagSav = colNew;
  return true;

}

//==========================================================================

// Class BrancherSplitFF, branch elemental for 2->3 gluon splittings.

//--------------------------------------------------------------------------

// Method to initialise data members specific to BrancherSplitFF.

void BrancherSplitFF::init() {
  branchType = 2; iAntSav = iGXsplitFF; swapped = false;}

//--------------------------------------------------------------------------

// Generate a new Q2 value .

double BrancherSplitFF::genQ2(int evTypeIn, double q2BegIn,
  Rndm* rndmPtr, const EvolutionWindow* evWindowIn,
  double, vector<double> headroomFlav,
  vector<double> enhanceFlav, int verboseIn) {

  // Initialise output value and save input parameters.
  q2NewSav    = 0.;
  evTypeSav   = evTypeIn;
  q2BegSav    = q2BegIn;
  evWindowSav = evWindowIn;

  // Total splitting weight summed over flavours
  double wtSum = 0.0;
  vector<double> wtFlav;
  unsigned int nFlav = headroomFlav.size();
  if (nFlav != enhanceFlav.size()) {
    if (verboseIn >=normal) {
      string errorMsg = "Error in "+__METHOD_NAME__
        +": inconsistent size of headroom and enhancement vectors.";
      cout<<errorMsg<<endl;
    }
    return 0.;
  }

  // First check if there is any phase space open for this flavour
  for (unsigned int iFlav = 0; iFlav < nFlav; ++iFlav) {
    double mFlav = evWindowSav->mass.at(iFlav+1);
    if (mAnt() - m0() - m1() < 2.*mFlav) {
      wtFlav.push_back(0.); continue;
    } else {
      double wt = headroomFlav[iFlav] * enhanceFlav[iFlav];
      wtFlav.push_back(wt);
      wtSum += wt;
    }
  }

  // Set flag that this call will produce a saved trial.
  hasTrialSav = true;
  // Overall normalisation factor common to all trial integrals.
  double normFac = kallenFac() * wtSum;
  // pT evolution.
  if (evTypeSav == 1) {
    double deltaZeta = 1.0;
    double coeff     = normFac * deltaZeta / 8. / M_PI;
    double xT1       = q2BegSav / sAnt();
    // Constant trial alphaS.
    if (evWindowSav->runMode <= 0) {
      double aStrial = evWindowSav->alphaSmax;
      double xT2     = xT1 * pow(rndmPtr->flat(), 1./aStrial/coeff);
      q2NewSav       = xT2 * sAnt();

    // Running trial alphaS.
    } else {
      double b0     = evWindowSav->b0;
      double facLam = evWindowSav->lambda2/evWindowSav->kMu2;
      double powR   = pow(rndmPtr->flat(),b0/coeff);
      q2NewSav      = facLam * pow(q2BegSav/facLam,powR);
    }
  }

  // Select flavour.
  double ranFlav = rndmPtr->flat() * wtSum;
  for (int iFlav = nFlav - 1; iFlav >= 0; --iFlav) {
    ranFlav -= wtFlav[iFlav];
    if (ranFlav < 0) {
      idFlavSav   = iFlav+1;
      // Set quark masses.
      mFlavSav    = evWindowSav->mass.at(idFlavSav);
      // Save corresponding headroom and enhancement factors.
      enhanceSav  = enhanceFlav[iFlav];
      headroomSav = headroomFlav[iFlav];
      break;
    }
  }
  if (q2NewSav > q2BegIn) {
    string errorMsg = "Error in "+__METHOD_NAME__
      +": Generated q2New > q2Beg."+" Returning 0.";
    cout<<errorMsg<<endl;
    q2NewSav = 0.;
  }
  return q2NewSav;

}

//--------------------------------------------------------------------------

// Generate complementary invariant(s) for saved trial scale
// for gluon splitting. Return false if no physical kinematics possible.

bool BrancherSplitFF::genInvariants(vector<double>& invariants,
  Rndm* rndmPtr, int) {

  // Clear output vector, and check if we have a sensible q2New scale.
  invariants.clear();
  if (q2NewSav <= 0.) return false;

  // pT evolution.
  double m2j = pow2(mFlavSav);
  if (evTypeSav == 1) {
    // Zeta is uniform in [0,1].
    double zetaTrial = rndmPtr->flat();
    double m2qq      = q2NewSav / zetaTrial;
    double sij       = m2qq - 2*m2j;
    if (sij < 0.) return false;

    // Here i=q, j=qbar is always the case, but change def for sjk,
    // sik depending on who is colour connected to the recoiler. G
    // anticolour side (XG antenna): pT = pT(qbar) = pTj; zeta = yjk +
    // mu2q.
    double sjk, sik;
    if (isXGsav) {
      sjk = sAnt()*zetaTrial - m2j;
      sik = sAnt() - m2qq - sjk;

    // G colour side (GX antenna): pT = pT(q) = pTi; zeta = yik + mu2q.
    } else {
      sik = sAnt()*zetaTrial - m2j;
      sjk = sAnt() - m2qq - sik;
    }
    if (sjk < 0.) return false;

    // Element 0 is the antenna sAnt.
    invariants.push_back(sAnt());
    invariants.push_back(sij);
    invariants.push_back(sjk);
    invariants.push_back(sik);

    // Save the generated invariants (for future query and use by pAccept).
    invariantsSav = invariants;
    setmPostVec();

    // Veto if point outside the available phase space.
    double det = gramDet(sij,sjk,sik,mPostSav[0],mPostSav[1],mPostSav[2]);
    if (det > 0.) return true;
    else return false;
  }
  else return false;

}

//--------------------------------------------------------------------------

// Compute antPhys/antTrial for gluon splittings, given antPhys.
// Note, antPhys should be normalised to include charge and coupling
// factors.

double BrancherSplitFF::pAccept(const double antPhys, int) {

  // Evaluate trial function with headroom factor.
  if (invariantsSav.size() <= 2) return 0.;
  double antTrial = headroomSav/(2.*(invariantsSav[1] + 2*pow2(mFlavSav)));

  // pT evolution (note evType =1 for pTj, evType -1 for pTi).
  if (evTypeSav == 1) {
    double alphaTrial = evWindowSav->alphaSmax;
    if (evWindowSav->runMode >= 1) {
      double b0         = evWindowSav->b0;
      double facLam     = evWindowSav->lambda2/evWindowSav->kMu2;
      double lnx        = log(q2NewSav/facLam);
      alphaTrial        = 1./b0/lnx;
    }
    antTrial *= alphaTrial;
  }
  return antPhys/antTrial;

}

//--------------------------------------------------------------------------

// Getter and setter methods.

double BrancherSplitFF::getQ2Max(int evType) {

  if      (evType == 1) return sAntSav/4.;
  else if (evType == 2) return sAntSav;
  else if (evType == 3) return sAntSav;
  else return 0.;

}

vector<double> BrancherSplitFF::setmPostVec() {

  mPostSav.clear();
  mPostSav.push_back(mFlavSav); // mi
  mPostSav.push_back(mFlavSav); // mj
  mPostSav.push_back(mSav[1]);  // mk
  return mPostSav;

}

void BrancherSplitFF::setidPost() {

  idPostSav.clear();
  idPostSav.push_back(idFlavSav);
  idPostSav.push_back(-idFlavSav);
  idPostSav.push_back(id1());

}

void BrancherSplitFF::setStatPost() {

  statPostSav.resize(iSav.size() + 1, 51);
  statPostSav[2] = 52;

}

void BrancherSplitFF::setMaps(int sizeOld) {

  // For updating the children of existing parents.
  mothers2daughters.clear();
  daughters2mothers.clear();
  mothers2daughters[i0()] = make_pair(sizeOld, sizeOld+1);
  mothers2daughters[i1()] = make_pair(sizeOld+2,sizeOld+2);

  // For adding mothers of new children.
  daughters2mothers[sizeOld] = make_pair(i0(),0);
  daughters2mothers[sizeOld+1] = make_pair(i0(),0);
  daughters2mothers[sizeOld+2] = make_pair(i1(),i1());

}

//--------------------------------------------------------------------------

// Generic getter method. Assumes setter methods called earlier.

bool BrancherSplitFF::getNewParticles(Event& event, vector<Vec4> momIn,
  vector<int> hIn, vector<Particle> &pNew, Rndm*, Colour*) {

  // Initialize.
  unsigned int nPost = iSav.size() + 1;
  pNew.clear();
  pNew.resize(nPost);
  setidPost();
  setStatPost();
  double scaleNew = sqrt(q2NewSav);
  setMaps(event.size());

  // Check everything set.
  if (momIn.size()!=nPost || hIn.size()!=nPost ||
      mPostSav.size() !=nPost || idPostSav.size() != nPost ||
      statPostSav.size() != nPost || invariantsSav.size() < 3) return false;
  vector<int> col(nPost,0);
  vector<int> acol(nPost,0);
  acol[0] = 0;
  col[0]  = event[i0()].col();
  acol[1] = event[i0()].acol();
  col[1]  = 0;
  acol[2] = event[i1()].acol();
  col[2]  = event[i1()].col();

  // Now populate particle vector.
  for (unsigned int ipart = 0; ipart < nPost; ++ipart) {
    pNew[ipart].status(statPostSav[ipart]);
    pNew[ipart].id(idPostSav[ipart]);
    pNew[ipart].pol(hIn[ipart]);
    pNew[ipart].p(momIn[ipart]);
    pNew[ipart].m(mPostSav[ipart]);
    pNew[ipart].setEvtPtr(&event);
    pNew[ipart].scale(scaleNew);
    pNew[ipart].daughters(0,0);
    pNew[ipart].col(col[ipart]);
    pNew[ipart].acol(acol[ipart]);
  }
  colTagSav = 0;
  return true;

}

//==========================================================================

// BrancherEmitRF class for storing information on antennae between a
// coloured resonance and final state parton, and generating a new
// emission.

//--------------------------------------------------------------------------

// Method to initialise data members specific to BrancherEmitRF.

void BrancherEmitRF::init(Event& event, vector<int> allIn,
  unsigned int posResIn, unsigned int posFIn, double Q2cut) {

  // Get Pythia indices of res and final.
  posRes      = posResIn;
  posFinal    = posFIn;
  int iRes    = allIn.at(posRes);
  int iFinal  = allIn.at(posFinal);
  colFlowRtoF = event[iRes].col() == event[iFinal].col() && event[iRes].col()
    != 0;

  // Extract the momenta of the rest.
  Vec4 recoilVec(0., 0., 0., 0.);
  for (vector<int>::iterator pos = allIn.begin(); pos != allIn.end(); ++pos) {
    if ((*pos == iRes) || (*pos == iFinal)) continue;
    recoilVec += event[*pos].p();
  }

  // This is not necesssarily p(res). In the case where one particle
  // always recieves the recoil e.g. W in t->bWX it is p_t -p_X.
  Vec4 resVec = recoilVec + event[iFinal].p();

  //Calculate the masses.
  mRes = resVec.mCalc();
  mFinal = event[iFinal].mCalc();
  mRecoilers = recoilVec.mCalc();
  sAK = getsAK(mRes, mFinal, mRecoilers);

  // Calculate common prefactor to trial integral.
  kallenFacSav = (2.0*sAK/(4.0*M_PI));
  kallenFacSav /= sqrt(KallenFunction(mRes*mRes, mFinal*mFinal,
      mRecoilers*mRecoilers));

  // Calculate zeta limits.
  zetaMin= zetaMinCalc(mRes, mFinal, mRecoilers, Q2cut);
  zetaMax = zetaMaxCalc(mRes, mFinal, mRecoilers);
  // Phase space is closed.
  if (zetaMax < zetaMin) zetaIntSave = 0.;
  // Calculate zeta integral (full phase space).
  else zetaIntSave = zetaIntegral(zetaMin, zetaMax);

  // Calculate Q2max.
  Q2MaxSav = calcQ2Max(mRes, mRecoilers, mFinal);
  branchType = 5;
  // TODO: swapped should be redundant since save posRes, posFinal.
  // R = Q.
  if (abs(colTypeSav[posRes]) == 1) {
    // F = Q.
    if (abs(colTypeSav[posFinal]) == 1) {
      iAntSav = iQQemitRF;
      swapped = false;
    // F = g.
    } else if (colTypeSav[posFinal] == 2) {
      iAntSav = iQGemitRF;
      swapped = posRes != 0;
    // Some other final state - don't know what to do with this yet!
    } else {
      iAntSav = -1;
      swapped = false;
    }
  // Some other resonance. Don't know what to do with this yet!
  } else {
    iAntSav = -1;
    swapped = false;
  }

}

//--------------------------------------------------------------------------

// Setter methods.

vector<double> BrancherEmitRF::setmPostVec() {
  mPostSav.clear();
  mPostSav.push_back(mRes);       // ma
  mPostSav.push_back(0.0);        // mj
  mPostSav.push_back(mFinal);     // mk
  mPostSav.push_back(mRecoilers); // mAK
  return mPostSav;
}

void BrancherEmitRF::setidPost() {
  idPostSav.clear();
  idPostSav = idSav;
  // Insert gluon in second position.
  idPostSav.insert(idPostSav.begin() + 1, 21);
}

void BrancherEmitRF::setStatPost() {
  statPostSav.resize(iSav.size() + 1, 52);
  statPostSav[posFinal] = 51;
  statPostSav[posFinal+1] = 51;
}

int BrancherEmitRF::iNew() {
  if (posFinal > 0 && iSav[posFinal] > 0
      && mothers2daughters.find(iSav[posFinal]) != mothers2daughters.end())
    return mothers2daughters[iSav[posFinal]].second;
  return 0;
}

void BrancherEmitRF::setMaps(int sizeOld) {
  mothers2daughters.clear();
  daughters2mothers.clear();
  posNewtoOld.clear();

  // For updating the children of existing parents.  Save children of
  // F (treat like 1->2 splitting).
  mothers2daughters[iSav[posFinal]] = make_pair(sizeOld, sizeOld + 1);
  daughters2mothers[sizeOld] = make_pair(iSav[posFinal], 0);
  daughters2mothers[sizeOld+1] = make_pair(iSav[posFinal], 0);

  //Save recoilers and insert the new emission at position 1.
  int iInsert = sizeOld + 2;
  unsigned int posNewEmit = 1;
  for (unsigned int pos = 0; pos < iSav.size(); pos++) {
    if (pos >= posNewEmit) posNewtoOld[pos + 1] = pos;
    else posNewtoOld[pos] = pos;
    if (pos == posRes || pos == posFinal) continue;
    else {
      mothers2daughters[iSav[pos]] = make_pair(iInsert, iInsert);
      daughters2mothers[iInsert] = make_pair(iSav[pos], iSav[pos]);
      iInsert++;
    }
  }
}

//--------------------------------------------------------------------------

// Generic method, assumes setter methods called earlier.

bool BrancherEmitRF::getNewParticles(Event& event, vector<Vec4> momIn,
  vector<int> hIn, vector<Particle> &pNew, Rndm* rndmPtr, Colour*) {

  // Initialize.
  unsigned int nPost = iSav.size() + 1;
  pNew.clear();
  setidPost();
  setStatPost();
  double scaleNew = sqrt(q2NewSav);
  setMaps(event.size());

  // Check everything set.
  if(momIn.size() != nPost || hIn.size() != nPost ||
     idPostSav.size() != nPost || statPostSav.size() != nPost) return false;

  // Generate new colour tag.
  int lastTag = event.lastColTag();
  int resTag = 0;
  int newTag = 0;
  if (colFlowRtoF) resTag = event[iSav[posRes]].col();
  else resTag = event[iSav[posRes]].acol();
  // New tag can't be same colour as neighbour.
  while (newTag%10 == resTag%10 || newTag%10 == 0)
    newTag = lastTag + 1 + rndmPtr->flat()*10;

  // Now populate particle vector.
  for (unsigned int ipart = 0; ipart < nPost; ++ipart) {
    Particle newPart;
    // Set mass and colours (we have repurposed mPost for antenna
    // function mass scales). This is new emission.
    if (posNewtoOld.find(ipart) == posNewtoOld.end()) {
      newPart.m(0.0);
      if (colFlowRtoF) newPart.cols(resTag, newTag);
      else newPart.cols(newTag, resTag);
    // Skip the resonance.
    } else if (posNewtoOld[ipart] == posRes) continue;
    else {
      newPart.m(mSav[posNewtoOld[ipart]]);
      int colNow  = event[iSav[posNewtoOld[ipart]]].col();
      int acolNow = event[iSav[posNewtoOld[ipart]]].acol();
      if (posNewtoOld[ipart] == posFinal) {
        if (colFlowRtoF) colNow = newTag;
        else acolNow = newTag;
      }
      newPart.cols(colNow,acolNow);
    }

    // Set other pre-determined particle properties.
    newPart.status(statPostSav[ipart]);
    newPart.id(idPostSav[ipart]);
    newPart.pol(hIn[ipart]);
    newPart.p(momIn[ipart]);
    newPart.setEvtPtr(&event);
    newPart.scale(scaleNew);
    newPart.daughters(0,0);
    if (abs(newPart.m() - newPart.mCalc()) > 0.001) return false;
    pNew.push_back(newPart);
  }
  colTagSav=newTag;
  return true;

}

//--------------------------------------------------------------------------

// Generate a new Q2 scale.

double BrancherEmitRF::genQ2(int evTypeIn, double Q2MaxNow, Rndm* rndmPtr,
  const EvolutionWindow* evWindowPtrIn, double colFac,
  vector<double> headroomIn, vector<double> enhanceIn,
  int verboseIn) {

  // Save headroom and enhancement factors.
  headroomSav = (headroomIn.size() >= 1) ?  headroomIn[0] : 1.0;
  enhanceSav  = (enhanceIn.size() >= 1) ?  enhanceIn[0] : 1.0;
  if (zetaIntSave <= 0.) {
    hasTrialSav = true;
    q2NewSav = 0.;
    return q2NewSav;
  }

  // pT evolution.
  if (evTypeIn == 1) {

    // Get multiplicative factors.
    double prefactor = kallenFacSav;
    prefactor *= colFac;
    prefactor *= headroomSav * enhanceSav;

    // Save info for later (to be used in pAccept).
    evTypeSav   = evTypeIn;
    q2BegSav    = Q2MaxNow;
    evWindowSav = evWindowPtrIn;
    colFacSav   = colFac;

    // Fixed alphaS.
    double logR = log(rndmPtr->flat());
    if (evWindowPtrIn->runMode <= 0){
      // Use max possible value for alphaS.
      prefactor *= evWindowPtrIn->alphaSmax;
      // Inverse of Q2 integral for fixed alphaS.
      q2NewSav = Q2MaxNow*exp(logR/(prefactor*zetaIntSave));
    // Running alphaS.
    } else{
      prefactor /= evWindowPtrIn->b0;
      double muRScaleMod = evWindowPtrIn->kMu2/evWindowPtrIn->lambda2;
      double logQ2Ratio = exp(logR/(prefactor*zetaIntSave));
      double logQ2maxFactor = log(Q2MaxNow*muRScaleMod);
      q2NewSav = exp(logQ2maxFactor*logQ2Ratio)/muRScaleMod;
    }
    if (q2NewSav > Q2MaxNow) {
      if (verboseIn >= superdebug)
        cout << "evolution mode = " << evWindowPtrIn->runMode << endl
             << "prefactor = " << prefactor << " zetaIntSave = " << zetaIntSave
             << " logR =  " << logR << endl
             << " kmu2 = " << evWindowPtrIn->kMu2
             << " lambda2 = " << evWindowPtrIn->lambda2 << endl;
      string errorMsg = "Error in "+__METHOD_NAME__
        +": Generated q2New > q2Max"+" Returning -1.";
      cout<<errorMsg<<endl;
      q2NewSav = -1.;
    }
  } else {
    if (verboseIn >= normal) {
      stringstream ss;
      ss << "evTypeIn = " << evTypeIn;
      string errorMsg = "Error in "+__METHOD_NAME__
        +": Unsupported Evolution Type."+" "+ss.str();
      cout<<errorMsg<<endl;
    }
    return 0.;
  }

  // Set flag that this call produces a saved trial.
  hasTrialSav = true;
  return q2NewSav;

}

//--------------------------------------------------------------------------

// Generate complementary invariant(s) for saved trial scale. Return
// false if no physical kinematics possible.

bool BrancherEmitRF::genInvariants(vector<double>& invariants,Rndm* rndmPtr,
  int verboseIn) {

  // Initialise and check we have a generated q2.
  invariants.clear();
  invariantsSav.clear();
  setmPostVec();
  if (q2NewSav <= 0.) return false;

  // Calculate invariants from zeta, q2.
  double zetaNext = getZetaNext(rndmPtr);
  double sjk = sAK*(zetaNext - 1.0);
  double saj = q2NewSav*(1.0 + sAK/sjk);
  double sak = sAK + sjk - saj;
  if (verboseIn >= louddebug) {
    stringstream ss;
    ss << "Phase space point: Q2next = " << q2NewSav << " zeta = " << zetaNext;
    printOut(__METHOD_NAME__, ss.str());
    ss.str("Scaled invariants: yaj = ");
    ss << saj/(sjk+sAK) << " yjk = " << sjk/(sjk+sAK);
    printOut(__METHOD_NAME__, ss.str());
  }

  //Save regardless.
  invariantsSav.push_back(sAK);
  invariantsSav.push_back(saj);
  invariantsSav.push_back(sjk);
  invariantsSav.push_back(sak);

  // Veto if the point is outside the available phase space.
  if (vetoPhSpPoint(saj, sjk, sak, verboseIn)) return false;
  else {invariants = invariantsSav; return true;}

}

//--------------------------------------------------------------------------

// Compute antPhys/antTrial, given antPhys.

double BrancherEmitRF::pAccept(const double antPhys, int verboseIn) {

  // Check q2.
  if (q2NewSav <= 0.) {
    if (verboseIn >= normal) {
      string errorMsg = "Error in "+__METHOD_NAME__+": q2NewSav not set."+
        " Returning 0.";
      cout<<errorMsg<<endl;
    }
    return 0.;
  }

  // Reduced trial antenna.
  double antTrial = 2.0/q2NewSav;

  // Multiply by headroom factor.
  antTrial *= headroomSav;

  // Multiply by col factor.
  antTrial *= colFacSav;

  //Multiply by alphaS, default is fixed alphaS.
  double alphaSTrial = evWindowSav->alphaSmax;
  if (evWindowSav->runMode >= 1) {
    double mu2 = q2NewSav;
    mu2 *= (evWindowSav->kMu2/evWindowSav->lambda2);
    alphaSTrial = 1.0/log(mu2)/evWindowSav->b0;
  }
  antTrial *= alphaSTrial;
  return antPhys/antTrial;

}

//--------------------------------------------------------------------------

// Protected helper methods for internal class use.

double BrancherEmitRF::KallenFunction(double x, double y, double z) {
  return x*x + y*y + z*z - 2.*(x*y + x*z + y*z);}

double BrancherEmitRF::zetaIntSingleLim(double zetaLim) {
  double x = zetaLim-1;  return x + log(x);}

double BrancherEmitRF::zetaIntegral(double zLow, double zHigh) {
  return zetaIntSingleLim(zHigh) - zetaIntSingleLim(zLow);}

double BrancherEmitRF::getsAK(double mA, double mK, double mAK) {
  return mA*mA +mK*mK - mAK*mAK;}

double BrancherEmitRF::zetaMinCalc(double mA, double mK, double mAK,
  double Q2cut) {
  return 1.0/(1.0 - Q2cut/(mA*mA -(mAK + mK)*(mAK + mK)));}

double BrancherEmitRF::zetaMaxCalc(double mA, double mK, double mAK) {
  return 1.0 + ((mA-mAK)*(mA-mAK) - mK*mK)/getsAK(mA, mK, mAK);}

double BrancherEmitRF::getZetaNext(Rndm* rndmPtr) {
  // Returns the solution for zeta to R =
  // I(zeta,zetamin)/I(zetamax,zetamin).
  double R = rndmPtr->flat();
  // I(zetamax,zetamin).
  double intZrange  = zetaIntegral(zetaMin, zetaMax);
  double intZMin    = zetaIntSingleLim(zetaMin);
  // exp(I(zeta)).
  double expIntZeta = exp(intZrange*R + intZMin);
  double lambWFact  = LambertW(expIntZeta);
  // Now invert to get zeta.
  return 1. + lambWFact;
}

double BrancherEmitRF::calcQ2Max(double mA, double mAK, double mK){
  double aM2 = (mA-mAK)*(mA-mAK) - mK*mK;
  double bM2 = mAK*(mA-mAK) + mK*mK;
  double cM = mA-mAK;
  return aM2*aM2*mA/(2.0*cM*bM2);
}

//--------------------------------------------------------------------------

// Veto point if outside available phase space.

bool BrancherEmitRF::vetoPhSpPoint(double saj, double sjk, double sak,
  int verboseIn) {

  // Make copies of masses (just for compactness of notation).
  double mAK = mRecoilers;
  double ma  = mPostSav[0];
  double mj  = mPostSav[1];
  double mk  = mPostSav[2];

  // Common sense: saj, sjk > 0. Not an error for splitters - mass
  // effects can make negative and push outside generated phase space.
  if (saj<0. || sjk<0.) {
    if (verboseIn >= louddebug) {
      stringstream ss;
      ss << "Negative invariants. saj = " << saj << " sjk = " << sjk;
      printOut(__METHOD_NAME__, ss.str());
    }
    return true;
  }

  // On-shell X condition.
  double invDiff = ma*ma + mj*mj + mk*mk - saj - sak + sjk - mAK*mAK;
  if (invDiff > 0.001) {
    if (verboseIn >= louddebug)
      printOut(__METHOD_NAME__, "Failed on-shell AK condition.");
    return true;
  }

  // On-shell j,k conditions.
  double Ek = sak/(2.0*ma);
  if (Ek*Ek < mk*mk) {
    if (verboseIn >= louddebug)
      printOut(__METHOD_NAME__, "Failed on-shell k condition.");
    return true;
  }
  double Ej = saj/(2.0*ma);
  if(Ej*Ej < mj*mj) {
    if (verboseIn >= louddebug)
      printOut(__METHOD_NAME__, "Failed on-shell j condition.");
    return true;
  }

  // When |cosTheta| < 1.
  double cosTheta = getCosTheta(Ej,Ek,mj,mk,sjk);
  if (abs(cosTheta) > 1.0) {
    if (verboseIn >= louddebug)
      printOut(__METHOD_NAME__, "Failed cos theta condition.");
    return true;
  }

  // This condition may be sufficient to remove above conditions.
  double det = saj*sjk*sak - saj*saj*mk*mk - sjk*sjk*ma*ma - sak*sak*mj*mj
    + 4.0*ma*ma*mj*mj*mk*mk;
  if (det <= 0.) {
    if (verboseIn >= louddebug)
      printOut(__METHOD_NAME__, "Gram det < 0 : Outside phase space");
  }
  return false;

}

//--------------------------------------------------------------------------

// Calculate maximum gluon energy in the centre of mass frame of res
// given cos theta.

double BrancherEmitRF::getEjMax(double cosTheta, double mA, double mAK,
  double mK) {

  double cos2Theta(cosTheta*cosTheta), sin2Theta(1. - cos2Theta),
    mA2(mA*mA), mK2(mK*mK), mAK2(mAK*mAK),
    tmp0(sqrt(sin2Theta*KallenFunction(mA2, mK2, mAK2) + 4.0*mAK2*mA2)),
    tmp1(mK/mA*tmp0), tmp2(cos2Theta*mK2 + mAK2), tmp3(mA2 - sin2Theta*mK2),
    tmp4((tmp2 + tmp1)/tmp3);
  double Emax = mA*(1 - tmp4)/2.0;
  double EabsoluteMax = mA/2.0 - (mK+mAK)*(mK+mAK)/(2.0*mA);
  return min(Emax,EabsoluteMax);

}

//==========================================================================

// BrancherSplitRF class for storing information on antennae between a
// coloured resonance and final state parton, and generating a new
// emission.

//--------------------------------------------------------------------------

void BrancherSplitRF::init(Event& event, vector<int> allIn,
  unsigned int posResIn, unsigned int posFIn, double Q2cut) {

  // Get Pythia indices of res and final.
  posRes     = posResIn;
  posFinal   = posFIn;
  int iRes   = allIn.at(posRes);
  int iFinal = allIn.at(posFinal);
  colFlowRtoF = event[iRes].col() == event[iFinal].col()
    && event[iRes].col() != 0;

  // Extract the momenta of the rest.
  Vec4 recoilVec(0., 0., 0., 0.);
  for (vector<int>::iterator pos=allIn.begin(); pos!=allIn.end(); ++pos) {
    if ((*pos == iRes) || (*pos == iFinal)) continue;
    recoilVec += event[*pos].p();
  }

  // This is not necesssarily p(res). In the case where one particle
  // always recieves the recoil, e.g. W in t->bWX it is p_t - p_X,
  Vec4 resVec = recoilVec + event[iFinal].p();
  mRes = resVec.mCalc();
  mFinal = 0.;
  mRecoilers = recoilVec.mCalc();
  sAK = getsAK(mRes, mFinal, mRecoilers);

  //Calculate common prefactor to trial integral.
  kallenFacSav = (0.5*sAK/(4.0*M_PI));
  kallenFacSav /= sqrt(KallenFunction(mRes*mRes, mFinal*mFinal,
      mRecoilers*mRecoilers));

  // Calculate zeta limits.
  zetaMin = zetaMinCalc(mRes, mFinal, mRecoilers,Q2cut);
  zetaMax = 1.0;

  // Calculate zeta integral (full phase space), for splitters is just flat.
  zetaIntSave= zetaMax-zetaMin;

  // Calculate Q2max.
  Q2MaxSav = calcQ2Max(mRes, mRecoilers,mFinal);
  branchType = 6;
  swapped = false;
  iAntSav = iXGsplitRF;

}

//--------------------------------------------------------------------------

// Setter methods.

vector<double> BrancherSplitRF::setmPostVec() {
  mPostSav.clear();
  mPostSav.push_back(mRes);
  mPostSav.push_back(mFlavSav);
  mPostSav.push_back(mFlavSav);
  mPostSav.push_back(mRecoilers);
  return mPostSav;
}

void BrancherSplitRF::setidPost(){
  idPostSav.clear();
  idPostSav = idSav;
  // Modify the splitting gluon to antiquark, insert quark in second position.
  if (colFlowRtoF) {
    idPostSav[posFinal] = -idFlavSav;
    idPostSav.insert(idPostSav.begin() + 1, idFlavSav);
  } else {
    idPostSav[posFinal] = idFlavSav;
    idPostSav.insert(idPostSav.begin() + 1, -idFlavSav);
  }
}

void BrancherSplitRF::setStatPost() {
  statPostSav.resize(iSav.size() + 1, 52);
  statPostSav[1] = 51;
  statPostSav[posFinal+1] = 51;
}

//--------------------------------------------------------------------------

// Generic method, assumes setter methods called earlier.

bool BrancherSplitRF::getNewParticles(Event& event, vector<Vec4> momIn,
  vector<int> hIn, vector<Particle>& pNew, Rndm*, Colour*){

  // Initialize.
  unsigned int nPost = iSav.size() + 1;
  pNew.clear();
  setidPost();
  setStatPost();
  double scaleNew = sqrt(q2NewSav);
  setMaps(event.size());
  if (momIn.size() != nPost || hIn.size() != nPost ||
      idPostSav.size() != nPost || statPostSav.size() != nPost ) return false;
  int resTag = 0;
  if (colFlowRtoF) resTag = event[iSav[posRes]].col();
  else resTag = event[iSav[posRes]].acol();

  // Now populate particle vector.
  for (unsigned int ipart = 0; ipart < nPost; ++ipart) {
    Particle newPart;
    // Set mass and colours, (we have repurposed mPost for antenna
    // function mass scales). This is new emission.
    if (posNewtoOld.find(ipart) == posNewtoOld.end()) {
      newPart.m(mFlavSav);
      if (colFlowRtoF) newPart.cols(resTag, 0);
      else newPart.cols(0, resTag);
    } else if (posNewtoOld[ipart] == posRes) {continue;
    } else {
      int colNow  = event[iSav[posNewtoOld[ipart]]].col();
      int acolNow = event[iSav[posNewtoOld[ipart]]].acol();
      if (posNewtoOld[ipart] == posFinal) {
        if (colFlowRtoF) colNow = 0;
        else acolNow = 0;
        newPart.m(mFlavSav);
      } else newPart.m(mSav[posNewtoOld[ipart]]);
      newPart.cols(colNow,acolNow);
    }

    //Set other pre-determined particle properties.
    newPart.status(statPostSav[ipart]);
    newPart.id(idPostSav[ipart]);
    newPart.pol(hIn[ipart]);
    newPart.p(momIn[ipart]);
    newPart.setEvtPtr(&event);
    newPart.scale(scaleNew);
    newPart.daughters(0,0);
    pNew.push_back(newPart);
  }
  colTagSav = 0;
  return true;

}

//--------------------------------------------------------------------------

// Generate a new Q2 scale.

double BrancherSplitRF::genQ2(int evTypeIn, double Q2MaxNow, Rndm* rndmPtr,
  const EvolutionWindow* evWindowPtrIn, double colFac,
  vector<double> headroomIn, vector<double> enhanceIn, int verboseIn) {

  if (zetaIntSave <= 0.) {
    hasTrialSav = true;
    q2NewSav = 0.;
    return q2NewSav;
  }

  // Total splitting weight summed over flavours
  double wtSum = 0.0;
  vector<double> wtFlav;
  unsigned int nFlav = headroomIn.size();
  if (nFlav != enhanceIn.size()) {
    if (verboseIn >= normal) {
      string errorMsg = "Error in "+__METHOD_NAME__
        +": Headroom and enhancement vectors have different sizes.";
      cout<<errorMsg<<endl;
    }
    return 0.;
  }
  for (unsigned int iFlav = 0; iFlav < nFlav; ++iFlav) {
    double wt = headroomIn[iFlav] * enhanceIn[iFlav];
    wtFlav.push_back(wt);
    wtSum += wt;
  }

  // pT evolution.
  if (evTypeIn == 1) {

    // Get multiplicative factors (sAK/(2*4pi*lambda^1/2)).
    double prefactor = kallenFacSav;
    prefactor  *= colFac;
    prefactor  *= wtSum;
    evTypeSav   = evTypeIn;
    q2BegSav    = Q2MaxNow;
    evWindowSav = evWindowPtrIn;
    colFacSav   = colFac;
    double logR = log(rndmPtr->flat());

    // Fixed alphaS.
    if (evWindowPtrIn->runMode <= 0) {
      // Use max possible value for alphaS.
      prefactor*= evWindowPtrIn->alphaSmax;
      // Inverse of Q2 integral for fixed alphaS.
      q2NewSav = Q2MaxNow*exp(logR/(prefactor*zetaIntSave));
    // Running alphaS.
    } else {
      prefactor /= evWindowPtrIn->b0;
      double muRScaleMod = evWindowPtrIn->kMu2/evWindowPtrIn->lambda2;
      double logQ2Ratio = exp(logR/(prefactor*zetaIntSave));
      double logQ2maxFactor = log(Q2MaxNow*muRScaleMod);
      q2NewSav = exp(logQ2maxFactor*logQ2Ratio)/muRScaleMod;
    }
  } else {
    if (verboseIn >= normal) {
      stringstream ss;
      ss << "evTypeIn = " << evTypeIn;
      string errorMsg = "Error in "+__METHOD_NAME__
        +": Unsupported Evolution Type."+" "+ss.str();
      cout<<errorMsg<<endl;
    }
    return 0.;
  }

  // Select flavour.
  double ranFlav = rndmPtr->flat() * wtSum;
  for (int iFlav = nFlav - 1; iFlav >= 0; --iFlav) {
    ranFlav -= wtFlav[iFlav];
    if (ranFlav < 0) {
      idFlavSav   = iFlav+1;
      mFlavSav    = evWindowSav->mass.at(idFlavSav);
      enhanceSav  = enhanceIn[iFlav];
      headroomSav = headroomIn[iFlav];
      break;
    }
  }

  // Debugging.
  if (verboseIn >= superdebug) {
    stringstream ss;
    ss << "Selected splitting flavour: " << idFlavSav;
    printOut(__METHOD_NAME__, ss.str());
  }
  if (q2NewSav > Q2MaxNow) {
    string errorMsg = "Error in "+__METHOD_NAME__
      +": Generated qq2New > q2Max"+" Returning -1.";
    cout<<errorMsg<<endl;
    q2NewSav = -1.;
  }
  hasTrialSav = true;
  return q2NewSav;

}

//--------------------------------------------------------------------------

// Generate complementary invariant(s) for saved trial scale. Return
// false if no physical kinematics possible.

bool BrancherSplitRF::genInvariants(vector<double>& invariants,Rndm* rndmPtr,
  int verboseIn) {

  //Initialize and check we have a generated q2.
  invariants.clear();
  invariantsSav.clear();
  if (q2NewSav <= 0.) return false;
  setmPostVec();

  // Get zeta.
  double zetaNext = getZetaNext(rndmPtr);
  if (zetaNext < 0.) cout << zetaMin<< "  " << zetaMax << endl;
  double m2q  = mFlavSav*mFlavSav;
  double sak  = zetaNext*sAK;
  double tmp1 = q2NewSav - (1.-zetaNext)*sAK + m2q;
  double tmp2 = sqrt(1.0 +4.0*sAK*q2NewSav/(tmp1*tmp1));
  double sjk  = tmp1*(1-tmp2)/2.0 -2.0*m2q ;
  double saj  = (1.0-zetaNext)*sAK + 2.0*m2q +sjk;
  if (verboseIn >= louddebug) {
    stringstream ss;
    ss << "Phase space point: Q2next = " << q2NewSav <<" zeta = " << zetaNext;
    printOut(__METHOD_NAME__, ss.str());
    ss.str("Scaled invariants: yaj = ");
    ss << saj/(sjk+sAK) << " yjk = " << sjk/(sjk+sAK);
    printOut(__METHOD_NAME__, ss.str());
  }

  //Save regardless.
  invariantsSav.push_back(sAK);
  invariantsSav.push_back(saj);
  invariantsSav.push_back(sjk);
  invariantsSav.push_back(sak);

  // Veto if the point is outside the available phase space.
  if (vetoPhSpPoint(saj, sjk, sak, verboseIn)) return false;
  else {invariants = invariantsSav; return true;}

}

//--------------------------------------------------------------------------

// Compute antPhys/antTrial, given antPhys. Note, antPhys should be
// normalised to include charge and coupling factors.

double BrancherSplitRF::pAccept(const double antPhys, int verboseIn) {

  if (q2NewSav <= 0.) {
    if (verboseIn >= normal) {
      string errorMsg = "Error in "+__METHOD_NAME__+": q2NewSav not set";
      cout<<errorMsg<<endl;
    }
    return 0.;
  } else if (invariantsSav.size() != 4) {
    if (verboseIn >= normal) {
      string errorMsg = "Error in "+__METHOD_NAME__+": invariants not set";
      cout<<errorMsg<<endl;
    }
    return 0.;
  }
  double saj = invariantsSav[1];
  double sjk = invariantsSav[2];
  double sak = invariantsSav[3];
  double m2q = mFlavSav*mFlavSav;
  double antTrial = 0.5/(sjk + 2.0*m2q);
  double norm = 1.0 + (sjk + 2.0*m2q)/(sAK+sjk+2.0*m2q)*(sak+m2q)/(saj-m2q);
  antTrial*=norm;

  // Multiply by headroom factor.
  antTrial *= headroomSav;

  // Multiply by col factor.
  antTrial *= colFacSav;

  // Multiply by alphaS, default is fixed alphaS.
  double alphaSTrial = evWindowSav->alphaSmax;

  // Running alphaS.
  if (evWindowSav->runMode >= 1) {
    double mu2 = q2NewSav;
    mu2 *= (evWindowSav->kMu2/evWindowSav->lambda2);
    alphaSTrial = 1.0/log(mu2)/evWindowSav->b0;
  }
  antTrial *= alphaSTrial;
  return antPhys/antTrial;

}

//==========================================================================

// The VinciaFSR class for resonant decays.

//--------------------------------------------------------------------------

// Initialize alphaStrong and related pTmin parameters (TimeShower).

void VinciaFSR::init( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn) {

  // Check if already initialized.
  if (isInit && settingsPtr->word("Merging:Process").compare("void") == 0)
    return;
  verbose = settingsPtr->mode("Vincia:verbose");
  if (verbose >= veryloud) printOut(__METHOD_NAME__, "begin --------------");
  allowforceQuit = false;
  forceQuit = false;
  nBranchQuit = -1;

  // Event counters and debugging. Set nAccepted to -1, so prepare
  // works properly.
  nAccepted       = -1;
  nSelected       = -1;
  nVetoUserHooks  = 0;
  nFailHadLevel   = 0;
  nCallPythiaNext = 0;

  // Showers on/off.
  doFF    = settingsPtr->flag("Vincia:doFF");
  doRF    = settingsPtr->flag("Vincia:doRF");
  doII    = settingsPtr->flag("Vincia:doII");
  doIF    = settingsPtr->flag("Vincia:doIF");
  doQED   = settingsPtr->flag("Vincia:doQED");

  // TODO: everything is evolved in PT in this version of VINCIA.
  evTypeEmit     = 1;
  evTypeSplit    = 1;

  // Store input pointers for future use.
  beamAPtr     = beamAPtrIn;
  beamBPtr     = beamBPtrIn;

  // Assume all events in same run have same beam-beam ECM.
  m2BeamsSav  = m2(beamAPtr->p(), beamBPtr->p());
  eCMBeamsSav = sqrt(m2BeamsSav);

  // Possibility to allow user veto of emission step.
  hasUserHooks    = (userHooksPtr != 0);
  canVetoEmission = (hasUserHooks && userHooksPtr->canVetoFSREmission());

  // Number of active quark flavours
  nGluonToQuark = settingsPtr->mode("Vincia:nGluonToQuark");

  // Number of flavours to be treated as massless (can be made
  // user-specifiable in future if desired).
  nFlavZeroMass = settingsPtr->mode("Vincia:nFlavZeroMass");

  // Global flag for helicity dependence.
  helicityShower = settingsPtr->flag("Vincia:helicityShower");

  // Global flag for sector showers on/off.
  sectorShower = settingsPtr->flag("Vincia:sectorShower");

  // Perturbative cutoff. Since emissions and splittings can have
  // different evolution measures, in principle allow for different
  // cutoff scales, for now forced same.
  q2CutoffEmit  = pow2(settingsPtr->parm("Vincia:cutoffScaleFF"));
  // Allow perturbative g->qq splittings to lower scales.
  q2CutoffSplit = pow2(settingsPtr->parm("Vincia:cutoffScaleFF"));

  // Set shower alphaS pointer.
  useCMW     = settingsPtr->flag("Vincia:useCMW");
  aSemitPtr  = &vinComPtr->alphaStrong;
  aSsplitPtr = &vinComPtr->alphaStrong;
  // Currently, CMW is applied to both emissions and splittings.
  if (useCMW) {
    aSemitPtr  = &vinComPtr->alphaStrongCMW;
    aSsplitPtr = &vinComPtr->alphaStrongCMW;
  }

  // AlphaS parameters.
  alphaSvalue    = settingsPtr->parm("Vincia:alphaSvalue");
  alphaSorder    = settingsPtr->mode("Vincia:alphaSorder");
  aSkMu2Emit     = settingsPtr->parm("Vincia:renormMultFacEmitF");
  aSkMu2Split    = settingsPtr->parm("Vincia:renormMultFacSplitF");
  alphaSmax      = settingsPtr->parm("Vincia:alphaSmax");
  alphaSmuFreeze = settingsPtr->parm("Vincia:alphaSmuFreeze");
  mu2freeze      = pow2(alphaSmuFreeze);

  // Smallest allowed scale for running alphaS.
  alphaSmuMin = 1.05 * max(aSemitPtr->Lambda3(), aSsplitPtr->Lambda3());
  mu2min      = pow2(alphaSmuMin);

  // For constant alphaS, set max = value (for efficiency).
  if (alphaSorder == 0) alphaSmax = alphaSvalue;
  initEvolutionWindows();

  // Settings for enhanced (biased) kernels.
  enhanceInHard   = settingsPtr->flag("Vincia:enhanceInHardProcess");
  enhanceInResDec = settingsPtr->flag("Vincia:enhanceInResonanceDecays");
  enhanceInMPI    = settingsPtr->flag("Vincia:enhanceInMPIshowers");
  enhanceAll      = settingsPtr->parm("Vincia:enhanceFacAll");
  // Explicitly allow heavy-quark enhancements only, not suppression
  enhanceBottom   = max(1., settingsPtr->parm("Vincia:enhanceFacBottom"));
  enhanceCharm    = max(1., settingsPtr->parm("Vincia:enhanceFacCharm"));
  enhanceCutoff   = settingsPtr->parm("Vincia:enhanceCutoff");

  // Resize pAccept to the maximum number of elements.
  pAccept.resize(max(weightsPtr->nWeights(), 1));

  // Statistics.
  nTrialsSum      = 0;
  int nAntPhysMax = 20;
  nTrials.resize(nAntPhysMax + 1);
  nTrialsAccepted.resize(nAntPhysMax + 1);
  nFailedVeto.resize(nAntPhysMax + 1);
  nFailedHull.resize(nAntPhysMax + 1);
  nFailedKine.resize(nAntPhysMax + 1);
  nFailedMass.resize(nAntPhysMax + 1);
  nFailedCutoff.resize(nAntPhysMax + 1);
  nClosePSforHQ.resize(nAntPhysMax + 1);
  nSectorReject.resize(nAntPhysMax + 1);

  // Initialize parameters for shower starting scale.
  pTmaxMatch     = settingsPtr->mode("Vincia:pTmaxMatch");
  pTmaxFudge     = settingsPtr->parm("Vincia:pTmaxFudge");
  pT2maxFudge    = pow2(pTmaxFudge);
  pT2maxFudgeMPI = pow2(settingsPtr->parm("Vincia:pTmaxFudgeMPI"));

  // Initialize the FSR antenna functions.
  if (verbose >= veryloud)
    printOut(__METHOD_NAME__, "initializing antennaSet");
  antSetPtr->init();
  kMapResEmit  = settingsPtr->mode("Vincia:kineMapRFemit");
  kMapResSplit = settingsPtr->mode("Vincia:kineMapRFsplit");

  // Initialise the QED shower module if not done already.
  if (!qedShowerPtr->isInit()) {
    if (verbose >= debug)
      printOut(__METHOD_NAME__, "initializing QED shower module");
    qedShowerPtr->init(beamAPtrIn, beamBPtrIn);
  }

  // Diagnostics.
  setDiagnostics(dynamic_pointer_cast<VinciaDiagnostics>(userHooksPtr));

  isInit=true;
  if (verbose >= veryloud) printOut(__METHOD_NAME__, "end --------------");

}

//--------------------------------------------------------------------------

// Possible limitation of first emission (TimeShower).

bool VinciaFSR::limitPTmax(Event& event, double, double) {

  // Check if limiting pT of first emission.
  if (pTmaxMatch == 1) return true;
  else if (pTmaxMatch == 2) return false;

  // Always restrict SoftQCD processes.
  else if (infoPtr->isNonDiffractive() || infoPtr->isDiffractiveA() ||
           infoPtr->isDiffractiveB() || infoPtr->isDiffractiveC())
    return true;

  // Look if jets or photons in final state of hard system (iSys = 0).
  else {
    const int iSysHard = 0;
    for (int i = 0; i < partonSystemsPtr->sizeOut(iSysHard); ++i) {
      int idAbs = event[partonSystemsPtr->getOut(iSysHard,i)].idAbs();
      if (idAbs <= 5 || idAbs == 21 || idAbs == 22) return true;
      else if (idAbs == 6 && nGluonToQuark == 6) return true;
    }
    // If no QCD/QED partons detected, allow to go to phase-space maximum
    return false;
  }

}

//--------------------------------------------------------------------------

// Top-level routine to do a full time-like shower in resonance decay
// (TimeShower).

int VinciaFSR::shower(int iBeg, int iEnd, Event& event, double pTmax,
  int nBranchMax) {

  // Verbose output.
  if (verbose >= debug) printOut("VinciaFSR::shower", "begin --------------");

  // Add new system, automatically with two empty beam slots.
  int iSys = partonSystemsPtr->addSys();

  // Verbose output.
  if (verbose >= 8) printOut("VinciaFSR::shower",
      "preparing to shower. System no. " + num2str(iSys));

  // Loop over allowed range to find all final-state particles.
  Vec4 pSum;
  for (int i = iBeg; i <= iEnd; ++i)
    if (event[i].isFinal()) {
      partonSystemsPtr->addOut( iSys, i);
      pSum += event[i].p();
    }
  partonSystemsPtr->setSHat( iSys, pSum.m2Calc() );

  // Now need to clear all systems of antennae. Should be fine
  // because shower() is only called in two cases: (1) showering off
  // hadronic resonances as they decay (e.g. upsilon) (2) by user, in
  // which case there should not be pre-existing systems.
  resEmitters.clear();
  resSplitters.clear();
  emitters.clear();
  splitters.clear();
  lookupBrancherRF.clear();
  lookupSplitterRF.clear();
  lookupBrancherFF.clear();
  lookupSplitter.clear();
  qedShowerPtr->iSystems.clear();
  qedShowerPtr->emitSystems.clear();
  qedShowerPtr->splitSystems.clear();
  qedShowerPtr->convSystems.clear();

  // Let prepare routine do the setup.
  prepare(iSys, event, false);

  // Begin evolution down in pT from hard pT scale.
  int nBranchNow = 0;
  do {
    // Do a final-state emission (if allowed).
    double pTtimes = pTnext(event, pTmax, 0.);
    if (pTtimes > 0.) {
      if (branch(event)) ++nBranchNow;
      pTmax = pTtimes;
    }

    // Keep on evolving until nothing is left to be done.
    else pTmax = 0.;
  } while (pTmax > 0. && (nBranchMax <= 0 || nBranchNow < nBranchMax));

  // Return number of emissions that were performed.
  return nBranchNow;

}

//--------------------------------------------------------------------------

// Method to add QED showers in hadron decays (TimeShower).

int VinciaFSR::showerQED(int iBeg, int iEnd, Event& event, double pTmax) {

  // Verbose output
  if(verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");

  // Construct a little QED system out of the given particles.
  partonSystemsPtr->addSys();
  int iSys = partonSystemsPtr->sizeSys()-1;
  // We could check if they all have the same mother and treat as
  // resonance decay, but currently do not.
  if (iBeg > iEnd) {
    partonSystemsPtr->addOut(iSys,iBeg);
    partonSystemsPtr->addOut(iSys,iEnd);
  } else {
    for (int i=iBeg; i<iEnd; ++i) partonSystemsPtr->addOut(iSys,i);
  }
  qedShowerPtr->prepare( iSys, event, true);
  double q2      = pow2(pTmax);
  double q2min   = qedShowerPtr->q2min();
  int nBranchNow = 0;
  while(true) {
    q2 = qedShowerPtr->generateTrialScale(event, q2);
    if (q2 < q2min) break;
    if (qedShowerPtr->checkVeto(event)) {
      qedShowerPtr->update(event, iSys);
      ++nBranchNow;
    }
  }
  return nBranchNow;

}

//--------------------------------------------------------------------------

// Method to add QED showers to partons below colour resolution scale
// (TimeShower).

int VinciaFSR::showerQEDafterRemnants(Event& event) {

  // Check if we are supposed to do anything.
  if (!doQED || infoPtr->getAbortPartonLevel()) return 0;
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");

  // Prepare for showering below hadronisation scale. Include partons
  // from all current systems (pass iSys = -1).
  qedShowerPtr->prepare( -1, event, true);
  // Retrieve actual iSys.
  int iSys       = partonSystemsPtr->sizeSys()-1;
  double q2      = qedShowerPtr->q2minColoured();
  double q2min   = max(qedShowerPtr->q2min(),1.e-13);
  int nBranchNow = 0;
  while(true) {
    q2 = qedShowerPtr->generateTrialScale(event, q2);
    if (q2 < q2min) break;
    if (qedShowerPtr->checkVeto(event)) {
      qedShowerPtr->update(event, iSys);
      ++nBranchNow;
    }
  }
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return nBranchNow;

}

//--------------------------------------------------------------------------

// Prepare system for evolution (TimeShower).

void VinciaFSR::prepare( int iSys, Event& event, bool){

  // Set isPrepared to false every time prepare is called.
  // Only reset back to true if method executes successfully.
  isPrepared = false;
  if (!isInit) return;

  // Check if we are supposed to do anything
  if (!(doFF || doRF)) return;
  if (infoPtr->getAbortPartonLevel()) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Received abort from PartonLevel().","Aborting.");
    return;
  }

  // Last chance to print header if not done already.
  if (!headerIsPrinted && verbose >= normal) header();
  if (verbose >= debug) {
    printOut(__METHOD_NAME__, "begin --------------");
    if (verbose >= louddebug) {
      stringstream ss;
      ss << "Preparing system " << iSys;
      printOut(__METHOD_NAME__, ss.str());
    }
    if (verbose >= superdebug) {
      event.list();
      partonSystemsPtr->list();
    }
  }

  // Statistics (zero info on qBranch scales).
  if (iSys < 4)
    for (int j = 0; j < 11; ++j) {
      qBranch[iSys][j] = 0.0;
      pTphys[iSys][j] = 0.0;
    }

  // Counter for accepted events.
  nAccepted = max(nAccepted, infoPtr->nAccepted());

  // First time prepare is called for an event.
  bool firstCallInEvent  = false;
  // getCounter(3), number of times a Pythia::next() call has been begun.
  int nCallPythiaNextNow = infoPtr->getCounter(3);
  if (nCallPythiaNextNow > nCallPythiaNext) {
    nCallPythiaNext  = nCallPythiaNextNow;
    nVetoUserHooks   = 0;
    nFailHadLevel    = 0;
    firstCallInEvent = true;
  }
  // Last event got accepted.
  long nSelectedNow = infoPtr->nSelected();
  if (nSelectedNow > nSelected) {
    nSelected        = nSelectedNow;
    nVetoUserHooks   = 0;
    nFailHadLevel    = 0;
    firstCallInEvent = true;
  }
  // getCounter(10), number of times the selection of a new hard
  // process has been begun. = 1 if no user hooks veto, > 1 if user
  // hooks veto.
  int nVetoUserHooksNow = (infoPtr->getCounter(10)-1);
  if (nVetoUserHooksNow > nVetoUserHooks) {
    nVetoUserHooks   = nVetoUserHooksNow;
    firstCallInEvent = true;
  }
  // getCounter(14), number of times loop over parton- and
  // hadron-level processing has begun for a hard process. = 1 if
  // everything after user hooks veto is succesful, > 1 eg if hadron
  // level fails.
  int nFailHadLevelNow = (infoPtr->getCounter(14) - 1);
  if (nFailHadLevelNow > nFailHadLevel) {
    nFailHadLevel    = nFailHadLevelNow;
    firstCallInEvent = true;
  }

  // Resetting for first time in new event.
  if (firstCallInEvent) {
    forceQuit = false;

    // Reset counters, weights in new events, and clear system information.
    vinComPtr->resetCounters();
    weightsPtr->resetWeights(infoPtr->nAccepted());
    clearContainers();
  }

  // Allow to quit after a certain number of emissions per event (just
  // for testing).
  if (forceQuit){
    if (verbose >= debug) printOut(__METHOD_NAME__, "User forced quit early.");
    return;
  }

  // Sanity check: at least two particles in system.
  int sizeSystem = partonSystemsPtr->sizeAll(iSys);
  if (sizeSystem <= 1) return;

  // Reset antenna list for first interaction and for resonance
  // decays. We don't have a starting scale for this system yet.
  Q2hat[iSys] = 0.0;
  // After prepare we always have zero branchings.
  nBranch[iSys] = 0;
  nBranchFSR[iSys] = 0;
  if (doDiagnostics) diagnosticsPtr->setnBranchSys(iSys,nBranch[iSys]);

  // Initialize polarisation flag (only consider final-state partons).
  bool checkIncoming = false;
  if (helicityShower)
    polarisedSys[iSys] = mecsPtr->isPolarised(iSys, event, checkIncoming);
  else polarisedSys[iSys] = false;
  stateChangeSys[iSys] = true;
  stateChangeLast      = true;
  iSysWin = iSys;
  iNewSav = 0;

  // Note that we have not yet binned the first branching scale.
  if (iSys == 0) firstQBranchBinned = false;

  // We are not creating new copies of the particles. Colour and
  // polarization information may be changed or added, respectively,
  // and masses may be set to zero for partons VINCIA wants to treat
  // as massless.
  bool makeNewCopies = false;

  // Note, for 2->2 systems, ISR::prepare() is called before
  // FRS::prepare() (if doISR) so ISR may already have done
  // everything.
  if ((doIF || doII) && isrPtr->prepared(iSys)) {
    makeNewCopies = false;

    // Ensure consistency between ISR + FSR lists.
    isHardSys[iSys]      = isrPtr->isHardSys[iSys];
    isResonanceSys[iSys] = false;
    doMECsSys[iSys]      = isrPtr->doMECsSys[iSys];
    polarisedSys[iSys]   = isrPtr->polarisedSys[iSys];

  // If ISR::prepare() not called for this system, prepare it now.
  } else {

    // Since ISR::prepare() not called, this will normally be a resonance
    // decay system, but still no harm in checking explicitly.
    if (partonSystemsPtr->hasInAB(iSys)) {
      isHardSys[iSys]      = ( iSys == 0 );
      isResonanceSys[iSys] = false;
    } else {
      isHardSys[iSys]      = false;
      isResonanceSys[iSys] = partonSystemsPtr->hasInRes(iSys);
      int iAncestor        = event[partonSystemsPtr->getOut(iSys,0)].mother1();
      // Make sure iAncestor is recorded as iInA in the parton system
      // (not standard usage in PYTHIA but useful in VINCIA).
      partonSystemsPtr->setInA(iSys, iAncestor);
    }

    // Make light quarks (and initial-state partons) explicitly massless.
    if (!vinComPtr->mapToMassless(iSys, event, partonSystemsPtr,
        makeNewCopies)) return;
    // Then see if we know how to compute matrix elements for this conf.
    doMECsSys[iSys] = mecsPtr->prepare(iSys, event);
    // Then see if and whether we can assign helicities.
    if (doMECsSys[iSys] && helicityShower)
      polarisedSys[iSys] = mecsPtr->polarise(iSys, event);
    // Decide if we should be doing ME corrections for next order.
    if (doMECsSys[iSys]) doMECsSys[iSys] = mecsPtr->doMEC(iSys, 1);
    // Then see if we should colourise this conf.
    colourPtr->colourise(iSys, event);
  }

  // Find antennae.
  if (verbose >= superdebug) printOut(__METHOD_NAME__, "Finding antennae....");
  if(!getAntennae(iSys, event)) return;

  // Set starting scale for this system.
  setStartScale(iSys, event);

  // Let others know we got to the end.
  isPrepared = true;
  if (verbose >= debug) {
    if (verbose >= superdebug) list();
    printOut(__METHOD_NAME__, "end --------------");
  }

}

//--------------------------------------------------------------------------

// Update antenna list after each ISR emission (TimeShower).

void VinciaFSR::update( int iSys, Event& event, bool) {

  // Do nothing if not prepared for FSR.
  if (!isPrepared) return;
  if (verbose >= debug) {
    printOut(__METHOD_NAME__, "begin --------------");
    if (verbose >= louddebug) event.list();
  }

  // Update QED system.
  qedShowerPtr->update(event, iSys);
  if (isResonanceSys[iSys]) {
    if (verbose >=normal) infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Update called unexpectedly in resonance shower.","Exiting.");
    return;
  }

  // Count number of branches.
  nBranch[iSys]++;

  // Particles in the list are already updated by ISR. Find and save
  // all colours and anticolours; find all FF antennae.
  map<int,int> indexOfAcol;
  map<int,int> indexOfCol;
  vector< pair<int,int> > antFF;
  const bool findFF = true;
  const bool findIX = false;
  colourPtr->makeColourMaps(iSys, event, indexOfAcol, indexOfCol,
    antFF, findFF, findIX);

  // In principle, the colour maps could here be used to look for any
  // unmatched tags -> junctions.

  // Sanity check: can only shower QCD systems with more than 1 FF
  // connection.
  if (antFF.size() <= 0) return;

  // Update any final-state antennae with partons changed by ISR
  // branching.
  for (int i = 0; i < (int)emitters.size(); i++) {
    Brancher* brancherPtr = &emitters[i];
    // Update any antennae with legs that were modified by the ISR
    // branching, i.e. whose current legs have been marked with status
    // < 0.
    int i0Old = brancherPtr->i0();
    int i1Old = brancherPtr->i1();
    int i0New = i0Old;
    int i1New = i1Old;

    if (event[i0Old].status() < 0 || event[i1Old].status() < 0) {
      // Get new positions from indexOfCol, indexOfAcol (could also
      // use daughter information from old i0, i1).
      i0New = indexOfCol[event[i0Old].col()];
      i1New = indexOfAcol[event[i1Old].acol()];
      // Update emitter.
      brancherPtr->reset(brancherPtr->system(), event, i0New, i1New);

      // Update lookup map and erase old keys.
      pair<int,bool> key = make_pair(i0Old, true);
      if (lookupBrancherFF.find(key)!=lookupBrancherFF.end())
        lookupBrancherFF.erase(key);
      key = make_pair(i1Old, false);
      if (lookupBrancherFF.find(key)!=lookupBrancherFF.end())
        lookupBrancherFF.erase(key);
      // Add new keys.
      key = make_pair(i0New,true);
      lookupBrancherFF[key] = i;
      key = make_pair(i1New,false);
      lookupBrancherFF[key] = i;

      // Update splitters.
      if (event[i0Old].isGluon()) {
        if (event[i0New].isGluon())
          updateSplitter(event,i0Old,i1Old,i0New,i1New,true);
        else removeSplitter(i0Old);
      }
      if (event[i1Old].isGluon()) {
        if (event[i1New].isGluon())
          updateSplitter(event,i1Old,i0Old,i1New,i0New,false);
        else removeSplitter(i1Old);
      }
    }

    // Remove the antennae out of the list. This way we can check
    // later if ISR added a new FF antenna i0/i1 is colour/anticolour.
    pair<int,int> pairNow = make_pair(i0New,i1New);
    vector< pair<int,int> >::iterator iter;
    iter = find (antFF.begin(), antFF.end(), pairNow);
    if (iter != antFF.end()) antFF.erase(iter);
  }


  // Is there a FF connection left?
  for (int i = 0; i < (int)antFF.size(); i++) {
    int i0 = antFF[i].first;  // i0/iNew[0] is colour.
    int i1 = antFF[i].second; // i1/iNew[2] is anticolour.
    // Don't include II or IF antennae.
    if (!event[i0].isFinal() || !event[i1].isFinal()) continue;
    if (verbose >= superdebug) {
      stringstream ss;
      ss << "Creating antenna between " << i0 << " , " << i1
         << " col = " << event[i0].col();
      printOut(__METHOD_NAME__, ss.str());
    }
    // Store new trial QCD gluon emission antenna.
    saveEmitter(iSys, event, i0, i1);
    // Store new trial QCD gluon splitting antenna(e).
    if (event[i0].isGluon()) saveSplitter(iSys, event, i0, i1, true);
    if (event[i1].isGluon()) saveSplitter(iSys, event, i1, i0, false);
  }

  // Sanity check.
  if (emitters.size() + splitters.size() <= 0) {
    if (verbose >= quiteloud)
      printOut(__METHOD_NAME__, "WARNING: Did not find any QCD antennae.");
    return;
  }
  if (!check(event)) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": failed update antennae.");
    list();
    if (verbose >= superdebug) printLookup();
    infoPtr->setAbortPartonLevel(true);
    return;
  }
  if (verbose >=debug) {
    if (verbose >= louddebug) list();
    if (verbose >= superdebug) printLookup();
    printOut(__METHOD_NAME__, "end --------------");
  }

}

//--------------------------------------------------------------------------

// Select next pT in downwards evolution (TimeShower).

double VinciaFSR::pTnext(Event& event, double pTevolBegAll,
  double pTevolEndAll, bool, bool) {

  // Check if we are supposed to do anything.
  if (infoPtr->getAbortPartonLevel() || !isPrepared) return 0.;
  if (forceQuit) {
    if (verbose >= debug) printOut(__METHOD_NAME__, "User forced quit early.");
    return 0.;
  }
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");

  // Denote VINCIA scales by "q", PYTHIA ones by "pTevol".
  double q2Begin  = pow2(pTevolBegAll);
  double q2EndAll = pow2(pTevolEndAll);

  // Initialize.
  q2WinSav  = 0.;
  winnerQED = false;
  winnerPtr = 0;

  // Generate next gluon-emission trial scale (above qEndAll).
  if (doFF && emitters.size() > 0) {
    if ( !q2NextEmit(q2Begin, q2EndAll) ) return 0.;
  }

  // Generate next gluon-splitting trial scale and compare to current qWin.
  if (doFF && splitters.size() > 0) {
    if ( !q2NextSplit(q2Begin, q2EndAll) ) return 0.;
  }

  // Generate next resonance gluon-emission trial and compare to current qWin.
  if (doRF && resEmitters.size() > 0) {
    if ( !q2NextResEmit(q2Begin, q2EndAll) ) return 0.;
  }

  // Generate nex resonance gluon-splitting trial and compare to current qWin.
  if (doRF && resSplitters.size() > 0) {
    if ( !q2NextResSplit(q2Begin, q2EndAll) ) return 0.;
  }

  // Generate next QED trial scale and compare to current qWin.
  if (doQED) {
    double q2QED = qedShowerPtr->generateTrialScale(event, q2Begin);
    if (q2QED >q2Begin) {
      stringstream ss;
      ss << "q2Begin = "<<q2Begin<<" q2QED = " << q2QED;
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Genereated q2QED > q2Begin.",ss.str());
      return 0.;
    }
    // Check for winning condition.
    if (q2QED > q2WinSav && q2QED > q2EndAll) {
      q2WinSav  = q2QED;
      winnerQED = true;
      winnerPtr = NULL;
    }
  }

  // If non-zero branching scale found: continue.
  if (q2WinSav > q2EndAll) {
    if (verbose >= debug) {
      stringstream ss;
      if (winnerPtr != 0)
        ss << " QCD Winner at scale qWinNow = "
           <<  sqrt(q2WinSav)
           << " col = " << event[winnerPtr->i0()].col()
           << " in System " << winnerPtr->system()
           << " qbegin = "<< pTevolBegAll;
      else
        ss << "=== QED Winner at scale qWinNow = "
           << sqrt(q2WinSav);
      printOut(__METHOD_NAME__, ss.str());
    }
    if (verbose >= superdebug) list();

  // No more branchings. Finalize.
  } else {
    q2WinSav = 0.0;
    if (verbose>=superdebug) event.list();

    // TODO: add back.
    // Need to make sure this is only done once per event.
    // FSR is done; set the weights.
    // weightsPtr->doWeighting();

    // If we have not yet binned a branching scale, bin 0 as the first
    // branching scale.
    if (verbose >= normal && !firstQBranchBinned) {
      if (vinciaHistos.find("1stBranchingQE/eCM") != vinciaHistos.end())
        vinciaHistos["1stBranchingQE/eCM"].fill(0);
      firstQBranchBinned = true;
    }
  }
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return sqrt(q2WinSav);

}

//--------------------------------------------------------------------------

// Branch event, including accept/reject veto (TimeShower).

bool VinciaFSR::branch(Event& event, bool ){

  // Verbose output.
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");

  // Hand off QED branchings to QED brancher.
  if (winnerQED) return branchQED(event);

  // Now handle QCD branchings.
  iSysWin                 = winnerPtr->system();
  stateChangeLast         = false;
  stateChangeSys[iSysWin] = false;
  iNewSav = 0;

  // Mark this trial as used so we do not risk reusing it.
  winnerPtr->needsNewTrial();
  // Find out which branching type we are doing.
  iAntWin   = winnerPtr->iAntPhys();

  // Count up global number of attempted trials.
  ++nTrialsSum;
  nTrials[iAntWin]++;

  // Decide whether to accept the trial. Store new particles in pNew
  // if keeping.
  if (!acceptTrial(event)) {
    if (verbose >= debug)
      printOut(__METHOD_NAME__, "Trial rejected (failed acceptTrial)");
    return false;
  }

  // Update event record, add new daughters. Make a copy of the event
  // to update (may want to veto)! Make a copy of junction info.
  Event newevent = event;
  resJunctionInfo junctionInfoCopy;
  if (hasResJunction[iSysWin]) junctionInfoCopy=junctionInfo[iSysWin];
  if (!updateEvent(newevent,junctionInfoCopy)) {
    if (verbose >= loud)
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Failed to update event.");
    return false;
  }

  // Possibility to allow user veto of emission step.
  if (canVetoEmission) {
    if (userHooksPtr->doVetoFSREmission(event.size(), newevent,
        iSysWin,isResonanceSys[iSysWin])) {
      if (verbose >= debug) printOut(__METHOD_NAME__, "Trial rejected "
          "(failed UserHooks::doVetoFSREmission)");
      return false;
    }
  }
  if (doDiagnostics) diagnosticsPtr->checkEvent(iSysWin,newevent,event.size());

  // Everything accepted -> overwrite original event.
  event = newevent;
  if (hasResJunction[iSysWin]) junctionInfo[iSysWin] = junctionInfoCopy;

  // Update partonSystems.
  updatePartonSystems();

  // Check momentum conservation.
  if (!vinComPtr->checkCoM(iSysWin,event,partonSystemsPtr)) {
    infoPtr->setAbortPartonLevel(true);
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Failed momentum conservation test.");
    return false;
  }

  // Update antennae.
  if(!updateAntennae(event)) {
    if (verbose >= loud)
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Failed to update antennae");
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // If diagnostic histograms are on, write in the branching scale.
  if (verbose >= normal && !firstQBranchBinned) {
    if (vinciaHistos.find("1stBranchingQE/eCM") != vinciaHistos.end())
      vinciaHistos["1stBranchingQE/eCM"].fill(sqrt(q2WinSav)/event[0].e());
    firstQBranchBinned = true;
  }

  // Count the number of branchings in the system.
  nBranch[iSysWin]++;
  nBranchFSR[iSysWin]++;

  // Do user-defined diagnostics.
  if(doDiagnostics) diagnosticsPtr->setnBranchSys(iSysWin,nBranch[iSysWin]);

  // Check the event after each branching.
  if (!vinComPtr->showerChecks(event, false)) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": Failed shower checks.");
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // Statistics for first 10 branchings in first 4 systems
  if (iSysWin < 4 && nBranchFSR[iSysWin] < 11 && nBranchFSR[iSysWin]>0) {
    qBranch[iSysWin][nBranchFSR[iSysWin]]=sqrt(q2WinSav);
    if(iNewSav>0){
      pTphys[iSysWin][nBranchFSR[iSysWin]]=event[iNewSav].pT();
    }
  }

  // Force stop by user (debugging only).
  if (allowforceQuit) {
    if (nBranchFSR[iSysWin] >= nBranchQuit && nBranchQuit > 0) {
      forceQuit = true;
      if (verbose >= debug) {
        stringstream ss;
        ss << "User forced quit after " << nBranchQuit << " emissions.";
        printOut(__METHOD_NAME__, ss.str());
      }
    }
  }

  // Done
  nTrialsAccepted[iAntWin]++;
  stateChangeSys[iSysWin] = true;
  stateChangeLast         = true;
  pTLastAcceptedSav       = sqrt(q2WinSav);
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return true;

}

//--------------------------------------------------------------------------

// Utility to print antenna list; for debug mainly (TimeShower).

void VinciaFSR::list() const {
  // Loop over antenna list and print it.
  if (resEmitters.size() + resSplitters.size() +
      emitters.size() + splitters.size() == 0) {
    cout << " --------  The list of FF antennae is empty -------------------"
      "------------------------------------------\n";
    return;
  }
  cout << endl << endl;
  for (unsigned int i = 0; i < resEmitters.size(); ++i) {
    if (i == 0) resEmitters[i].list("Gluon Resonance Emission Antennae");
    else resEmitters[i].list();
  }
  for (unsigned int i = 0; i < resSplitters.size(); ++i) {
    if (i == 0) resSplitters[i].list("Gluon Resonance Splitting Antennae");
    else resSplitters[i].list();
  }
  for (int i = 0; i < (int)emitters.size(); ++i) {
    if (i == 0) emitters[i].list("Gluon Emission Antennae");
    else emitters[i].list();
  }
  for (int i = 0; i < (int)splitters.size(); ++i) {
    if (i == 0) splitters[i].list("Gluon Splitting Antennae");
    else splitters[i].list();
  }
  cout << " --------  End VINCIA FF Antenna Listing ------------------------"
    "----------------------------------\n";

}

//--------------------------------------------------------------------------

// Initialise pointers to Vincia objects.

void VinciaFSR::initVinciaPtrs(Colour* colourPtrIn,
  shared_ptr<VinciaISR> isrPtrIn, QEDShower* qedPtrIn,MECs* mecsPtrIn,
  Resolution* resolutionPtrIn, VinciaCommon* vinComPtrIn,
  VinciaWeights* vinWeightsPtrIn) {
  colourPtr     = colourPtrIn;
  isrPtr        = isrPtrIn;
  qedShowerPtr  = qedPtrIn;
  mecsPtr       = mecsPtrIn;
  resolutionPtr = resolutionPtrIn;
  vinComPtr     = vinComPtrIn;
  weightsPtr    = vinWeightsPtrIn;
}

//--------------------------------------------------------------------------

// Print header information (version, settings, parameters, etc.).

void VinciaFSR::header() {

  // Must be initialised before printing header.
  if (!isInit) return;

  // Avoid printing header several times.
  if (headerIsPrinted) return;
  headerIsPrinted = true;

  cout <<setprecision(3);
  cout.setf(ios::left);
  cout << "\n";
  cout << " *-------  VINCIA Global Initialization  ------"
       << "-------------------------------------------------*\n";

  // Print header information about shower.
  cout << " |\n";
  cout << " | QCD Shower:     doII,doIF,doFF,doRF       =   "
       << bool2str(doII,3) <<","<<bool2str(doIF,3)
       <<","<<bool2str(doFF,3)<<","<<bool2str(doRF,3)
       <<"\n";
  cout << " |                 nGluonToQuark (FSR)       = "
       << num2str(settingsPtr->mode("Vincia:nGluonToQuark"),9)<<"\n";
  cout << " |                 convertGluonToQuark (ISR) = "
       << bool2str(settingsPtr->flag("Vincia:convertGluonToQuark"),9)<<"\n";
  cout << " |                 convertQuarkToGluon (ISR) = "
       << bool2str(settingsPtr->flag("Vincia:convertQuarkToGluon"),9)<<"\n";
  cout << " |                 helicityShower            = "
       << bool2str(settingsPtr->flag("Vincia:helicityShower"),9)<<"\n";
  cout << " |                 sectorShower              = "
       << bool2str(settingsPtr->flag("Vincia:sectorShower"),9)<<"\n";

  // Print header information about alphaS
  cout << " |\n"
       << " | Alpha_s:        alphaS(mZ)|MSbar          = "
       << num2str(alphaSvalue,9)<<"\n"
       << " |                 order                     = "
       << num2str(alphaSorder,9)<<"\n";
  if (alphaSorder >= 1) {
    if (useCMW) {
      cout << " |                 LambdaQCD[nF]|MSbar       = "
           << num2str(vinComPtr->alphaStrong.Lambda3(),9)<<"[3] "
           << num2str(vinComPtr->alphaStrong.Lambda4(),7)<<"[4] "
           << num2str(vinComPtr->alphaStrong.Lambda5(),7)<<"[5] "
           << num2str(vinComPtr->alphaStrong.Lambda6(),7)<<"[6]\n";
      cout << " |                 LambdaQCD[nF]|CMW         = "
           << num2str(vinComPtr->alphaStrongCMW.Lambda3(),9)<<"[3] "
           << num2str(vinComPtr->alphaStrongCMW.Lambda4(),7)<<"[4] "
           << num2str(vinComPtr->alphaStrongCMW.Lambda5(),7)<<"[5] "
           << num2str(vinComPtr->alphaStrongCMW.Lambda6(),7)<<"[6]\n";
    } else {
      cout << " |                 LambdaQCD[nF]            = "
           << num2str(vinComPtr->alphaStrong.Lambda3(),9)<<"[3] "
           << num2str(vinComPtr->alphaStrong.Lambda4(),7)<<"[4] "
           << num2str(vinComPtr->alphaStrong.Lambda5(),7)<<"[5] "
           << num2str(vinComPtr->alphaStrong.Lambda6(),7)<<"[6]\n";
    }
    cout << " |                 useCMW                    = "
         << bool2str(settingsPtr->flag("Vincia:useCMW"),9)<<"\n";
    cout << " |                 renormMultFacEmitF        = "
         << num2str(settingsPtr->parm("Vincia:renormMultFacEmitF"),9)
         <<" (muR prefactor for FSR emissions)\n";
    cout << " |                 renormMultFacSplitF       = "
         << num2str(settingsPtr->parm("Vincia:renormMultFacSplitF"),9)
         <<" (muR prefactor for FSR splittings)\n";
    cout << " |                 renormMultFacEmitI        = "
         << num2str(settingsPtr->parm("Vincia:renormMultFacEmitI"),9)
         <<" (muR prefactor for ISR emissions)\n";
    cout << " |                 renormMultFacSplitI       = "
         << num2str(settingsPtr->parm("Vincia:renormMultFacSplitI"),9)
         <<" (muR prefactor for ISR splittings)\n";
    cout << " |                 renormMultFacConvI        = "
         << num2str(settingsPtr->parm("Vincia:renormMultFacConvI"),9)
         <<" (muR prefactor for ISR conversions)\n";

    cout << " |                 alphaSmuFreeze            = "
         << num2str(alphaSmuFreeze,9)<<"\n";
    cout << " |                 alphaSmax                 = "
         << num2str(alphaSmax,9)<<"\n";
  }

  // Print header information about IR regularization.
  cout << " |\n"
       << " |   IR Reg.:      cutoffScaleEmitFF         = "
       << num2str(sqrt(q2CutoffEmit),9)<<"\n"
       << " |                 cutoffScaleSplitFF        = "
       << num2str(sqrt(q2CutoffSplit),9)<<"\n"
       << " |                 cutoffScaleII             = "
       << num2str(settingsPtr->parm("Vincia:cutoffScaleII"),9)<<"\n"
       << " |                 cutoffScaleIF             = "
       << num2str(settingsPtr->parm("Vincia:cutoffScaleIF"),9)<<"\n";

  // Information about QED shower, so far main switches only.
  cout << " |\n"
       << " | QED Shower:     doQED                     = "
       << bool2str(doQED,9)<<endl;
  if (doQED) {
    cout << " |                 nGammaToQuark             = "
         <<num2str(settingsPtr->mode("Vincia:nGammaToQuark"),9)<<"\n"
         << " |                 nGammaToLepton            = "
         <<num2str(settingsPtr->mode("Vincia:nGammaToLepton"),9)<<"\n"
         << " |                 convertGammaToQuark       = "
         <<bool2str(settingsPtr->flag("Vincia:convertGammaToQuark"),9)<<"\n"
         << " |                 convertQuarkToGamma       = "
         <<bool2str(settingsPtr->flag("Vincia:convertQuarkToGamma"),9)<<"\n";
  }

  // Print header information about antenna functions.
  if (verbose >= 2) {
    cout<<" |\n"
        <<" | AntennaFunctions:         "
        <<"                      chargeFactor   kineMap"<<endl;
    vector<int> iAnt = antSetPtr->getIant();
    int modeSLC      = settingsPtr->mode("Vincia:modeSLC");

    // FF and RF antennae.
    for (int i=0;i<(int)iAnt.size();++i) {
      int iAntPhys = iAnt[i];
      AntennaFunction* antPtr = antSetPtr->getAnt(iAntPhys);
      if (antPtr == nullptr) continue;
      // Print antenna name.
      cout.setf(ios::left);
      cout << setprecision(2);
      string antName = antPtr->vinciaName()+" ["+antPtr->humanName()+"]";
      cout << " |                 " << left << setw(32) << antName << "    ";
      // Print colour/charge factor.
      double chargeFac = antPtr->chargeFac();
      cout<<fixed<<setw(6)<<chargeFac;
      // Put asterisk next to QG colour factor if using -1/NC2 correction.
      if (modeSLC == 2) {
        if (antPtr->vinciaName() == "Vincia:QGEmitFF" ||
            antPtr->vinciaName() == "Vincia:GQEmitFF" ||
            antPtr->vinciaName() == "Vincia:QGEmitRF" ||
            antPtr->vinciaName() == "Vincia:GQEmitRF") cout << "*";
        else cout << " ";
      } else cout << " ";
      int kineMap = antPtr->kineMap();
      cout << "    " << right << setw(5) << kineMap << left << "\n";
    }

    // II and IF antennae.
    AntennaSetISR* antSetISRPtr = isrPtr->antSetPtr;
    if (antSetISRPtr != nullptr) {
      vector<int> iAntISR = antSetISRPtr->getIant();
      for (int i = 0; i < (int)iAntISR.size(); ++i) {
        int iAntPhys = iAntISR[i];
        AntennaFunctionIX* antPtr = antSetISRPtr->getAnt(iAntPhys);
        if (antPtr == nullptr) continue;
        // Print antenna name.
        cout.setf(ios::left);
        cout << setprecision(2) << " |                 " << left << setw(32)
             << antPtr->vinciaName() + " [" + antPtr->humanName() + "]"
             << "    ";
        // Print colour/charge factor.
        double chargeFac = antPtr->chargeFac();
        cout << fixed << setw(6) << chargeFac;
        if (modeSLC == 2) {
          if (antPtr->vinciaName() == "Vincia:QGEmitII" ||
              antPtr->vinciaName() == "Vincia:GQEmitII" ||
              antPtr->vinciaName() == "Vincia:QGEmitIF" ||
              antPtr->vinciaName() == "Vincia:GQEmitIF") cout<<"*";
          else cout << " ";
        } else cout << " ";
        int kineMap = antPtr->kineMap();
        cout << "    " << right << setw(5) << kineMap << left << "\n";
      }
      if (modeSLC == 2)
        cout << " |                 *: GQ antennae interpolate between "
             << "CA and 2CF (modeSLC = 2)" << endl;
    }
  }
  // Print header information about matrix-element Corrections.
  mecsPtr->header();

  // Print references.
  cout << " |\n";
  cout << " |-------------------------------------------"
       << "------------------------------------------\n |\n";
  cout << " | References :"<<endl;
  cout << " |    VINCIA     : Fischer, Prestel, Ritzmann, Skands, "
       << "EPJC76(2016)589, arXiv:1605.06142" << endl;
  cout << " |    PYTHIA 8   : Sjostrand et al., CPC191(2015)159, "
       << "arXiv:1410.3012" << endl;
  cout << " |\n *-------  End VINCIA Initialization  "
       << "----------------------------------------------------*\n\n";
  cout.setf(ios::right);

}

//--------------------------------------------------------------------------

// Print final statistics information.

void VinciaFSR::printInfo(bool pluginCall) {

  // Weight info.
  if (!isInit) return;
  int nTotWeightsInt         = weightsPtr->nTotWeights;
  int nNonunityInitialWeight = weightsPtr->nNonunityInitialWeight;
  int nNegativeInitialWeight = weightsPtr->nNegativeInitialWeight;
  int nNonunityWeight        = weightsPtr->nNonunityWeight;
  int nNegativeWeight        = weightsPtr->nNegativeWeight;
  double nTotWeights         = max(1.0,((double)(nTotWeightsInt)));
  vector<double> weightSum   = weightsPtr->weightSum;
  vector<double> weightMax   = weightsPtr->weightsMax;
  vector<double> weightMin   = weightsPtr->weightsMin;

  cout << "\n";
  if (pluginCall)
    cout << " *--------  VINCIA FSR and ISR Statistics  ----------------"
         <<"--------------------------------------------------------*\n";
  else
    cout << " *--------  VINCIA FSR Statistics  ------------------------"
         <<"--------------------------------------------------------*\n";
  cout << " |                                                           "
       <<"                                                      |\n";
  cout << " |                                                           "
       <<"             " << setw(40) << " " << " |\n";
  cout << " | Total Number of (accepted) events                           = "
       << num2str(nTotWeightsInt, 9) << setw(40) << " " << " |\n";
  cout << " | Number of events with nonunity initial weight               = "
       << ((nNonunityInitialWeight <= 0) ?
         "     none                                " :
         num2str(nNonunityInitialWeight, 9)
         + " <-- (INITIALLY) WEIGHTED EVENTS")
       << setw(8) << " " << " |\n";
  cout << " | Number of events with negative initial weight               = "
       << ((nNegativeInitialWeight <= 0) ? "     none"
         : num2str(nNegativeWeight, 9))
       << setw(40) << " " << " |\n";
  cout << " | Number of events with nonunity reweight                     = "
       << ((nNonunityWeight <= 0) ? "     none                      "
         : num2str(nNonunityWeight, 9)+" <-- REWEIGHTED EVENTS" )
       << setw(18) << " " << " |\n";
  cout << " | Number of events with negative reweight                     = "
       << ((nNegativeWeight <= 0) ? "     none"
         : num2str(min(nNegativeWeight,nNonunityWeight), 9))
       << setw(40) << " " << " |\n";
  cout << " |                                                            "
       << "            " << setw(40) << " " << " |\n";
  cout << " |                      weight(i)          Avg Wt   Avg Dev"
       << "  rms(dev)      kUnwt          Expected effUnw          |\n";
  cout << " | This run               i =     IsUnw       <w>     <w-1>   "
       << "             1/<w>        Max Wt  <w>/MaxWt   Min Wt |\n";
  cout << " |"<<setw(4)<<" "<<"User settings            0 ";
  cout << ((abs(1.0-weightSum[0]/nTotWeights) < TINY) ?
           "   yes " : "    no " );
  cout << num2str(weightSum[0]/nTotWeights,9) << " ";
  cout << num2str(weightSum[0]/nTotWeights-1.0,9) << " ";
  cout << setw(9) << "-" << "  ";
  cout << ((weightSum[0] != 0.0) ?
           num2str(nTotWeights/weightSum[0], 9) : num2str(0.0, 9));
  cout << num2str(weightMax[0],14)<<" ";
  cout << ((weightMax[0] != 0.0) ? num2str(weightSum[0]/nTotWeights/
            weightMax[0],9) : num2str(0.0,9)) << " ";
  cout << num2str(weightMin[0],9) << " ";
  cout << "|\n";
  if (pluginCall) {
    cout << " | - - - -  FSR only Statistics  - - - - - - - - - - - - "
         << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - |\n";
    cout << " |                                                       "
         << "                 " << setw(40) << " " << " |\n";
  }
  if (verbose >= 2 && nTrialsSum > 0) {
    cout << " |                                               ----------"
         << "----------- P(Reject) --------------------"
         << setw(15) << "|" << endl;
    cout << " | Trial Veto Rates:   nTrials     nAcc   eff    Zeta  kinMap"
         << "    Veto  Sector      IR    mMin    HQPS"
         << setw(15) << "|" << endl;
    for (int iAntPhys = 0; iAntPhys < (int)nTrials.size(); ++iAntPhys) {
      double nTot = nTrials[iAntPhys];
      if (nTot <= 0) continue;
      cout << " | " << setw(17);
      if (iAntPhys == iQQemitFF) cout << "QQemitFF";
      else if (iAntPhys == iQGemitFF)  cout << "QGemitFF";
      else if (iAntPhys == iGQemitFF)  cout << "GQemitFF";
      else if (iAntPhys == iGGemitFF)  cout << "GGemitFF";
      else if (iAntPhys == iGXsplitFF) cout << "GXsplitFF";
      else if (iAntPhys == iQQemitRF)  cout << "QQemitRF";
      else if (iAntPhys == iQGemitRF)  cout << "QGemitRF";
      else if (iAntPhys == iXGsplitRF) cout << "XGsplitRF";
      cout << fixed<<setprecision(2)
           << " "<< num2str(int(nTot+0.5), 9)
           << " "<< num2str(int(nTrialsAccepted[iAntPhys]+0.5), 8)
           << "  " << nTrialsAccepted[iAntPhys]/nTot
           << "    " << nFailedHull[iAntPhys]*1.0/max(1., nTot);
      nTot -= nFailedHull[iAntPhys];
      cout << "    " << nFailedKine[iAntPhys]*1.0/max(1., nTot);
      nTot -= nFailedKine[iAntPhys];
      cout << "    " << nFailedVeto[iAntPhys]*1.0/max(1., nTot);
      nTot -= nFailedVeto[iAntPhys];
      cout << "    " << nSectorReject[iAntPhys]*1.0/max(1., nTot);
      nTot -= nSectorReject[iAntPhys];
      cout << "    " << nFailedCutoff[iAntPhys]*1.0/max(1., nTot);
      nTot -= nFailedCutoff[iAntPhys];
      cout << "    " << nFailedMass[iAntPhys]*1.0/max(1., nTot);
      nTot -= nFailedMass[iAntPhys];
      cout << "    " << nClosePSforHQ[iAntPhys]*1.0/max(1., nTot);
      cout << setw(16) << "|\n";
    }
    cout << " |" << setw(113) << " " << "|\n";
  }
  if (!pluginCall)
    cout << " *--------  End VINCIA FSR Statistics --------------------"
         << "---------------------------------------------------------*\n\n";

}

//--------------------------------------------------------------------------

// Print internal and diagnostic histrograms.

void VinciaFSR::printHistos() {
  for (map<string,Hist>::iterator iH = vinciaHistos.begin();
       iH != vinciaHistos.end(); ++iH) {
    string Hname=iH->first;
    if (vinciaHistos[Hname].getEntries() >= 1)
      cout << Hname << vinciaHistos[Hname] << endl;
  }
}

//--------------------------------------------------------------------------

// Write internal and diagnostic histrograms to file.

void VinciaFSR::writeHistos(string fileName, string suffix) {
  for (map<string,Hist>::const_iterator iH = vinciaHistos.begin();
       iH != vinciaHistos.end(); ++iH) {
    string Hname=iH->first;
    if (vinciaHistos[Hname].getEntries() >= 1) {
      string file = sanitizeFileName(
        fileName + "-Hist-" + Hname + "." + suffix);
      cout << "Writing " << file << endl;
      iH->second.table(file);
    }
  }
}

//--------------------------------------------------------------------------

// Get number of branchings in a system (return -1 if no such
// system). If iSys < 0, sum over all.

int VinciaFSR::getNbranch(int iSys) {
  int n = 0;
  if (iSys < 0)
    for (int i = 0; i < (int)nBranchFSR.size(); ++i) n += nBranchFSR[i];
  else if (iSys < (int)nBranchFSR.size()) n = nBranchFSR[iSys];
  else n = -1;
  return n;
}

//--------------------------------------------------------------------------

// Check event.

bool VinciaFSR::check(Event &event) {
  stringstream ss;
  for (int i = 0; i < (int)emitters.size(); ++i) {
    if (!event[emitters[i].i0()].isFinal()) {
      if (verbose > normal){
        ss << "Emitter " << i
           << " i0 = " << emitters[i].i0() << " not final.";
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed to update emitter (not final).", ss.str());
      }
      return false;
    } else if (!event[emitters[i].i1()].isFinal()) {
      if (verbose > normal) {
        ss << "Emitter " << i
           << " i1 = " << emitters[i].i1() << " not final.";
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed to update emitter (not final).", ss.str());
      }
      return false;
    }
  }
  for (int i = 0; i < (int)splitters.size(); ++i) {
    if(!event[splitters[i].i0()].isFinal()){
      if (verbose > normal) {
        ss << "Splitter " << i
           << " i0 = " << splitters[i].i0() << " not final.";
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed to update splitter (not final).", ss.str());
      }
      return false;
    } else if (!event[splitters[i].i1()].isFinal()) {
      if (verbose > normal) {
        ss << "Splitter " << i
           << " i1 = " << splitters[i].i1() << " not final.";
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed to update splitter (not final).", ss.str());
      }
      return false;
    }
  }
  if (verbose > superdebug)
    printOut(__METHOD_NAME__, "Passed all checks on antennae.");
  return true;

}

//--------------------------------------------------------------------------

// Initialize evolution windows.

void VinciaFSR::initEvolutionWindows(void){

  evWindowsEmit.clear();
  evWindowsSplit.clear();
  EvolutionWindow window;
  window.alphaSmax = alphaSmax;
  window.mass[1]   = 0.;
  window.mass[2]   = 0.;
  window.mass[3]   = 0.;
  window.mass[4]   = (nFlavZeroMass >= 4) ? 0.0 : particleDataPtr->m0(4);
  window.mass[5]   = (nFlavZeroMass >= 5) ? 0.0 : particleDataPtr->m0(5);
  window.mass[6]   = (nFlavZeroMass >= 6) ? 0.0 : particleDataPtr->m0(6);

  for(int iWindow = 0; iWindow < 4; iWindow++) {
    //Get minimum boundaries of window.
    double qMinNowEmit=getQ2Window(iWindow, q2CutoffEmit);
    double qMinNowSplit=getQ2Window(iWindow ,q2CutoffSplit);

    // Lowest window, use constant trial alphaS for scales below charm mass.
    if (iWindow == 0) {
      window.b0        = 0.;
      window.lambda2   = 0.;
      window.kMu2      = 1.;
      window.runMode   = 0 ;
      window.qMin = qMinNowEmit;
      evWindowsEmit[qMinNowEmit]   = window;
      window.qMin = qMinNowSplit;
      evWindowsSplit[qMinNowSplit] = window;
    } else {
      // Emissions.
      window.runMode = alphaSorder;
      int nFnow = 5;
      if (qMinNowEmit < particleDataPtr->m0(4)) nFnow = 3;
      else if (qMinNowEmit < particleDataPtr->m0(5)) nFnow = 4;
      else if (qMinNowEmit >= particleDataPtr->m0(6)) nFnow = 6;
      window.b0 = (33.0 - 2.0*nFnow) / (12.0 * M_PI);
      double lambdaNow = getLambda(nFnow,aSemitPtr);
      window.lambda2 = (lambdaNow*lambdaNow);
      window.kMu2 = aSkMu2Emit;
      window.qMin = qMinNowEmit;
      evWindowsEmit[qMinNowEmit]=window;
      // Splittings.
      nFnow = 5;
      if (qMinNowSplit < particleDataPtr->m0(4)) nFnow = 3;
      else if (qMinNowSplit < particleDataPtr->m0(5)) nFnow = 4;
      else if (qMinNowSplit >= particleDataPtr->m0(6)) nFnow = 6;
      window.b0 = (33.0 - 2.0*nFnow) / (12.0 * M_PI);
      lambdaNow = getLambda(nFnow,aSsplitPtr);
      window.lambda2 = (lambdaNow*lambdaNow);
      window.kMu2 = aSkMu2Split;
      window.qMin = qMinNowSplit;
      evWindowsSplit[qMinNowSplit]=window;
    }
  }

}

//--------------------------------------------------------------------------

// Return window Q2.

double VinciaFSR::getQ2Window(int iWindow, double q2cutoff){
  double qMinNow = 0.;
  switch (iWindow) {
  case 0:
    // [cutoff, mc]
    qMinNow = min(sqrt(q2cutoff),particleDataPtr->m0(4));
    break;
  case 1:
    // [mc, mb] with 4-flavour running trial alphaS.
    qMinNow = max(1.0,particleDataPtr->m0(4));
    break;
  case 2:
    // [mb, mt] with 5-flavour running trial alphaS.
    qMinNow = max(3.0,particleDataPtr->m0(5));
    break;
  default:
    // [>mt] with 6-flavour running trial alphaS.
    qMinNow = max(100.0,particleDataPtr->m0(6));
    break;
  }
  return qMinNow;
}

//--------------------------------------------------------------------------

// Return Lambda value.

double VinciaFSR::getLambda(int nFin, AlphaStrong* aSptr) {
  if (nFin <= 3) return 0.;
  else if (nFin == 4) return aSptr->Lambda4();
  else if (nFin == 5) return aSptr->Lambda5();
  else return aSptr->Lambda6();
}

//--------------------------------------------------------------------------

// Method to return renormalisation-scale prefactor.

double VinciaFSR::getkMu2(bool isEmit){
  double kMu2 = 1.;
  if (isEmit) {
    kMu2 = aSkMu2Emit;
    bool muSoftCorr = false;
    if (useCMW && muSoftCorr) {
      // TODO: generalize.
      double xj =  winnerPtr->getXj();
      double muSoftInvPow = 4;
      double a  = 1./muSoftInvPow;
      kMu2      = pow(xj,a) * kMu2 + (1.-pow(xj,a));
    }
  } else kMu2 = aSkMu2Split;
  return kMu2;
}

//--------------------------------------------------------------------------

// Method to return renormalisation scale. Default scale is kMu *
// evolution scale.

double VinciaFSR::getMu2(bool isEmit) {
  double mu2 = winnerPtr->q2Trial();
  double kMu2 = getkMu2(isEmit);
  mu2 = max(mu2min, mu2freeze + mu2*kMu2);
  return mu2;
}

//--------------------------------------------------------------------------

// Reset (or clear) sizes of all containers.

void VinciaFSR::clearContainers() {
  headroomSav.clear();
  enhanceSav.clear();
  Q2hat.clear();
  isHardSys.clear();
  isResonanceSys.clear();
  doMECsSys.clear();
  polarisedSys.clear();
  stateChangeSys.clear();
  nBranch.clear();
  nBranchFSR.clear();
  mSystem.clear();
  nG.clear();
  nQ.clear();
  nLep.clear();
  nGam.clear();
}

//--------------------------------------------------------------------------

// Method to set up antennae, called in prepare.

bool VinciaFSR::getAntennae(int iSys, Event& event){

  // Sanity check.
  if (partonSystemsPtr == nullptr) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": partonSystems pointer is NULL!");
    return false;
  }
  // Reset antenna list for first interaction and for resonance decays
  // (don't do this if this is a new system due to MPI).
  if (iSys == 0 || partonSystemsPtr->hasInRes(iSys) ||
      !partonSystemsPtr->hasInAB(iSys)) {
    resEmitters.clear();
    resSplitters.clear();
    emitters.clear();
    splitters.clear();
    lookupBrancherRF.clear();
    lookupSplitterRF.clear();
    lookupBrancherFF.clear();
    lookupSplitter.clear();
  }
  // Check that iSys is a valid value.
  if (iSys > partonSystemsPtr->sizeSys()) return false;
  // Check that system contains some final state partons.
  if (partonSystemsPtr->sizeOut(iSys) == 0) return false;

  // Fetch index of resonance if any.
  int iMother     = -1;
  // Colour information of resonance (will remain negative if no resonance).
  int resCol      = -1;
  int resACol     = -1;
  int colPartner  = -1;
  int acolPartner = -1;
  if (partonSystemsPtr->hasInRes(iSys)) {
    iMother = partonSystemsPtr->getInRes(iSys);
    //Check that mother no longer is in system.
    if (event[iMother].status() > 0) return false;
    resCol = event[iMother].col();
    resACol = event[iMother].acol();
  }

  // Map of colour index to decay product (Pythia index).
  map<int, int> coltoDecID;
  // Map of anticolour index to decay product (Pythia index).
  map<int, int> aColtoDecID;
  // List of all decay products.
  vector<int> daughters;

  //Loop over members of current system and get colour information.
  Vec4 pSum;
  for (int iPart = 0; iPart < partonSystemsPtr->sizeOut(iSys); iPart++) {

    // Sum total FS momentum.
    int iOut = partonSystemsPtr->getOut(iSys, iPart);
    pSum += event[iOut].p();
    // Require final state.
    if (!event[iOut].isFinal()) return false;
    // Check colour.
    if (event[iOut].col() !=0 ) {
      // Check if colour partner of a resonance.
      if (event[iOut].col() == resCol) colPartner = iOut;
      // Otherwise save.
      else coltoDecID[event[iOut].col()] = iOut;
    }
    if (event[iOut].acol() !=0 ) {
      // Check if colour partner of a resonance.
      if (event[iOut].acol()==resACol) acolPartner = iOut;
      // Otherwise save.
      else aColtoDecID[event[iOut].acol()] = iOut;
    }
    // Save all.
    if (iOut != colPartner && iOut != acolPartner) daughters.push_back(iOut);
  }
  double mSys = m(pSum);

  // Check momentum conservation.
  Vec4 total(0., 0., 0., 0.);
  if (!isResonanceSys[iSys]) {
    if (partonSystemsPtr->getInA(iSys) > 0)
      total += event[partonSystemsPtr->getInA(iSys)].p();
    if (partonSystemsPtr->getInB(iSys) > 0)
      total += event[partonSystemsPtr->getInB(iSys)].p();
  } else total += event[partonSystemsPtr->getInRes(iSys)].p();
  total -= pSum;
  total /= mSys;
  if (abs(total.e()) > SMALL || abs(total.px()) > SMALL ||
      abs(total.py()) > SMALL || abs(total.pz()) > SMALL) {
    event.list();
    cout << "total = " << setprecision(10) << total.e() << " "
         << total.px() << " " << total.py() << " " <<total.pz() << endl;
    stringstream ss;
    ss << "iSys = " << iSys;
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Failed momentum conservation test.", ss.str());
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // Store total invariant mass of the final state.
  mSystem[iSys] = mSys;

  // Prepare QED shower (above colour resolution scale).
  if (doQED) qedShowerPtr->prepare(iSys, event, false);

  // Find any resonance antennae.
  if (colPartner > 0) {
    // Get a copy of daughters.
    vector<int> resSysAll = daughters;
    if (acolPartner != colPartner && acolPartner > 0)
      resSysAll.push_back(acolPartner);

    // Insert col partner and res at front (just convention).
    resSysAll.insert(resSysAll.begin(), colPartner);
    resSysAll.insert(resSysAll.begin(), iMother);
    unsigned int posRes(0), posPartner(1);
    saveResEmitter(iSys, event, resSysAll, posRes, posPartner, true);
    if (event[colPartner].isGluon())
      saveResSplitter(iSys, event, resSysAll, posRes, posPartner, true);
  }
  if (acolPartner > 0) {
    // Get a copy of daughters.
    vector<int> resSysAll = daughters;
    if (acolPartner != colPartner && colPartner > 0)
      resSysAll.push_back(colPartner);

    // Insert col partner and res at front (just convention).
    resSysAll.insert(resSysAll.begin(), acolPartner);
    resSysAll.insert(resSysAll.begin(), iMother);
    unsigned int posRes(0), posPartner(1);
    saveResEmitter(iSys, event, resSysAll, posRes, posPartner, false);
    if (event[acolPartner].isGluon())
      saveResSplitter(iSys, event, resSysAll, posRes, posPartner, false);
  }

  // Find any f-f that are colour connected, but not directly to a
  // resonance create normal branchers for these.
  for (map<int,int>::iterator it = coltoDecID.begin(); it != coltoDecID.end();
       ++it) {
    int col = it->first;
    int i0  = it->second;
    int i1  = aColtoDecID[col];

    // Exclude antennae that are not FF.
    if (!event[i0].isFinal() || !event[i1].isFinal()) continue;

    // Add to list of QCD gluon emission trial antennae.
    saveEmitter(iSys, event, i0, i1);

    // Add gluon-splitting antennae. Default, same 2->3 antenna
    // structure as for gluon emissions.
    if (event[i0].isGluon()) saveSplitter(iSys, event, i0, i1, true);
    if (event[i1].isGluon()) saveSplitter(iSys, event, i1, i0, false);
  }


  // Deal with any resonance junctions, n.b. assumes that these are
  // colour junctions not anticolour.
  hasResJunction[iSys] = false;
  if (isResonanceSys[iSys] && resCol > 0 && colPartner >0) {
    // Loop over junctions.
    for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
      // Loop over ends.
      for (int iLeg = 0; iLeg < 3; ++iLeg) {
        if (event.endColJunction(iJun, iLeg) == resCol) {
          // Found a resonance junction.
          hasResJunction[iSys] = true;
          junctionInfo[iSys].iJunction=iJun;
          junctionInfo[iSys].iEndCol=iLeg;
          junctionInfo[iSys].iEndColTag=resCol;
          junctionInfo[iSys].iEndQuark=colPartner;
          junctionInfo[iSys].colours.clear();
          junctionInfo[iSys].colours.push_back(resCol);

          // In context of matching might have many partons in systems
          // already.
          while (!event[junctionInfo[iSys].iEndQuark].isQuark()) {
            int colNow = event[junctionInfo[iSys].iEndQuark].acol();
            if (aColtoDecID.find(colNow) != aColtoDecID.end()) {
              int newPart = coltoDecID[colNow];
              junctionInfo[iSys].colours.push_back(colNow);
              junctionInfo[iSys].iEndQuark=newPart;
              junctionInfo[iSys].iEndColTag=colNow;
            } else {
              infoPtr->errorMsg("Error in "+__METHOD_NAME__
                +": Resonance involved in junction that cannot be traced.");
              hasResJunction[iSys] = false;
              break;
            }
          }
          if (event[junctionInfo[iSys].iEndQuark].col() == 0 ||
              !event[junctionInfo[iSys].iEndQuark].isFinal()) {
            infoPtr->errorMsg("Error in "+__METHOD_NAME__
              +": Failed to find end quark in resonance junction.");
            hasResJunction[iSys] = false;
            break;
          }
        }
      }
    }
  }

  // Count up number of gluons, quarks, and photons.
  nG[iSys]   = 0;
  nQ[iSys]   = 0;
  nGam[iSys] = 0;
  nLep[iSys] = 0;
  for (int i = 0; i < partonSystemsPtr->sizeAll(iSys); ++i) {
    Particle* partonPtr = &event[partonSystemsPtr->getAll(iSys, i)];
    if (partonPtr->id()==21) nG[iSys]++;
    else if (abs(partonPtr->id()) < 7) nQ[iSys]++;
    else if (abs(partonPtr->id()) == 22) nGam[iSys]++;
    else if (partonPtr->isLepton()) nLep[iSys]++;
  }

  // Sanity checks.
  if (verbose >= debug) {
    if (resEmitters.size() + resSplitters.size() +
        emitters.size() + splitters.size() <= 0)
      printOut(__METHOD_NAME__, "did not find any QCD antennae.");
    else if (verbose >= louddebug) {
      list();
      if (verbose >= superdebug) printLookup();
    }
  }
  return true;

}

//--------------------------------------------------------------------------

// Set starting scale of shower (power vs wimpy) for system iSys.

void VinciaFSR::setStartScale(int iSys, Event& event){

  // Set nIn: 1->n or 2->n.
  int nIn = 0;
  if (isResonanceSys[iSys]) nIn = 1;
  else if (partonSystemsPtr->hasInAB(iSys)) nIn = 2;

  // Set FSR starting scale of this system (can be different from qFac).
  // Resonance decay systems always start at Q2 = m2..
  if (isResonanceSys[iSys]) {
    Q2hat[iSys] = pow2(mSystem[iSys]);
    return;

  // Hard system: start at phase-space maximum or factorisation scale.
  } else if (isHardSys[iSys]) {
    if (verbose >= superdebug)
      printOut(__METHOD_NAME__, "Setting FSR starting scale for hard system");
    // pTmaxMatch = 1 : always start at QF (modulo kFudge).
    if (pTmaxMatch == 1) Q2hat[iSys] = pT2maxFudge * infoPtr->Q2Fac();
    // pTmaxMatch = 2 : always start at eCM.
    else if (pTmaxMatch == 2) Q2hat[iSys] = m2BeamsSav;
    // Else check if this event has final-state jets or photons.
    else {
      bool hasRad = false;
      for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
        int idAbs = event[partonSystemsPtr->getOut(iSys, i)].idAbs();
        if (idAbs <= 5 || idAbs == 21 || idAbs == 22) hasRad = true;
        if (idAbs == 6 && nGluonToQuark == 6) hasRad = true;
        if (hasRad) break;
      }
      // If no QCD/QED partons detected, allow to go to phase-space maximum.
      if (hasRad) Q2hat[iSys] = pT2maxFudge * infoPtr->Q2Fac();
      else Q2hat[iSys] = m2BeamsSav;
    }
  } else if (nIn == 2) {
    if (verbose >= superdebug)
      printOut(__METHOD_NAME__, "Setting FSR starting scale of MPI system");
    // Set starting scale for MPI systems: min of incoming parton
    // scales. Find positions of incoming colliding partons.
    int in1 = partonSystemsPtr->getInA(iSys);
    int in2 = partonSystemsPtr->getInB(iSys);
    Q2hat[iSys] = pT2maxFudgeMPI
      * pow2(min(event[in1].scale(),event[in2].scale()));
  } else {
    // Assume hadron -> partons decay. Starting scale = mSystem.
    Q2hat[iSys] = pow2(mSystem[iSys]);
  }

}

//--------------------------------------------------------------------------

// Auxiliary methods to generate trial scales for various shower
// components.

bool VinciaFSR::q2NextResEmit(const double q2Begin, const double q2End) {
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "begin --------------");
  double q2EndNow = max(q2End, q2CutoffEmit);
  bool gen = q2NextBranch<BrancherEmitRF>(resEmitters, evWindowsEmit,
    evTypeEmit, q2Begin, q2EndNow, true);
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "end --------------");
  return gen;
}

bool VinciaFSR::q2NextResSplit(const double q2Begin, const double q2End) {
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "begin --------------");
  double q2EndNow = max(q2End, q2CutoffSplit);
  bool gen = q2NextBranch<BrancherSplitRF>(resSplitters, evWindowsSplit,
    evTypeSplit, q2Begin, q2EndNow, false);
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__,"end --------------");
  return gen;
}

bool VinciaFSR::q2NextEmit(const double q2Begin, const double q2End) {
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "begin --------------");
  double q2EndNow = max(q2End, q2CutoffEmit);
  bool gen = q2NextBranch<BrancherEmitFF>(emitters, evWindowsEmit, evTypeEmit,
    q2Begin, q2EndNow, true);
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__,"end --------------");
  return gen;
}

bool VinciaFSR::q2NextSplit(const double q2Begin, const double q2End) {
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__,"begin --------------");
  double q2EndNow = max(q2End, q2CutoffSplit);
  bool gen = q2NextBranch<BrancherSplitFF>(splitters, evWindowsSplit,
    evTypeSplit, q2Begin, q2EndNow, false);
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__,"end --------------");
  return gen;
}

//--------------------------------------------------------------------------

// Return the Q2 for the next branching.

template <class Brancher> bool VinciaFSR::q2NextBranch(
  vector<Brancher>& brancherVec, const map<double, EvolutionWindow> &evWindows,
  const int evType, const double q2Begin, const double q2End, bool isEmit) {


  // Sanity check
  if (verbose  >= superdebug) {
    stringstream ss;
    ss << "qBegin = " << sqrt(q2Begin);
    printOut(__METHOD_NAME__, ss.str());
  }
  if (q2Begin <= q2End) {
    if (verbose >= louddebug)
      printOut(__METHOD_NAME__, "q2Begin below cutoff. Nothing to do.");
    return true;
  } else if (!isEmit && nGluonToQuark == 0) return true;

  // Loop over resonance antennae.
  unsigned int numAnt = brancherVec.size();
  if (verbose >= superdebug) {
    stringstream ss;
    ss << "Looping over " << numAnt << " antennae.";
    printOut(__METHOD_NAME__, ss.str());
  }
  unsigned int iAnt = 0;
  for (typename vector<Brancher>::iterator ibrancher = brancherVec.begin();
       ibrancher!=brancherVec.end(); ++ibrancher) {
    iAnt++;
    if (verbose >= superdebug) {
      stringstream ss;
      ss << "Antenna " << iAnt <<" / " << numAnt;
      printOut(__METHOD_NAME__,ss.str());
    }

    //Check if there is any phase space left for current antenna.
    double Q2MaxNow = min(q2Begin, ibrancher->getQ2Max(evType));
    if (Q2MaxNow < q2End) {
      if (verbose >= louddebug) printOut(__METHOD_NAME__,
          "No phase space left for current antenna, continuing.");
      continue;
    }

    // Check if a saved trial exists for this brancher.
    double q2Next = 0.;
    if (ibrancher->hasTrial()) {
      q2Next = ibrancher->q2Trial();
      if (verbose >= louddebug) {
        stringstream ss;
        ss << "Retrieving saved trial Q=" << sqrt(q2Next);
        printOut(__METHOD_NAME__, ss.str());
      }
    // Else generate new trial scale.
    } else {
      if (verbose >= superdebug)
        printOut(__METHOD_NAME__, "Generating new trial");

      // Fetch system and colour factor for current brancher.
      int iSys   = ibrancher->system();
      double colFac = getAnt(ibrancher->iAntPhys())->chargeFac();
      if (verbose >= louddebug) {
        stringstream ss;
        ss << "Starting shower for current brancher at Q=" << sqrt(Q2MaxNow);
        printOut(__METHOD_NAME__, ss.str());
      }

      // Impose evolution windows (for alphaS running); fetch the
      // current window.
      map<double, EvolutionWindow>::const_iterator
        it = evWindows.lower_bound(sqrt(Q2MaxNow));
      // Cast as a reverse iterator to go downwards in q2.
      map<double, EvolutionWindow>::const_reverse_iterator itWindowNow(it);

      // Go through regions.
      if (verbose >= superdebug)
        printOut(__METHOD_NAME__, "Looping over Q2 windows...");
      while(itWindowNow != evWindows.rend()) {

        // Bottom of current window.
        double Q2MinWindow = pow2(itWindowNow->first);
        const EvolutionWindow* windowPtr = &(itWindowNow->second);

        // Set headroom and enhancement factors.
        vector<double> headroomVec = getHeadroom(iSys, isEmit, Q2MaxNow);
        // For sector showers, use more headroom.
        if (sectorShower && isEmit) {
          if (ibrancher->colType0() == 2) headroomVec[0] *= 5;
          if (ibrancher->colType1() == 2) headroomVec[0] *= 5;
        }
        vector<double> enhanceVec = getEnhance(iSys, isEmit, Q2MaxNow);
        double Q2NextWindow = ibrancher->genQ2(evType, Q2MaxNow, rndmPtr,
          windowPtr, colFac, headroomVec, enhanceVec, verbose);
        if (Q2NextWindow < 0.) {
          infoPtr->setAbortPartonLevel(true);
          return false;
        }
        if (verbose >= superdebug) {
          stringstream ss;
          ss << "Generated QNextWindow = " << sqrt(Q2NextWindow)
             << " (QMinWindow = " << itWindowNow->first << " )";
          printOut(__METHOD_NAME__, ss.str());
        }

        // Check if Q2next is in the current window.
        if (Q2NextWindow > Q2MinWindow || Q2NextWindow <= 0.) {
          q2Next=Q2NextWindow;
          break;
        } else {
          if (verbose >= superdebug) printOut(__METHOD_NAME__,
              "QNext below window threshold. Continuing to next window.");
        }
        // Else go straight to next window.
        Q2MaxNow = Q2MinWindow;
        // Increment reverse iterator (go down in scale).
        itWindowNow++;
      } // End loop over evolution windows.
      if (verbose >= superdebug && itWindowNow == evWindows.rend())
        printOut(__METHOD_NAME__, "Out of windows. Continuing to "
          "next antenna.");
    } // End generate new trial for this antenna.

    // Check for winning condition.
    if (q2Next > q2WinSav && q2Next > q2End) {
      q2WinSav = q2Next;
      winnerPtr = &(*ibrancher);
    }
  } // End loop over QCD antennae.
  return true;

}

//--------------------------------------------------------------------------

// Perform a QED branching.

bool VinciaFSR::branchQED(Event& event) {

  // QED trial accept.
  if (verbose >= debug) printOut(__METHOD_NAME__, "begin --------------");
  int sizeOld = event.size();
  bool updated = false;
  if (qedShowerPtr->checkVeto(event)) {
    if (verbose >= debug)
      printOut(__METHOD_NAME__, "QED trial accepted. About to update.");
    qedShowerPtr->update(event, qedShowerPtr->sysWin());

    // Check momentum conservation.
    if (!vinComPtr->checkCoM(qedShowerPtr->sysWin(),event,partonSystemsPtr)) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Failed (E,p) conservation check.");
      infoPtr->setAbortPartonLevel(true);
      return false;
    }

    // Update QCD branchers after QED emissions.
    updated = updateAfterQED(event, sizeOld);

    // Check PartonSystems in debug mode.
    if (verbose > quiteloud) {
      if (partonSystemsPtr->hasInAB(iSysWin)) {
        int inA = partonSystemsPtr->getInA(iSysWin);
        int inB = partonSystemsPtr->getInB(iSysWin);
        if (inA <= 0 || inB <= 0 ) {
          stringstream ss;
          ss << "iSysWin = "<<iSysWin << " non-positive. inA = "<< inA
             << " inB = " << inB;
          infoPtr->errorMsg("Error in "+__METHOD_NAME__
            +": Non-positive incoming parton.", ss.str());
          infoPtr->setAbortPartonLevel(true);
          return false;
        } else if (event[inA].mother1() > 2 || event[inB].mother1() > 2) {
          stringstream ss;
          ss << "iSysWin = "<<iSysWin;
          infoPtr->errorMsg("Error in "+__METHOD_NAME__
            +": Failed to update incoming particles after QED branching.",
            ss.str());
          infoPtr->setAbortPartonLevel(true);
          return false;
        }
      }
    }

  // Else QED trial failed.
  } else {
    if (verbose >= debug) printOut(__METHOD_NAME__, "QED trial failed.");
    return false;
  }

  // Update saved scale of last branching.
  pTLastAcceptedSav = sqrt(qedShowerPtr->q2Trial);
  if (verbose >= debug) printOut(__METHOD_NAME__, "end --------------");
  return updated;

}

//--------------------------------------------------------------------------

// Perform an early antenna rejection.

bool VinciaFSR::rejectEarly(AntennaFunction* &antFunPtr, bool doMEC) {

  bool reject = true;
  if (winnerPtr->getBranchType() < 0) {
    if (verbose >= veryloud)
      printOut(__METHOD_NAME__, "WARNING: could not identify branching type.");
    return reject;
  }

  if (doDiagnostics) diagnosticsPtr->setBranchType(winnerPtr->getBranchType());

  // If enhancement was applied to the trial function but branching is
  // below enhancement cutoff, we do an early accept/reject here with
  // probability trial/enhanced-trial to get back to unenhanced trial
  // probability.

  // Trials only enhanced for enhanceFac > 1.
  if (winnerPtr->enhanceFac() > 1.0 &&
      winnerPtr->q2Trial() <= pow2(enhanceCutoff)) {
    if (rndmPtr->flat() > 1./winnerPtr->enhanceFac()) {
      if (verbose >= debug) printOut(__METHOD_NAME__,
          "Trial rejected (enhance applied below enhanceCutoff)");
      return reject;
    }
    // If passed, save that enhancement factor has now been canceled.
    winnerPtr->resetEnhanceFac(1.0);
  }

  // Generate post-branching invariants. Can check some vetos already
  // at this level, without full kinematics.
  vector<double> invariants;
  if (!winnerPtr->genInvariants(invariants, rndmPtr, verbose)) {
    if (verbose >= debug)
      printOut(__METHOD_NAME__, "Trial rejected (failed genInvariants)");
    if (doDiagnostics) diagnosticsPtr->checkInvariants(
        iSysWin,winnerPtr->iAntPhys(), winnerPtr->getInvariants(),false);
    ++nFailedHull[iAntWin];
    return reject;
  } else {
    if (doDiagnostics) diagnosticsPtr->checkInvariants(iSysWin,
        winnerPtr->iAntPhys(), invariants, true);
  }

  // Impose g->QQ mass thresholds for flavours treated as massless.
  if (iAntWin == iGXsplitFF && winnerPtr->idNew() <= nFlavZeroMass) {
    // m2(qq) > 4m2q => s(qq) > 2m2q, but allow for larger factor.
    // Factor 4 roughly matches n(g->bb) for massive b quarks.
    double facM = 4;
    if (invariants[1] < facM*pow2(particleDataPtr->m0(winnerPtr->idNew()))) {
      ++nFailedMass[iAntWin];
      return reject;
    }
  }

  // Compute physical antenna function (summed over possible helicities).
  double antPhys = getAntPhys(antFunPtr);
  // Get accept probability.
  pAccept[0]=pAcceptCalc(antPhys);
  // Do user diagnostics.
  if (doDiagnostics) diagnosticsPtr->checkpAccept(iSysWin, pAccept[0]);

  // If doing ME corrections, don't allow to reject yet.
  if (!doMEC) {
    // Check if rejecting the trial.
    double R = rndmPtr->flat();
    if (R > pAccept[0]) {
      // TODO: Note, here we want to put a call to something which computes
      // uncertainty variations for pure-shower branchings. We also
      // may want to take into account if there was a enhancement
      // applied to this branching.
      if (verbose >= debug)
        printOut(__METHOD_NAME__, "Trial rejected (failed R<pAccept)");
      ++nFailedVeto[iAntWin];
      return reject;
    }

    // Set accept probability to 1, so no later opportunity to reject
    // unles we apply an enhancement factor.
    pAccept[0] = 1.;
  }

  //Trial accepted so far, n.b. proper acccept/reject condition later.
  return false;

}

//--------------------------------------------------------------------------

// Compute physica antenna function.
double VinciaFSR::getAntPhys(AntennaFunction* &antFunPtr) {

  // Set antenna function pointer and check if this antenna is "on".
  antFunPtr= getAnt(iAntWin);
  if (antFunPtr->chargeFac() <= 0.) {
    if (verbose >= veryloud)
      printOut(__METHOD_NAME__, "Trial rejected (chargeFac <= 0)");
    return 0.;
  }
  bool isEmit = (iAntWin == iQQemitFF || iAntWin == iQGemitFF ||
                 iAntWin == iGQemitFF || iAntWin == iGGemitFF ||
                 iAntWin == iQQemitRF || iAntWin == iQGemitRF);

  // AlphaS, impose default choice. Can differ slighly from trial even
  // when running inside trial integral, due to flavor
  // thresholds. Here, alphaS(mu) is returned directly, with the
  // number of flavors active at mu, whereas the number of flavors in
  // the trial integral is controlled by the value of the trial scale.
  double alphaSNow = alphaSmax;
  if (alphaSorder >= 1) {
    double mu2 = getMu2(isEmit);
    AlphaStrong * alphaSptr = aSemitPtr;
    if(!isEmit) alphaSptr = aSsplitPtr;
    alphaSNow = min(alphaSmax, alphaSptr->alphaS(mu2));
  }

  // Compute physical antenna function (summed over final state
  // helicities). Note, physical antenna function can have swapped
  // labels (eg GQ -> GGQ).
  vector<double> mPost = winnerPtr->getmPostVec();
  vector<double> invariants = winnerPtr->getInvariants();
  unsigned int nPre = winnerPtr->iVec().size();
  vector<int> hPre = ( helicityShower && polarisedSys[iSysWin] ) ?
    winnerPtr->hVec() : vector<int>(nPre, 9);
  vector<int> hPost(nPre+1,9);
  double antPhys = antFunPtr->antFun(invariants, mPost, hPre, hPost);
  if (antPhys < 0.) {
    if (verbose > normal) infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Negative Antenna Function.", num2str(iAntWin));
    return 0.;
  }
  antPhys *= antFunPtr->chargeFac();

  // Check antenna function before multiplying by alphaS.
  if (doDiagnostics) diagnosticsPtr->checkAnt(iSysWin, antPhys);
  antPhys*=alphaSNow;
  return antPhys;

}

//--------------------------------------------------------------------------

// Calculate acceptance probability.

double VinciaFSR::pAcceptCalc(double antPhys) {
  double prob = winnerPtr->pAccept(antPhys,verbose);
  if (verbose >= louddebug)
    printOut(__METHOD_NAME__,"Shower pAccept = " + num2str(prob));
  return prob;
}

//--------------------------------------------------------------------------

// Generate the full kinematics.

bool VinciaFSR::genFullKinematics(int kineMap, Event event,
  vector<Vec4> &pPost) {

  // Generate branching kinematics, starting from antenna parents.
  vector<Vec4> pPre;
  vector<int> iPre          = winnerPtr->iVec();
  int nPre                  = iPre.size();
  int nPost                 = winnerPtr->iVec().size() + 1;
  vector<double> invariants = winnerPtr->getInvariants();
  vector<double> mPost      = winnerPtr->getmPostVec();
  int branchType            = winnerPtr->getBranchType();
  double phi                = 2 * M_PI * rndmPtr->flat();
  for (int i = 0; i < nPre; ++i) pPre.push_back(event[iPre[i]].p());

  // Special case for resonance decay.
  if (branchType == 5 || branchType == 6) {
    if (!vinComPtr->map2toNRFmassive(pPost, pPre, winnerPtr->posR(),
          winnerPtr->posF(), invariants,phi,mPost)) {
      if (verbose >= debug)
        printOut(__METHOD_NAME__, "Trial rejected (failed map2toNRF)");
      ++nFailedKine[iAntWin];
      return false;
    }
  } else {
    // 2->3 kinematics.
    if (nPre == 2 && nPost == 3) {
      if (!vinComPtr->map2to3FF(pPost, pPre, kineMap, invariants, phi,
          mPost)) {
        if (verbose >=  debug)
          printOut(__METHOD_NAME__, "Trial rejected (failed map2to3)");
        ++nFailedKine[iAntWin];
        return false;
      }
    // 2->4 kinematics
    } else if (nPre == 2 && nPost == 4) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": 2->4 kinematics map not implemented yet.");
      return false;
    // 3->4 kinematics
    } else if (nPre == 3 && nPost == 4) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": 3->4 kinematics map not implemented yet.");
      return false;
    }
  }
  return true;

}

//--------------------------------------------------------------------------

// Check if a trial is accepted.

bool VinciaFSR::acceptTrial(Event& event) {
  bool accept = false;
  bool doMEC  = doMECsSys[iSysWin];
  AntennaFunction* antFunPtr;

  // Check to see if we veto early before generating full kinematics,
  // i.e. just based on invariants.
  if (rejectEarly(antFunPtr,doMEC)) return accept;
  if (!getNewParticles(event,antFunPtr,pNew)) return accept;

  // Check sector veto.
  if (sectorShower) {
    vector<Particle> stateNew;
    stateNew = mecsPtr->makeParticleList(iSysWin,event,pNew,winnerPtr->iVec());
    double q2sector = resolutionPtr->q2sector2to3(&pNew[0],&pNew[2],&pNew[1]);
    vector<int> iSctDum;
    if (q2sector > resolutionPtr->findSector(iSctDum, stateNew)) {
      ++nSectorReject[iAntWin];
      return accept;
    }
  }

  // Check if phase space is closed for getting rid of heavy quarks.
  vector<Particle> stateOld;
  if (!isrPtr->checkHeavyQuarkPhaseSpace(stateOld,iSysWin)) {
    stateOld = mecsPtr->makeParticleList(iSysWin, event);
    if (verbose >= debug) printOut(__METHOD_NAME__, "Trial rejected (failed "
        "checkHeavyQuarkPhaseSpace)");
    // Mark this trial as "used", will need to generate a new one.
    ++nClosePSforHQ[iAntWin];
    return accept;
  }

  // TODO: matrix element corrections. If we want to compute a MEC,
  // we make a temporary new event record, and a temporary new
  // PartonSystem, with the tentative new state.
  double pMEC = 1.0;
  if (doMEC) pAccept[0] *= pMEC;

  // Count number of shower-type partons (for diagnostics and headroom
  // factors).
  int nQbef(0), nGbef(0), nBbef(0);
  for (int i = 0; i < partonSystemsPtr->sizeOut(iSysWin); ++i) {
    if (event[partonSystemsPtr->getOut(iSysWin,i)].id() == 21) ++nGbef;
    else if (event[partonSystemsPtr->getOut(iSysWin,i)].idAbs() <= 4) ++nQbef;
    else if (event[partonSystemsPtr->getOut(iSysWin,i)].idAbs() == 5) ++nBbef;
  }

  // Fill diagnostics histograms.
  if (verbose >= quiteloud || (verbose >= normal && doMEC)) {
    string state;
    if (nQbef >= 1) state += num2str(nQbef,1) + "q";
    if (nBbef >= 1) state += num2str(nBbef,1) + "b";
    if (nGbef >= 1) state += num2str(nBbef,1) + "g";
    if (pNew[1].colType() == 2) state += "Emit";
    else if (pNew[1].colType() == -1) state += "SplitA";
    else if (pNew[1].colType() == 1) state += "SplitB";
    string HPacc = "Log10(ME/AntTrial):" + state;
    if (vinciaHistos.find(HPacc) != vinciaHistos.end())
      vinciaHistos[HPacc].fill(log10(max(TINY,pAccept[0])));
    string HqTrial = "Ln(q2trial/sSystem):" + state;
    if (vinciaHistos.find(HqTrial) != vinciaHistos.end())
      vinciaHistos[HqTrial].fill(
        log(max(TINY, q2WinSav/pow2(mSystem[iSysWin]))));
  }

  // Print MC violations.
  if (doMEC && verbose >= verylouddebug) {
    stringstream ss;
    ss << " MEC pAccept = " << pAccept[0];
    printOut(__METHOD_NAME__, ss.str());
  }
  if (verbose >= loud ) {
    bool violation  = (pAccept[0] > 1.0 + TINY);
    bool negPaccept = (pAccept[0] < 0.0);
    if (violation) infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": pAccept > 1");
    if (negPaccept) infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": pAccept < 0");
    //Print more info for bad cases.
    if ((violation || negPaccept) && verbose > debug) winnerPtr->list();
  }

  // TODO: MC reweighting and uncertainty bands.

  // Enhance factors < 1 (radiation inhibition) treated by modifying pAccept.
  const double enhanceFac = winnerPtr->enhanceFac();
  if (rndmPtr->flat() > min(1.0,enhanceFac)*pAccept[0]) {
    if (verbose >= debug) printOut(__METHOD_NAME__ , "Trial rejected at veto "
        "step. wPhys/wTrial = " + num2str(pAccept[0])
        + " * enhanceFac = "+num2str(enhanceFac));

    // Reweighting to account for enhancement factor (reject).
    if (enhanceFac != 1.0)
      weightsPtr->scaleWeightEnhanceReject(pAccept[0],enhanceFac);

    // Count up number of vetoed branchings
    ++nFailedVeto[iAntWin];
  } else {
    if (verbose >= louddebug) printOut(__METHOD_NAME__, "Trial accepted.");

    // Reweighting to account for enhancement factor (accept).
    if (enhanceFac != 1.0) weightsPtr->scaleWeightEnhanceAccept(enhanceFac);

    // Fill diagnostics histos.
    if (verbose >= 3 || (verbose >= 2 && doMEC)) {
      string state;
      if (nQbef >= 1) state += num2str(nQbef,1) + "q";
      if (nBbef >= 1) state += num2str(nBbef,1) + "b";
      if (nGbef >= 1) state += num2str(nBbef,1) + "g";
      if (pNew[1].colType() == 2) state += "Emit";
      else if (pNew[1].colType() == -1) state += "SplitA";
      else if (pNew[1].colType() == 1) state += "SplitB";
      string HqPhys = "Ln(q2/sSystem):" + state;
      if (vinciaHistos.find(HqPhys) != vinciaHistos.end())
        vinciaHistos[HqPhys].fill(
          log(max(TINY, q2WinSav/pow2(mSystem[iSysWin]))));
      string HalphaS = "alphaS:" + state;
      string HalphaSratio = "alphaSratio:" + state;
      if (doMEC) {
        string HPacc = "Log10(ME/AntPhys):" + state;
        if (vinciaHistos.find(HPacc) != vinciaHistos.end()) {
          vinciaHistos[HPacc].fill(log10(max(TINY,pMEC)));
        }
      }
    }
    accept = true;
  }
  return accept;

}

//--------------------------------------------------------------------------

// Generate new particles for the antenna.

bool VinciaFSR::getNewParticles(Event& event, AntennaFunction* antFunPtr,
  vector<Particle>& newParts) {

  // Generate full kinematics.
  if (antFunPtr == nullptr) {
    if (verbose >= normal)
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": antFunPtr is NULL pointer.");
    return false;
  }
  newParts.clear();
  vector<Vec4> pPost;
  int maptype = antFunPtr->kineMap();
  if (!genFullKinematics(maptype, event, pPost)) {
    if (verbose >= debug)
      printOut(__METHOD_NAME__, "Failed to generate kinematics");
    return false;
  }

  // Generate new helicities.
  vector<int> hPost = genHelicities(antFunPtr);
  if (pPost.size() != hPost.size()) {
    if (verbose >= normal) {
      stringstream ss;
      ss << " pPost.size() = "
         << pPost.size() <<"  hPost.size() = " << hPost.size();
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Wrong size containers.", ss.str());
    }
    return false;
  } else if (!winnerPtr->getNewParticles(event, pPost, hPost, newParts,
      rndmPtr, colourPtr)) {
    if (verbose >= debug)
      printOut(__METHOD_NAME__, "Failed to generate new particles");
    return false;
  } else return true;

}

//--------------------------------------------------------------------------

// Generate new helicities for the antenna.

vector<int> VinciaFSR::genHelicities(AntennaFunction* antFunPtr) {

  vector<int> hPre = winnerPtr->hVec();
  vector<int> hPost = hPre;
  hPost.insert(hPost.begin() + 1, 9);
  if (hPost.size() >=3) {
    if (helicityShower && polarisedSys[iSysWin]) {
      vector<double> mPost = winnerPtr->getmPostVec();
      vector<double> invariants = winnerPtr->getInvariants();
      double helSum = antFunPtr->antFun(invariants, mPost, hPre, hPost);
      double randHel = rndmPtr->flat() * helSum;
      double aHel = 0.0;
      // Select helicity, n.b. positions hard-coded. hPost may be
      // larger than 3 in case of resonance decays but meaning of
      // first 3 positions is same (rest are unchanged).
      for(int iHel = 0; iHel < 8; ++iHel) {
        hPost[0] = ( (iHel%2)   )*2 -1;
        hPost[1] = ( (iHel/2)%2 )*2 -1;
        hPost[2] = ( (iHel/4)%2 )*2 -1;
        aHel = antFunPtr->antFun(invariants, mPost, hPre, hPost);
        randHel -= aHel;
        if (verbose >= verylouddebug) printOut(__METHOD_NAME__, "antPhys(" +
            num2str(int(hPre[0])) + " " + num2str(int(hPre[1])) + "  -> " +
            num2str(hPost[0]) + " " + num2str(hPost[1]) + " " +
            num2str(hPost[2]) + ") = " + num2str(aHel) + ", m(IK,ij,jk) = " +
            num2str(sqrt(invariants[0])) + ", " +
            num2str(sqrt(invariants[1])) + ", " +
            num2str(sqrt(invariants[2])) + "; sum = "+num2str(helSum));
        if (randHel < 0.) break;
      }
      if (doDiagnostics)
        diagnosticsPtr->checkAntHel(iSysWin, aHel, hPre, hPost);
    }
    if (verbose >= louddebug)
      printOut(__METHOD_NAME__, "selected"+num2str((int)(hPre[0]))
        + " " + num2str(int(hPre[1])) + "  -> " + num2str(hPost[0]) + " "
        + num2str(hPost[1]) + " " + num2str(hPost[2]) + "\n");
  }
  return hPost;

}

//--------------------------------------------------------------------------

// Update the event.

bool VinciaFSR::updateEvent(Event& event, resJunctionInfo& junctionInfoIn) {

  for (unsigned int i = 0; i < pNew.size(); ++i) event.append(pNew[i]);
  map<int, pair<int,int> >::iterator it;
  for (it = winnerPtr->mothers2daughters.begin();
       it != winnerPtr->mothers2daughters.end(); ++it) {
    int mother    = it->first;
    int daughter1 = (it->second).first;
    int daughter2 = (it->second).second;
    if (mother<event.size() && mother > 0) {
      event[mother].daughters(daughter1,daughter2);
      event[mother].statusNeg();
    } else return false;
  }

  // Add mothers to new daughters.
  for(it = winnerPtr->daughters2mothers.begin();
      it != winnerPtr->daughters2mothers.end(); ++it) {
    int daughter = it->first;
    int mother1  = (it->second).first;
    int mother2  = (it->second).second;
    if (daughter<event.size() && daughter > 0)
      event[daughter].mothers(mother1, mother2);
    else return false;
  }

  // Tell Pythia if we used a colour tag.
  if (winnerPtr->colTag() != 0) {
    int lastTag = event.nextColTag();
    int colMax  = winnerPtr->colTag();
    while (colMax > lastTag) lastTag = event.nextColTag();
  }
  iNewSav = winnerPtr->iNew();

  // Deal with any resonance junctions.
  if (hasResJunction[iSysWin]) {
    vector<int>* colours = &junctionInfoIn.colours;
    if (!event[junctionInfoIn.iEndQuark].isQuark()) {
      infoPtr->errorMsg("Error in "+__METHOD_NAME__
        +": Can't update junction. iEndQuark is not a quark!");
      hasResJunction[iSysWin]=false;
      return false;
    }

    // Check if resonance splitting.
    int branchType = winnerPtr->getBranchType();
    if (branchType == 2 || branchType == 6) {
      int splitter = winnerPtr->i0();
      if(branchType == 6) splitter = winnerPtr->iVec().at(winnerPtr->posF());
      // First update list of colours.
      int colLeft = event[splitter].col();
      // Find position col index.
      vector<int>::iterator pos = find(colours->begin(), colours->end(),
        colLeft);
      // Check if emission in string.
      if (pos != colours->end()) {
        // Remove part of string that has split off.
        colours->erase(pos + 1, colours->end());
        // Now update the junction info.
        int d1 = event[splitter].daughter1();
        int d2 = event[splitter].daughter2();
        if (event[d1].isQuark() && event[d1].col() > 0) {
          junctionInfoIn.iEndQuark  = d1;
          junctionInfoIn.iEndColTag = event[d1].col();
        } else if(event[d2].isQuark() && event[d2].col() > 0) {
          junctionInfoIn.iEndQuark  = d2;
          junctionInfoIn.iEndColTag = event[d2].col();
        }
        // Update junction.
        event.endColJunction(junctionInfoIn.iJunction, junctionInfoIn.iEndCol,
          junctionInfoIn.iEndColTag);
      }
    } else if (branchType == 1 || branchType == 5){
      //First update list of colours.
      int iNew = winnerPtr->iNew();

      // Find radiator (parton whose colours changed).
      int iRad = event[iNew].mother1();
      if(branchType == 1) {
        // Need to test both mothers.
        int m2 = event[iNew].mother2();
        if (m2 !=0) {
          // Get daughter that isn't iNew.
          int m2after=event[m2].daughter1();
          if (m2after==iNew) m2after = event[m2].daughter2();
          //Check, did this mother change colours or was it the
          // recoiler?
          int colBef    = event[m2].col();
          int acolBef   = event[m2].acol();
          int colAfter  = event[m2after].col();
          int acolAfter = event[m2after].acol();
          if(colBef != colAfter || acolBef != acolAfter) iRad = m2;
        }
      }

      //Find new colour to insert and old colour.
      int colNew = 0;
      int colLeft = event[iRad].col();
      if (event[iNew].col() == colLeft) colNew = event[iNew].acol();
      else colNew = event[iNew].col();
      if (colNew == 0) {
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Couldn't find colour for updating junction info.");
        return false;
      }

      // Find position of radiator col index.
      vector<int>::iterator pos = find(colours->begin(), colours->end(),
        colLeft);
      if (pos!=colours->end()) colours->insert(pos+1,colNew);

      // Now check if end quark has changed colour.
      int iEndQuark = junctionInfoIn.iEndQuark;
      if (!event[iEndQuark].isFinal()) {
        int d1 = event[iEndQuark].daughter1();
        int d2 = event[iEndQuark].daughter2();
        if (event[d1].isQuark() && event[d1].col() > 0) {
          junctionInfoIn.iEndQuark  = d1;
          junctionInfoIn.iEndColTag = event[d1].col();
        } else if (event[d2].isQuark() && event[d2].col() > 0) {
          junctionInfoIn.iEndQuark  = d2;
          junctionInfoIn.iEndColTag = event[d2].col();
        } else {
          infoPtr->errorMsg("Error in "+__METHOD_NAME__
            +": Couldn't update junction.");
          return false;
        }
        //Update junction.
        event.endColJunction(junctionInfoIn.iJunction,
          junctionInfoIn.iEndCol, junctionInfoIn.iEndColTag);
      }
    }
  }
  if (verbose >= louddebug) {
    printOut(__METHOD_NAME__, "Succesfully updated event after emission.");
    event.list();
  }
  return true;

}

//--------------------------------------------------------------------------

// Update the parton systems.

void VinciaFSR::updatePartonSystems() {

  // List parton systems.
  if (verbose >= verylouddebug) {
    printOut(__METHOD_NAME__, "Parton systems before update: ");
    partonSystemsPtr->list();
  }

  // Loop over mothers.
  vector<int> newpartons;
  for (map<int, pair<int,int> >::iterator it =
         winnerPtr->mothers2daughters.begin();
       it != winnerPtr->mothers2daughters.end(); ++it) {
    int mother    = it->first;
    int daughter1 = (it->second).first;
    int daughter2 = (it->second).second;
    // Two identical non-zero daughters -> recoilers, just update.
    if (daughter1 == daughter2 && daughter1 != 0 && daughter2 != 0) {
      partonSystemsPtr->replace(iSysWin, mother, daughter1);
      newpartons.push_back(daughter1);
    }
    // Two non-identical daughters -> brancher.
    else if (daughter1 != daughter2 && daughter1 != 0 && daughter2 != 0) {
      // Check if we have already added either daughter.
      bool found1 = false;
      bool found2 = false;
      vector<int>::iterator findit = find(newpartons.begin(), newpartons.end(),
        daughter1);
      if (findit != newpartons.end()) found1 = true;
      findit = find(newpartons.begin(), newpartons.end(), daughter2);
      if (findit != newpartons.end()) found2=true;
      // Both added already. Just continue.
      if (found1 && found2) continue;
      // 1 in record already - just update mother with 2.
      else if (found1 && !found2) {
        partonSystemsPtr->replace(iSysWin, mother, daughter2);
        newpartons.push_back(daughter2);
      // 2 in record already - just update mother with 1
      } else if (!found1 && found2) {
        partonSystemsPtr->replace(iSysWin, mother, daughter1);
        newpartons.push_back(daughter1);
      }
      // Neither in record, update mother with 1, add 2.
      else {
        partonSystemsPtr->replace(iSysWin, mother, daughter1);
        partonSystemsPtr->addOut(iSysWin, daughter2);
        newpartons.push_back(daughter1);
        newpartons.push_back(daughter2);
      }
    }
  }
  if (verbose >= verylouddebug) {
    printOut(__METHOD_NAME__, "Parton systems after update: ");
    partonSystemsPtr->list();
  }

}

//--------------------------------------------------------------------------

// Create a new emission brancher.

void VinciaFSR::saveEmitter(int iSysIn, Event& event, int i0, int i1) {
  if (event[i0].col() == event[i1].acol()) {
    emitters.push_back(BrancherEmitFF(iSysIn,event,i0,i1));
    lookupBrancherFF[make_pair(i0,true)]=(emitters.size()-1);
    lookupBrancherFF[make_pair(i1,false)]=(emitters.size()-1);
  }
}

//--------------------------------------------------------------------------

// Create a new resonance emission brancher.

void VinciaFSR::saveResEmitter(int iSysIn, Event& event, vector<int> allIn,
  unsigned int posResIn, unsigned int posFIn, bool colMode) {

  int iRes = allIn[posResIn];
  if (kMapResEmit == 2 && allIn.size() > 3) {
    // Save radiator.
    allIn.clear();
    int iRad = allIn[posFIn];
    int iRec = 0;
    int d1   = event[iRes].daughter1();
    int d2   = event[iRes].daughter2();

    // Find original colour connected.
    if (colMode) {
      // d2 was original recoiler.
      if (event[d1].col() > 0 && event[iRes].col() == event[d1].col())
        iRec = event[d2].iBotCopy();
      // d1 was original recoiler.
      else iRec = event[d1].iBotCopy();
    } else {
      // d2 was original recoiler.
      if(event[d1].acol() > 0 && event[iRes].acol() == event[d1].acol() )
        iRec = event[d2].iBotCopy();
      // d1 was original recoiler.
      else iRec = event[d1].iBotCopy();
    }
    allIn.push_back(iRes);
    allIn.push_back(iRad);
    allIn.push_back(iRec);
    posResIn=0;
    posFIn=1;
  }

  // Discriminate between colour and anticolour res antennae to avoid
  // degeneracy in lookupBrancherRF0 if res is colour octet.
  if (!colMode) iRes *= -1;

  // TODO: how to set zeta cut?
  BrancherEmitRF temp(iSysIn,event,allIn,posResIn,posFIn,q2CutoffEmit);
  resEmitters.push_back(temp);
  lookupBrancherRF[make_pair(iRes,true)]=(resEmitters.size() - 1);
  lookupBrancherRF[make_pair(allIn[posFIn],false)]=(resEmitters.size() - 1);

}

//--------------------------------------------------------------------------

// Create a new resonance splitter.

void VinciaFSR::saveResSplitter(int iSysIn, Event& event, vector<int> allIn,
  unsigned int posResIn, unsigned int posFIn,bool colMode) {

  int iRes = allIn[posResIn];
  if (kMapResSplit == 2 && allIn.size() > 3) {
    // Save radiator.
    allIn.clear();
    int iRad = allIn[posFIn];
    int iRec = 0;
    int d1   = event[iRes].daughter1();
    int d2   = event[iRes].daughter2();

    // Find original colour connected.
    if (colMode) {
      // d2 was original recoiler.
      if(event[d1].col() > 0 && event[iRes].col() == event[d1].col())
        iRec = event[d2].iBotCopy();
      // d1 was original recoiler.
      else iRec = event[d1].iBotCopy();
    } else {
      // d2 was original recoiler.
      if(event[d1].acol() > 0 && event[iRes].acol() == event[d1].acol())
        iRec = event[d2].iBotCopy();
      // d1 was original recoiler.
      else iRec = event[d1].iBotCopy();
    }
    allIn.push_back(iRes);
    allIn.push_back(iRad);
    allIn.push_back(iRec);
    posResIn=0;
    posFIn=1;
  }

  // Discriminate between colour and anticolour res antennae to avoid
  // degeneracy in lookupBrancherRF0 if res is colour octet.
  if (!colMode) iRes*=-1;
  BrancherSplitRF temp(iSysIn,event,allIn,posResIn,posFIn,q2CutoffSplit);
  resSplitters.push_back(temp);
  lookupSplitterRF[make_pair(iRes,true)]=(resSplitters.size() -1);
  lookupSplitterRF[make_pair(allIn[posFIn],false)]=(resSplitters.size() -1);

}

//--------------------------------------------------------------------------

// Create a new splitter brancher.

void VinciaFSR::saveSplitter(int iSysIn, Event& event, int i0, int i1,
  bool col2acol) {
  BrancherSplitFF temp(iSysIn, event, i0, i1, col2acol);
  splitters.push_back(temp);
  if (event[i0].isGluon()) {
    // Colour to anti-colour.
    if (col2acol) {
      lookupSplitter[make_pair(i0,true)]=(splitters.size()-1);
      lookupSplitter[make_pair(i1,false)]=(splitters.size()-1);
    // Anti-colour to colour.
    } else {
      lookupSplitter[make_pair(-i0,true)]=(splitters.size()-1);
      lookupSplitter[make_pair(-i1,false)]=(splitters.size()-1);
    }
  }
}

//--------------------------------------------------------------------------

// Update the branchers.

template <class Brancher> void VinciaFSR::updateBranchers(
  vector<Brancher>& brancherVec, map<pair<int, bool>,
  unsigned int>& lookupBrancher, Event& event, int iOld, int iNew) {

  pair<int,bool> key = make_pair(iOld, true);
  if (lookupBrancher.find(key) != lookupBrancher.end()) {
    unsigned int pos = lookupBrancher[key];
    int iRec         = brancherVec[pos].i1();
    int iSysNow      = brancherVec[pos].system();
    brancherVec[pos].reset(iSysNow, event, abs(iNew), iRec);
    lookupBrancher.erase(key);
    lookupBrancher[make_pair(iNew,true)] = pos;
  }
  key = make_pair(iOld, false);
  if (lookupBrancher.find(key) != lookupBrancher.end()) {
    unsigned int pos = lookupBrancher[key];
    int iEmit        = brancherVec[pos].i0();
    int iSysNow      = brancherVec[pos].system();
    brancherVec[pos].reset(iSysNow, event, iEmit, abs(iNew));
    lookupBrancher.erase(key);
    lookupBrancher[make_pair(iNew,false)]=pos;
  }

}

//--------------------------------------------------------------------------

// Update a single brancher.

template <class Brancher> void VinciaFSR::updateBrancher(
  vector<Brancher>& brancherVec, map<pair<int, bool>,
  unsigned int>& lookupBrancher, Event& event, int iOld1, int iOld2,
  int iNew1, int iNew2) {

  pair<int,bool> key1 = make_pair(iOld1,true);
  pair<int,bool> key2 = make_pair(iOld2,false);
  if (lookupBrancher.find(key1) != lookupBrancher.end()) {
    unsigned int pos = lookupBrancher[key1];
    if (lookupBrancher.find(key2) != lookupBrancher.end()) {
      unsigned int pos2=lookupBrancher[key2];
      if (pos == pos2) {
        lookupBrancher.erase(key1);
        lookupBrancher.erase(key2);
        int iSysNow = brancherVec[pos].system();
        brancherVec[pos].reset(iSysNow, event, abs(iNew1), abs(iNew2));
        lookupBrancher[make_pair(iNew1,true)]=pos;
        lookupBrancher[make_pair(iNew2,false)]=pos;
      }
    }
  }

}

//--------------------------------------------------------------------------

// Update emission branchers due to a recoiled parton.

void VinciaFSR::updateEmitters(Event& event, int iOld, int iNew){
  updateBranchers<BrancherEmitFF>(emitters,lookupBrancherFF,event,iOld,iNew);
}

//--------------------------------------------------------------------------

// Update emission brancher due to an emission.

void VinciaFSR::updateEmitter(Event& event,int iOld1, int iOld2,
  int iNew1, int iNew2) {
  updateBrancher<BrancherEmitFF>(emitters, lookupBrancherFF, event,
    iOld1, iOld2, iNew1, iNew2);
}

//--------------------------------------------------------------------------

// Update splitter branchers due to a recoiled parton.

void VinciaFSR::updateSplitters(Event& event, int iOld, int iNew) {
  updateBranchers<BrancherSplitFF>(splitters, lookupSplitter, event,
    iOld, iNew);
  updateBranchers<BrancherSplitFF>(splitters, lookupSplitter, event,
    -iOld, -iNew);
}

//--------------------------------------------------------------------------

// Update splitter brancher due to an emission.

void VinciaFSR::updateSplitter(Event& event,int iOld1, int iOld2, int iNew1,
  int iNew2, bool col2acol) {
  if (col2acol) updateBrancher<BrancherSplitFF>(splitters, lookupSplitter,
      event, iOld1, iOld2, iNew1, iNew2);
  else updateBrancher<BrancherSplitFF>(splitters, lookupSplitter, event,
      -iOld1, -iOld2, -iNew1, -iNew2);
}

//--------------------------------------------------------------------------

// Remove a splitter due to a gluon that has branched.

void VinciaFSR::removeSplitter(int iRemove) {

  for (int isign = 0; isign < 2; isign++) {
    int sign = 1 - 2*isign;
    pair<int,bool> key = make_pair(sign*iRemove, true);

    // Update map.
    if (lookupSplitter.find(key) != lookupSplitter.end()) {
      unsigned int pos = lookupSplitter[key];
      lookupSplitter.erase(key);
      int iRec = splitters[pos].i1();
      pair<int,bool> recoilkey = make_pair(sign*iRec, false);
      if (lookupSplitter.find(recoilkey) != lookupSplitter.end())
        lookupSplitter.erase(recoilkey);
      if (pos < splitters.size()) {
        splitters.erase(splitters.begin()+pos);

        // Update map with modified positions.
        for(; pos < splitters.size(); pos++) {
          //Get brancher at current pos.
          BrancherSplitFF splitter = splitters.at(pos);
          // Find indices.
          int i0(splitter.i0()), i1(splitter.i1());
          // Update lookup map to new pos.
          if (!splitter.isXG()){
            lookupSplitter[make_pair(i0,true)]=pos;
            lookupSplitter[make_pair(i1,false)]=pos;
          } else{
            lookupSplitter[make_pair(-i0,true)]=pos;
            lookupSplitter[make_pair(-i1,false)]=pos;
          }
        } // End loop over splitters.
      }
    }
  } // End loop over signs.

}

//--------------------------------------------------------------------------

// Update resonance emitter due to changed downstream decay products.

bool VinciaFSR::updateResBranchers(int iSysRes, Event& event, int iRes) {

  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "begin --------------");

  // Look up decay products using partonSystems, assumed to be updated
  // already. Colour information of resonance.
  int resCol      = event[iRes].col();
  int resACol     = event[iRes].acol();
  int colPartner  = -1;
  int acolPartner = -1;
  vector<int> daughters;

  // Loop over members of current decay system and get colour information.
  int sizeOut = partonSystemsPtr->sizeOut(iSysRes);
  for (int iDecPart = 0; iDecPart < sizeOut; iDecPart++) {
    int iOut = partonSystemsPtr->getOut(iSysRes,iDecPart);

    // Check if colouredm partner of the resonance.
    if (event[iOut].col() != 0 && event[iOut].col() == resCol)
      colPartner = iOut;
    if (event[iOut].acol() != 0 && event[iOut].acol() == resACol)
      acolPartner = iOut;
    if (iOut != colPartner && iOut != acolPartner)
      daughters.push_back(iOut);
  }
  if (verbose >= verylouddebug) {
    stringstream ss;
    ss << "col partner = " << colPartner << " acol partner = " << acolPartner;
    printOut(__METHOD_NAME__,ss.str());
  }

  if (colPartner > 0) {
    // Get a copy of daughters.
    vector<int> resSysAll = daughters;
    if (acolPartner != colPartner && acolPartner > 0)
      resSysAll.push_back(acolPartner);
    // Insert col partner and res at front (just convention).
    resSysAll.insert(resSysAll.begin(),colPartner);
    resSysAll.insert(resSysAll.begin(),iRes);
    unsigned int posRes(0), posPartner(1);
    updateResBranchers(iSysRes,event,resSysAll,posRes,posPartner,true);
  }
  if (acolPartner > 0) {
    // Get a copy of daughters.
    vector<int> resSysAll = daughters;
    if (acolPartner != colPartner && colPartner > 0)
      resSysAll.push_back(colPartner);
    // Insert col partner and res at front (just convention).
    resSysAll.insert(resSysAll.begin(),acolPartner);
    resSysAll.insert(resSysAll.begin(),iRes);
    unsigned int posRes(0), posPartner(1);
    updateResBranchers(iSysRes,event,resSysAll,posRes,posPartner,false);
  }
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "end --------------");

  return true;

}

//--------------------------------------------------------------------------

// Update resonance emitter due to changed downstream decay products.

void VinciaFSR::updateResBranchers(int iSysRes, Event& event,
  vector<int> resSysAll, unsigned int posRes, unsigned int posPartner,
  bool isCol) {

  if (posRes >= resSysAll.size() || posPartner >= resSysAll.size()) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__+": Invalid positions.");
    infoPtr->setAbortPartonLevel(true);
    return;
  }
  int iRes      = resSysAll[posRes];
  int iPartner  = resSysAll[posPartner];
  int posREmit  = posRes;
  int posFEmit  = posPartner;
  int posRSplit = posRes;
  int posFSplit = posPartner;
  vector<int> resSysAllEmit;
  vector<int> resSysAllSplit;

  // If "bad" recoil map need to update recoiler system resSysAll.
  if (kMapResEmit == 2 && resSysAll.size() > 3) {
    // Fetch daughters of res.
    int iRec = 0;
    int d1 = event[iRes].daughter1();
    int d2 = event[iRes].daughter2();

    // Find original colour connected.
    if (isCol) {
      // d2 was original recoiler.
      if (event[d1].col() > 0 && event[iRes].col() == event[d1].col())
        iRec = event[d2].iBotCopy();
      // d1 was original recoiler.
      else iRec = event[d1].iBotCopy();
    } else {
      // d2 was original recoiler.
      if (event[d1].acol() > 0 && event[iRes].acol() == event[d1].acol())
        iRec = event[d2].iBotCopy();
      // d1 was original recoiler.
      else iRec = event[d1].iBotCopy();
    }
    resSysAllEmit.push_back(iRes);
    resSysAllEmit.push_back(iPartner);
    resSysAllEmit.push_back(iRec);
    posREmit=0;
    posFEmit=1;
  } else resSysAllEmit = resSysAll;
  if (kMapResSplit == 2) {
    resSysAllSplit = resSysAllEmit;
    posRSplit=0;
    posFSplit=1;
  } else resSysAllSplit = resSysAll;
  if (!isCol) iRes*=-1;

  // First update emission brancher -> always need to update because
  // downstream recoilers will have changed.
  pair<int,bool> branchkey = make_pair(iRes, true);
  if (lookupBrancherRF.find(branchkey) != lookupBrancherRF.end()) {
    unsigned int pos =lookupBrancherRF[branchkey];
    int iRec = (resEmitters[pos].iVec())[resEmitters[pos].posF()];
    pair<int,bool> recoilkey=make_pair(iRec,false);
    // Delete map to recoiler.
    if (lookupBrancherRF.find(recoilkey) != lookupBrancherRF.end())
      lookupBrancherRF.erase(recoilkey);
    // Reset brancher.
    resEmitters[pos].resetResBrancher(iSysRes, event, resSysAllEmit, posREmit,
      posFEmit,q2CutoffEmit);
    // Add new map.
    recoilkey = make_pair(iPartner,false);
    lookupBrancherRF[recoilkey] = pos;
  }

  // Splitters - treatement depends on latest emission.
  if (lookupSplitterRF.find(branchkey) != lookupSplitterRF.end()){
    unsigned int pos = lookupSplitterRF[branchkey];
    int iSplit = (resSplitters[pos].iVec())[resSplitters[pos].posF()];
    pair<int,bool> splitkey=make_pair(iSplit,false);

    // Delete map to recoiler.
    if (lookupSplitterRF.find(splitkey) != lookupSplitterRF.end())
      lookupSplitterRF.erase(splitkey);

    //Do we need to remove this splitter, is the splitter still a gluon?
    if (!event[iPartner].isGluon()) {
      lookupSplitterRF.erase(branchkey);
      resSplitters.erase(resSplitters.begin()+pos);
      // Update any other splitters' positions in lookup map.
      for (unsigned int ipos = pos; ipos < resSplitters.size(); ipos++) {
        BrancherSplitRF splitter = resSplitters.at(ipos);
        int itmpSplit = (resSplitters[ipos].iVec())[resSplitters[ipos].posF()];
        // Update lookup map to new pos.
        lookupSplitterRF[make_pair(iRes,true)] = ipos;
        lookupSplitterRF[make_pair(itmpSplit,false)] = ipos;
      }
    // Otherwise just update.
    } else {
      resSplitters[pos].resetResBrancher(iSysRes, event, resSysAllSplit,
        posRSplit, posFSplit, q2CutoffSplit);
      // Add new map.
      splitkey = make_pair(iPartner,false);
      lookupSplitterRF[splitkey]=pos;
    }
  // Else if last branch type was res branch add new res splitter.
  } else if (winnerPtr!= nullptr) {
    if (winnerPtr->getBranchType()==5 && event[iPartner].isGluon())
      saveResSplitter(iSysRes,event,resSysAllSplit,posRSplit,posFSplit,isCol);
  }

}

//--------------------------------------------------------------------------

// Update the antennae.

bool VinciaFSR::updateAntennae(Event& event) {

  if (verbose >= debug) {
    printOut(__METHOD_NAME__, "begin --------------");
    if (verbose >= superdebug) printLookup();
  }
  if (winnerPtr == nullptr) {
    if (verbose >= normal) infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": winnerPtr is NULL");
    return false;
  }

  // Update QED system(s), then QCD.
  if (doQED) qedShowerPtr->update(event, iSysWin);

  // Was this g->ffbar?
  int branchType = winnerPtr->getBranchType();
  if (branchType == 2) {
    // Remove old splitters where g = i0 and update splitters where g is rec.
    int splitOld = winnerPtr->i0();
    int recOld = winnerPtr->i1();
    removeSplitter(splitOld);

    // Get daughters.
    int iColSplit = event[splitOld].daughter1();
    int iaColSplit = event[splitOld].daughter2();
    if(event[iColSplit].col() == 0 && event[iaColSplit].acol() == 0 &&
       event[iColSplit].acol() != 0 && event[iaColSplit].col() != 0) {
      iColSplit = event[splitOld].daughter2();
      iaColSplit = event[splitOld].daughter1();
    }

    //Find colour connected partner.
    int iColPartner = 0;
    int iaColPartner = 0;
    pair<int,bool> testkey = make_pair(splitOld,true);
    if (lookupBrancherFF.find(testkey) != lookupBrancherFF.end()) {
      int iTest = emitters[lookupBrancherFF[testkey]].i1();
      if(event[iTest].acol() == event[splitOld].col()) iaColPartner = iTest;
    }
    testkey = make_pair(splitOld,false);
    if (lookupBrancherFF.find(testkey) != lookupBrancherFF.end()) {
      int iTest = emitters[lookupBrancherFF[testkey]].i0();
      if (event[iTest].col() == event[splitOld].acol()) iColPartner = iTest;
    }

    //Update splitters where g is (anti-)colour-connected recoiler/emitter.
    updateSplitter(event,iaColPartner,splitOld,iaColPartner,iColSplit,false);
    updateSplitter(event,iColPartner,splitOld,iColPartner,iaColSplit,true);
    updateEmitter(event,iColPartner,splitOld,iColPartner,iaColSplit);
    updateEmitter(event,splitOld,iaColPartner,iColSplit,iaColPartner);

    // Update recoiler.
    int recNew = event[recOld].daughter1();
    updateSplitters(event,recOld,recNew);
    updateEmitters(event,recOld,recNew);
  }

  // Emission.
  else if (branchType == 1) {
    // Update old splitters.
    int iOld1 = winnerPtr->i0();
    int iOld2 = winnerPtr->i1();
    int iNew1 = event[iOld1].daughter1();
    int iNew2 = event[iOld1].daughter2();
    int iNew3 = event[iOld2].daughter1();

    // Switch 1<->2 so that 2 is repeated daughter.
    if (iNew3 == iNew1) {
      iNew1=iNew2;
      iNew2=iNew3;
      iNew3=event[iOld2].daughter2();
    } else if (iNew3 == iNew2) iNew3=event[iOld2].daughter2();

    // Update emitters, determine antenna to preserve.
    // ab->12.
    if (event[iOld1].col() == event[iNew1].col()) {
      updateEmitter(event,iOld1,iOld2,iNew1,iNew2);
      if(event[iNew2].col() == event[iNew3].acol())
        saveEmitter(iSysWin,event,iNew2,iNew3);
    // ab->23.
    } else {
      updateEmitter(event,iOld1,iOld2,iNew2,iNew3);
      if(event[iNew1].col()==event[iNew2].acol())
        saveEmitter(iSysWin,event,iNew1,iNew2);
    }
    if (event[iNew1].isGluon())
      updateSplitter(event,iOld1,iOld2,iNew1,iNew2,true);
    if (event[iNew3].isGluon())
      updateSplitter(event,iOld2,iOld1,iNew3,iNew2,false);

    // New splitters.
    if (event[iNew2].isGluon()) {
      saveSplitter(iSysWin,event,iNew2,iNew3,true);
      saveSplitter(iSysWin,event,iNew2,iNew1,false);
    }

    // Update other connected-connected antenna, excluding antenna
    // which branched.
    updateEmitters(event,iOld1,iNew1);
    updateEmitters(event,iOld2,iNew3);
    updateSplitters(event,iOld1,iNew1);
    updateSplitters(event,iOld2,iNew3);

  // Resonance emission.
  } else if (branchType == 5 || branchType == 6) {
    // Update emitters and splitters.
    for (map<int, pair<int, int> >::iterator it =
         winnerPtr->mothers2daughters.begin();
         it!= winnerPtr->mothers2daughters.end(); ++it){
      int mother    = it->first;
      int daughter1 = (it->second).first;
      int daughter2 = (it->second).second;
      // Recoiler -> just update.
      if (daughter1 == daughter2) {
        updateEmitters(event,mother,daughter1);
        updateSplitters(event,mother,daughter1);
      // Resonance emitter.
      } else {
        // Convention of res emission: daughter1 is new emission but
        // check anyway.
        if (branchType == 5 && event[daughter1].isGluon()) {
          if (event[daughter1].col()==event[daughter2].acol())
            saveEmitter(iSysWin,event,daughter1,daughter2);
          else if(event[daughter1].acol()==event[daughter2].col())
            saveEmitter(iSysWin,event,daughter2,daughter1);
          // TODO: check colour condition here.
          bool col2acol = false;
          if (event[daughter1].col() == event[daughter2].acol())
            col2acol = true;
          saveSplitter(iSysWin,event,daughter1,daughter2,col2acol);
          updateEmitters(event,mother,daughter2);
          updateSplitters(event,mother,daughter2);
        // Resonant splitter.
        } else if (branchType == 6 && event[mother].isGluon()
          && !event[daughter1].isGluon() && !event[daughter2].isGluon()) {
          removeSplitter(mother);
          int iColSplit  = daughter1;
          int iaColSplit = daughter2;
          if(event[mother].col() != event[daughter1].col()) {
            iColSplit  = daughter2;
            iaColSplit = daughter1;
          }

          // Find colour connected partner.
          int iColPartner(0), iaColPartner(0);
          pair<int,bool> testkey = make_pair(mother,true);
          if (lookupBrancherFF.find(testkey) != lookupBrancherFF.end()) {
            int iTest = emitters[lookupBrancherFF[testkey]].i1();
            if (event[iTest].acol() == event[mother].col())
              iaColPartner=iTest;
          }
          testkey = make_pair(mother,false);
          if (lookupBrancherFF.find(testkey) != lookupBrancherFF.end()) {
            int iTest = emitters[lookupBrancherFF[testkey]].i0();
            if (event[iTest].col() == event[mother].acol())
              iColPartner=iTest;
          }

          // Update splitters where mother was a
          // (anti)colour-connected recoiler/emitter.
          updateSplitter(event, iaColPartner, mother, iaColPartner, iColSplit,
            false);
          updateSplitter(event, iColPartner, mother, iColPartner, iaColSplit,
            true);
          updateEmitter(event, iColPartner, mother, iColPartner, iaColSplit);
          updateEmitter(event, mother, iaColPartner, iColSplit, iaColPartner);

        }
      } // End found branching in mothers2daughters.
    } // End for loop over mothers2daughters.
  } // End resonance brancher case.

  // If system containing resonance has branched, must always update
  // (because must update downstream recoilers regardless of if last
  // branch was res emission).
  if (isResonanceSys[iSysWin]) {
    if (!updateResBranchers(iSysWin, event,
        partonSystemsPtr->getInRes(iSysWin))){
      if (verbose >= normal)
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed updateResEmitters.");
      return false;
    }
  }
  if (verbose >= debug) {
    if (verbose >= louddebug) list();
    if (verbose >= superdebug) {
      printLookup();
      if (!check(event)) {
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed update antennae.");
        return false;
      }
    }
    printOut(__METHOD_NAME__, "end --------------");
  }
  return true;

}

//--------------------------------------------------------------------------

// Update systems of QCD antennae after a QED branching.

bool VinciaFSR::updateAfterQED(Event& event, int sizeOld) {

  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "begin --------------");
  bool updated = false;
  iSysWin = qedShowerPtr->sysWin();

  // Create colour-anticolour map for post-branching partons.
  map<int,int> iOfCol, iOfAcol;
  // Also check for coloured partons that were created in the
  // splitting, i.e. with status 51, used to create new emission
  // antennae below.
  vector<int>  status51coloured;
  vector<pair<int,int> > colouredrecoilers;
  for (int i = sizeOld; i < event.size(); ++i) {
    int col  = event[i].col();
    int acol = event[i].acol();
    if (col != 0) iOfCol[col] = i;
    if (acol != 0) iOfAcol[acol] = i;
    // Check which were "created" (as opposed to recoiling) - to see
    // if we need to create splitter.
    if (event[i].colType() != 0 && event[i].status() == 51)
      status51coloured.push_back(i);
    else if (event[i].colType() != 0 &&
             (event[i].status() == 52 || event[i].status() == 43 ||
              event[i].status() == 44)) {
      int moth = event[i].mother1();
      if (moth > 0) colouredrecoilers.push_back(make_pair(moth, i ));
    }
  }

  if (status51coloured.size() == 2) {
    int i1    = status51coloured[0];
    int i2    = status51coloured[1];
    int iCol  = event[i1].colType() > 0 ? i1 : i2;
    int iAcol = event[i1].colType() > 0 ? i2 : i1;
    // If this was a splitting to coloured partons, create new
    // emission antenna. Create a QCD emission antenna between the two
    // status-51 partons.
    if (qedShowerPtr->isTrialSplit) saveEmitter(iSysWin,event,iCol,iAcol);
    // Need to update existing QCD antennae.
    else {
      int moth1 = event[i1].mother1();
      colouredrecoilers.push_back(make_pair(moth1, i1));
      int moth2 = event[i2].mother1();
      colouredrecoilers.push_back(make_pair(moth2, i2));
      if(event[moth1].col() == event[moth2].acol())
        updateEmitter(event,moth1,moth2,i1,i2);
    }
  } else if (status51coloured.size() == 1) {
    int i1    = status51coloured[0];
    int moth1 = event[i1].mother1();
    colouredrecoilers.push_back(make_pair(moth1, i1));
  } else if (status51coloured.size() > 2){
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Too many status 51 particles");
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // Reset any emission antennae involving quarks that have recoiled.
  for(vector<pair<int,int> >::iterator it = colouredrecoilers.begin();
      it!= colouredrecoilers.end(); ++it) {
    int recOld = it->first;
    int recNew = it->second;
    updateEmitters(event,recOld,recNew);
    updateSplitters(event,recOld,recNew);
  }

  // Update resonance antennae.
  if (isResonanceSys[iSysWin]){
    if (!updateResBranchers(iSysWin, event,
        partonSystemsPtr->getInRes(iSysWin))) {
      if (verbose >= normal)
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Failed updateResEmitters.");
      return updated;
    }
  }

  // Check the event.
  if (!check(event)) {
    infoPtr->errorMsg("Error in "+__METHOD_NAME__
      +": Failed update antennae.");
    list();
    if (verbose >= superdebug) printLookup();
    infoPtr->setAbortPartonLevel(true);
    return false;
  }

  // Now check if end quark has changed.
  if (hasResJunction[iSysWin]) {
    int iEndQuark = junctionInfo[iSysWin].iEndQuark;
    if (!event[iEndQuark].isFinal()) {
      int d1 = event[iEndQuark].daughter1();
      int d2 = event[iEndQuark].daughter2();
      if (event[d1].isQuark() && event[d1].col() > 0)
        junctionInfo[iSysWin].iEndQuark = d1;
      else if(event[d2].isQuark() && event[d2].col() > 0)
        junctionInfo[iSysWin].iEndQuark = d2;
      else {
        infoPtr->errorMsg("Error in "+__METHOD_NAME__
          +": Couldn't update junction information");
        return false;
      }
    }
  }
  if (verbose >= verylouddebug)
    printOut(__METHOD_NAME__, "end --------------");

  return true;

}

//--------------------------------------------------------------------------

// Print a brancher lookup.

void VinciaFSR::printLookup(map< pair<int, bool>, unsigned int>
  lookupBrancher, string name) {
  for (map< pair<int, bool>, unsigned int >::iterator ilook =
         lookupBrancher.begin(); ilook != lookupBrancher.end(); ++ilook)
    cout << "  lookup" << name << "[" << (ilook->first).first
         << "," << (ilook->first).second << "] = " << ilook->second << endl;
}

//--------------------------------------------------------------------------

// Print the brancher lookup maps.

void VinciaFSR::printLookup() {
  cout << endl << "  --------" << "  Brancher lookup maps"
       << "  -------------------------------------------------------------"
       << endl;
  printLookup(lookupBrancherRF,"BrancherRF");
  printLookup(lookupSplitterRF,"SplitterRF");
  printLookup(lookupBrancherFF,"BrancherFF");
  printLookup(lookupSplitter,"SplitterFF");
  cout << "  --------" << "       End lookup     "
       <<"  -------------------------------------------------------------"
       << endl << endl;
}

//--------------------------------------------------------------------------

// Calculate the headroom factor.

vector<double> VinciaFSR::getHeadroom(int iSys, bool isEmit, double) {

  // TODO: ensure a decent number of failed trials if doing uncertainties.
  bool doUncert = false;

  // Check if we have we encountered this headroom criterion before.
  pair<int, pair<bool,bool> > headroomKey =
    make_pair(iSys,make_pair(isEmit,doUncert));
  if (headroomSav.find(headroomKey) != headroomSav.end())
    return headroomSav[headroomKey];
  // Otherwise calculate, and save for next time we need it.
  else {
    // Emissions.
    vector<double> headroomVec;
    if (isEmit) {
      double headroomFac = 1.0;
      // Increased headroom factor if doing MECs and/or uncertainties.
      if (doMECsSys[iSys] && mecsPtr->doMEC(iSys,nBranch[iSys]+1)) {
        headroomFac = 1.5;
        // More headroom for 2->2 than for resonance decays.
        if (!isResonanceSys[iSys]) headroomFac *= 2.;
        // More headroom for helicity dependence.
        if (helicityShower && polarisedSys[iSys]) headroomFac *= 1.5;
      }
      if (doUncert) headroomFac *= 1.33;
      headroomVec.push_back(headroomFac);

    // Other.
    } else {
      for (int iFlav = 1; iFlav <= nGluonToQuark; ++iFlav) {
        double headroomFac = 1.0;
        // For sector showers, trial probability should be twice as large.
        if (sectorShower) headroomFac *= 2;
        // Heavy flavours get 50% larger trial (since mass correction > 0).
        if (iFlav > nFlavZeroMass) headroomFac *= 1.5;
        // MECs also get increased headroom.
        if (doMECsSys[iSys] && mecsPtr->doMEC(iSys,nBranch[iSys]+1)) {
          headroomFac *= 2.;
          // More headroom for 2->2 than for resonance decays.
          if (!isResonanceSys[iSys]) headroomFac *= 2.;
          // More headroom for helicity dependence.
          if (helicityShower && polarisedSys[iSys]) headroomFac *= 2.;
        }
        headroomVec.push_back(headroomFac);
      }
    }
    headroomSav[headroomKey] = headroomVec;
    return headroomVec;
  }

}

//--------------------------------------------------------------------------

// Calculate the enhancement factor.

vector<double> VinciaFSR::getEnhance(int iSys, bool isEmit, double q2In) {

  bool doEnhance = false;
  if (q2In > pow2(enhanceCutoff)) {
    if (isHardSys[iSys] && enhanceInHard) doEnhance = true;
    else if (isResonanceSys[iSys] && enhanceInResDec) doEnhance = true;
    else if (!isHardSys[iSys] && !isResonanceSys[iSys] &&
             partonSystemsPtr->hasInAB(iSys) && enhanceInMPI) doEnhance = true;
  }

  // Check if we have encountered this enhancement criterion before.
  pair<int,pair<bool,bool> > enhanceKey =
    make_pair(iSys,make_pair(isEmit,doEnhance));
  vector<double> enhanceVec;
  if (enhanceSav.find(enhanceKey) != enhanceSav.end())
    enhanceVec = enhanceSav[enhanceKey];
  else {
    double enhanceFac = 1.0;
    // Emissions.
    if (isEmit) {
      if (doEnhance) enhanceFac *= enhanceAll;
      enhanceVec.push_back(enhanceFac);
    // Other.
    } else {
      for (int iFlav = 1; iFlav <= nGluonToQuark; ++iFlav) {
        if (doEnhance) {
          enhanceFac = enhanceAll;
          // Optional extra enhancement for g->cc, g->bb.
          if (iFlav == 4) enhanceFac *= enhanceCharm;
          else if (iFlav == 5) enhanceFac *= enhanceBottom;
        }
        enhanceVec.push_back(enhanceFac);
      }
    }
    enhanceSav[enhanceKey] = enhanceVec;
  }
  return enhanceVec;
}

//==========================================================================

} // end namespace Pythia8
