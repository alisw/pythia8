// MiniStringFragmentation.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the .
// MiniStringFragmentation class

#include "Pythia8/MiniStringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The MiniStringFragmentation class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Since diffractive by definition is > 1 particle, try hard.
const int MiniStringFragmentation::NTRYDIFFRACTIVE = 200;

// After one-body fragmentation failed, try two-body once more.
const int MiniStringFragmentation::NTRYLASTRESORT  = 100;

// Loop try to combine available endquarks to valid hadron.
const int MiniStringFragmentation::NTRYFLAV        = 10;

//--------------------------------------------------------------------------

// Initialize and save pointers.

void MiniStringFragmentation::init(StringFlav* flavSelPtrIn,
  StringPT* pTSelPtrIn, StringZ* zSelPtrIn) {

  // Save pointers.
  flavSelPtr      = flavSelPtrIn;
  pTSelPtr        = pTSelPtrIn;
  zSelPtr         = zSelPtrIn;

  // Calculation and definition of hadron space-time production vertices.
  hadronVertex    = mode("HadronVertex:mode");
  setVertices     = flag("Fragmentation:setVertices")
                 || flag("HadronLevel:Rescatter");
  kappaVtx        = parm("HadronVertex:kappa");
  smearOn         = flag("HadronVertex:smearOn");
  xySmear         = parm("HadronVertex:xySmear");
  constantTau     = flag("HadronVertex:constantTau");

  // Charm and bottom quark masses used for space-time offset.
  mc              = particleDataPtr->m0(4);
  mb              = particleDataPtr->m0(5);

  // Initialize the MiniStringFragmentation class proper.
  nTryMass        = mode("MiniStringFragmentation:nTry");

  // Initialize the b parameter of the z spectrum, used when joining jets.
  bLund           = zSelPtr->bAreaLund();

}

//--------------------------------------------------------------------------

// Do the fragmentation: driver routine.

bool MiniStringFragmentation::fragment(int iSub, ColConfig& colConfig,
  Event& event, bool isDiff, bool systemRecoil) {

  // Check for junction topologies
  iParton = colConfig[iSub].iParton;
  isJunctionSystem = colConfig[iSub].hasJunction;
  SaveJunctionState saveJunctionState(*this, event);

  // First check if iParton is a junction system.
  if (iParton.front() < 0) {

    // If so, then save some information and reduce the junction to a
    // simple junction by absorbing any gluon momenta into the (di)quarks.
    saveJunctionState.saveMomenta();
    if (iParton.size() > 6 && !reduce2SimpleJunction(event)) {
      loggerPtr->ERROR_MSG("failed to reduce the size of junction system "
        "containing gluons");
      return false;
    }

    // Then save some more information.
    flav1  = FlavContainer( event[iParton[1]].id() );
    flav2  = FlavContainer( event[iParton[3]].id() );
    flavj3 = FlavContainer( event[iParton[5]].id() );
    pSum   = colConfig[iSub].pSum;
    mSum   = colConfig[iSub].mass;
    m2Sum  = mSum*mSum;

    // If all the ends are quarks, simply join two of them into a diquark
    // and continue with the normal ministring machinery.
    if (!flav1.isDiquark() && !flav2.isDiquark() && !flavj3.isDiquark())
      reduce2SimpleString(event);

    // If a diquark is involved we need to create two hadrons.
    else {
      if (minijunction2two( nTryMass, event)) return true;
      loggerPtr->ERROR_MSG("minijunction2Hadrons failed");
      return false;
    }
  }

  // Read in info on system to be treated.
  flav1    = FlavContainer( event[ iParton.front() ].id() );
  flav2    = FlavContainer( event[ iParton.back() ].id() );
  pSum     = colConfig[iSub].pSum;
  mSum     = colConfig[iSub].mass;
  m2Sum    = mSum*mSum;
  isClosed = colConfig[iSub].isClosed;

  // Do not want diffractive systems to easily collapse to one particle.
  int nTryFirst = (isDiff) ? NTRYDIFFRACTIVE : nTryMass;

  // First try to produce two hadrons from the system.
  if (ministring2two( nTryFirst, event, false)) return true;

  // If this fails, then form one hadron and shuffle momentum.
  if (ministring2one( iSub, colConfig, event, false)) return true;

  // If also this fails, try to produce two hadrons with lower mass.
  if (ministring2two( NTRYLASTRESORT, event, true)) return true;

  // If also this fails, try to form a hadron with lower mass.
  if (ministring2one( iSub, colConfig, event, true)) return true;

  // For low-energy systems may also search for a single hadron recoiler.
  if (!systemRecoil) {
    if (ministring2one( iSub, colConfig, event, false, false)) return true;
    if (ministring2one( iSub, colConfig, event, true,  false)) return true;
  }

  // Else complete failure.
  loggerPtr->ERROR_MSG("no 1- or 2-body state found above mass threshold");
  return false;

}

//--------------------------------------------------------------------------

// Reduce the junction size by absorbing the gluons into quarks/diquarks.

bool MiniStringFragmentation::reduce2SimpleJunction(Event& event) {

  // Find summed momentum and first and last parton on each leg.
  vector<Vec4> pleg;
  vector<int> legj, legq;
  for (int i = 0; i < int(iParton.size()); ++i) {
    if (iParton[i] < 0) {
      pleg.push_back(Vec4());
      legq.push_back(0);
      legj.push_back(iParton[i]);
    } else {
      legq.back() = iParton[i];
      pleg.back() += event[iParton[i]].p();
    }
  }

  // Failure if system contains more than three endpoints.
  if ( pleg.size() != 3 ) {
    loggerPtr->ERROR_MSG("cannot process multi-junction system");
    return false;
  }

  // Put summed momentum in endpoint quark and simplify parton system.
  for (int i = 0; i < 3; ++i) event[legq[i]].p(pleg[i]);
  iParton = { legj[0], legq[0], legj[1], legq[1], legj[2], legq[2] };
  return true;

}

//--------------------------------------------------------------------------

// Reduce the junction to a string by merging the lightest q + q combination
// into a diquark.

void MiniStringFragmentation::reduce2SimpleString(Event& event){

  // Parton pair (squared) invariant masses.
  Vec4 p1 = event[iParton[1]].p();
  Vec4 p2 = event[iParton[3]].p();
  Vec4 p3 = event[iParton[5]].p();
  double m12 = (p1 + p2).m2Calc();
  double m13 = (p1 + p3).m2Calc();
  double m23 = (p2 + p3).m2Calc();

  // Identify pair with smallest mass. Switch index from iParton to event.
  int ip1 = 1, ip2 = 3, ip3 = 5;
  if (m13 > m12) {ip1 = 1; ip2 = 5; ip3 = 3;}
  if (m23 > max(m12, m13)) {ip1 = 3; ip2 = 5; ip3 = 1;}
  ip1 = iParton[ip1];
  ip2 = iParton[ip2];
  ip3 = iParton[ip3];

  // Construct the diquark and add to event record with correct colour index.
  Vec4 pqq = event[ip1].p() + event[ip2].p();
  int idqq = flavSelPtr->makeDiquark(event[ip1].id(), event[ip2].id());
  int ipqq = event.append( idqq, 78, 0, 0, 0, 0, 0, 0, pqq, pqq.mCalc());
  if (idqq > 0) event[ipqq].acol(event[ip3].col());
  else          event[ipqq].col(event[ip3].acol());
  iParton = {ip3, ipqq};

  // Production vertex of diquark as average.
  if (setVertices) {
    Vec4 vqq = 0.5 * (event[ip1].vProd() + event[ip2].vProd());
    event[ipqq].vProd( vqq);
  }
  return;

}

//--------------------------------------------------------------------------

// Attempt to produce two hadrons from the minijunction.

bool MiniStringFragmentation::minijunction2two( int nTry, Event& event) {

  // Properties of the produced hadrons.
  pair<int,int> idHad;
  double mHad1   = 0.;
  double mHad2   = 0.;
  double mHadSum = 0.;

  // Order junction ends in increasing idAbs (and thereby likely mass).
  int ip[3]    = {iParton[1], iParton[3], iParton[5]};
  int idAbs[3] = {abs(flav1.id), abs(flav2.id), abs(flavj3.id)};
  if (idAbs[0] > idAbs[1]) {swap( ip[0], ip[1]); swap(idAbs[0], idAbs[1]);}
  if (idAbs[1] > idAbs[2]) {swap( ip[1], ip[2]); swap(idAbs[1], idAbs[2]);}
  if (idAbs[0] > idAbs[1]) {swap( ip[0], ip[1]); swap(idAbs[0], idAbs[1]);}

  // Allow a few attempts to find a particle pair with low enough masses.
  for (int iTry = 0; iTry < nTry; ++iTry) {
    // Obtain the identities of two hadrons. Can fail, with id = 0.
    idHad = flavSelPtr->combineDiquarkJunction( flav1.id, flav2.id, flavj3.id);
    if (idHad.first == 0 || idHad.second == 0) mHadSum = mSum + 1.;
    // Check whether the mass sum fits inside the available phase space.
    else {
      mHad1 = particleDataPtr->mSel(idHad.first);
      mHad2 = particleDataPtr->mSel(idHad.second);
      mHadSum = mHad1 + mHad2;
    }
    if (mHadSum < mSum) break;
  }
  // Give up if all trials fail.
  if (mHadSum >= mSum) return false;

  // Define momentum of each hadron. Treat it as a -> b + c particle decay.
  pair<Vec4, Vec4> ps = rndmPtr->phaseSpace2(mSum, mHad1, mHad2);
  ps.first.bst(pSum, mSum);
  ps.second.bst(pSum, mSum);

  // Add produced particles to the event record. Note that the mothers
  // will be overwritten later. Note also the new status code 89.
  int iHad1 = event.append(idHad.first, 89, ip[0], ip[1], 0, 0, 0, 0,
    ps.first, mHad1);
  int iHad2 = event.append(idHad.second, 89, ip[0], ip[2], 0, 0, 0, 0,
    ps.second, mHad2);

  // Set production vertex of the two hadrons.
  if (setVertices) {
    Vec4 vHad1, vHad2;

    // Weighted average of the quarks entering each hadron, so 2 for a diquark.
    // Logic closely matches StringFlav::combineDiquarkJunction.
    if (event[iParton[1]].hasVertex()) {
      // When all ends are diquarks: split first to form two baryons.
      if (idAbs[0] > 10) {
        vHad1 = (event[ip[0]].vProd() + 2. * event[ip[1]].vProd()) / 3.;
        vHad2 = (event[ip[0]].vProd() + 2. * event[ip[2]].vProd()) / 3.;

      // When last two are diquarks: split middle to form a meson + baryon.
      } else if (idAbs[1] > 10) {
        vHad1 = (event[ip[0]].vProd() + event[ip[1]].vProd()) / 2.;
        vHad2 = (event[ip[1]].vProd() + 2. * event[ip[2]].vProd()) / 3.;

      // Else only third a diquark: split it to form two mesons.
      } else {
        vHad1 = (event[ip[0]].vProd() + event[ip[2]].vProd()) / 2.;
        vHad2 = (event[ip[1]].vProd() + event[ip[2]].vProd()) / 2.;
      }
    }

    // Add some "fragmentation time" offset based on hadron momenta.
    // Cruder approach than in setHadronVertices, but good enough.
    double vRel = sqrtpos( pow2( mSum*mSum - mHad1*mHad1 - mHad2*mHad2)
      - pow2( 2. * mHad1 * mHad2) ) / (mSum*mSum);
    vHad1 += 0.5 * (vRel / kappaVtx) * (ps.first / mHad1) * FM2MM;
    vHad2 += 0.5 * (vRel / kappaVtx) * (ps.second / mHad2) * FM2MM;

    // Set hadron production vertices.
    event[iHad1].vProd( vHad1);
    event[iHad2].vProd( vHad2);
  }

  // Set lifetime of hadrons and return.
  event[iHad1].tau( event[iHad1].tau0() * rndmPtr->exp() );
  event[iHad2].tau( event[iHad2].tau0() * rndmPtr->exp() );
  return true;

}

//--------------------------------------------------------------------------

// Attempt to produce two particles from the ministring.
// Note that popcorn baryons are not allowed.

bool MiniStringFragmentation::ministring2two( int nTry, Event& event,
  bool findLowMass) {

  // Properties of the produced hadrons.
  int    iFront  = iParton.front();
  int    iBack   = iParton.back();
  int    idHad1  = 0;
  int    idHad2  = 0;
  double mHad1   = 0.;
  double mHad2   = 0.;
  double mHadSum = 0.;

  // Allow a few attempts to find a particle pair with low enough masses.
  for (int iTry = 0; iTry < nTry; ++iTry) {

    // For closed gluon loop need to pick an initial flavour.
    if (isClosed) do {
      int idStart = flavSelPtr->pickLightQ();
      FlavContainer flavStart(idStart, 1);
      flavStart = flavSelPtr->pick( flavStart, -1., 0., false);
      flav1 = flavSelPtr->pick( flavStart, -1., 0., false);
      flav2.anti(flav1);
    } while (flav1.id == 0 || flav1.nPop > 0);

    // Create a new q qbar flavour to form two hadrons.
    // Start from a diquark, if any.
    do {
      FlavContainer flav3 =
        (flav1.isDiquark() || (!flav2.isDiquark() && rndmPtr->flat() < 0.5) )
        ? flavSelPtr->pick( flav1, -1., 0., false)
        : flavSelPtr->pick( flav2, -1., 0., false).anti();
      if (findLowMass) {
        idHad1 = flavSelPtr->combineToLightest( flav1.id, flav3.id);
        idHad2 = flavSelPtr->combineToLightest( flav2.id, -flav3.id);
      } else {
        idHad1 = flavSelPtr->combine( flav1, flav3);
        idHad2 = flavSelPtr->combine( flav2, flav3.anti());
      }
    } while (idHad1 == 0 || idHad2 == 0);

    // Check whether the mass sum fits inside the available phase space.
    mHad1 = particleDataPtr->mSel(idHad1);
    mHad2 = particleDataPtr->mSel(idHad2);
    mHadSum = mHad1 + mHad2;
    if (mHadSum < mSum) break;
  }

  // If not enough mass to create baryon-antibaryon pair in diquark-antidiquark
  // system, force reconnect to mesonic topology.
  if (mHadSum >= mSum && flav1.isDiquark() && flav2.isDiquark()) {
    // Split up diquark into individual flavours.
    int idTmp1 = flav2.id / 1000;
    int idTmp2 = (flav2.id % 1000) / 100;
    int idTmp3 = flav1.id / 1000;
    int idTmp4 = (flav1.id % 1000) / 100;

    // If findLowMass, select smallest mSum pairing, otherwise pick random one.
    int idHad13, idHad14, idHad23, idHad24;
    do {
      if (findLowMass) {
        idHad13 = flavSelPtr->combineToLightest( idTmp1, idTmp3);
        idHad14 = flavSelPtr->combineToLightest( idTmp1, idTmp4);
        idHad23 = flavSelPtr->combineToLightest( idTmp2, idTmp3);
        idHad24 = flavSelPtr->combineToLightest( idTmp2, idTmp4);
      } else {
        FlavContainer flavTmp1(idTmp1), flavTmp2(idTmp2), flavTmp3(idTmp3),
          flavTmp4(idTmp4);
        idHad13 = flavSelPtr->combine( flavTmp1, flavTmp3 );
        idHad14 = flavSelPtr->combine( flavTmp1, flavTmp4);
        idHad23 = flavSelPtr->combine( flavTmp2, flavTmp3);
        idHad24 = flavSelPtr->combine( flavTmp2, flavTmp4);
      }
      double mHad13 = particleDataPtr->mSel(idHad13);
      double mHad14 = particleDataPtr->mSel(idHad14);
      double mHad23 = particleDataPtr->mSel(idHad23);
      double mHad24 = particleDataPtr->mSel(idHad24);
      if ( ( findLowMass && mHad13 + mHad24 < mHad14 + mHad23)
        || (!findLowMass && rndmPtr->flat() > 0.5 ) ) {
        idHad1 = idHad13;
        idHad2 = idHad24;
        mHad1  = mHad13;
        mHad2  = mHad24;
      } else {
        idHad1 = idHad14;
        idHad2 = idHad23;
        mHad1  = mHad14;
        mHad2  = mHad23;
      }
    } while (idHad1 == 0 || idHad2 == 0);
    // Randomise which is considered coming from + side and which from -.
    if (rndmPtr->flat() > 0.5) {
      swap(idHad1, idHad2);
      swap(mHad1, mHad2);
    }
    mHadSum = mHad1 + mHad2;
  }

  // As last resort keep original flavours and split off pi0. Else fail.
  if (mHadSum >= mSum && findLowMass && !isClosed) {
    idHad1 = flavSelPtr->combineToLightest( flav1.id, flav2.id);
    if (idHad1 == 0) return false;
    idHad2 = 111;
    mHad1 = particleDataPtr->mSel(idHad1);
    mHad2 = particleDataPtr->mSel(idHad2);
    mHadSum = mHad1 + mHad2;
  }
  if (mHadSum >= mSum) return false;

  // Define an effective two-parton string, by splitting intermediate
  // gluon momenta in proportion to their closeness to either endpoint.
  Vec4 pSum1 = event[ iFront ].p();
  Vec4 pSum2 = event[ iBack ].p();
  if (iParton.size() > 2) {
    Vec4 pEnd1 = pSum1;
    Vec4 pEnd2 = pSum2;
    Vec4 pEndSum = pEnd1 + pEnd2;
    for (int i = 1; i < int(iParton.size()) - 1 ; ++i) {
      Vec4 pNow = event[ iParton[i] ].p();
      double ratio = (pEnd2 * pNow) / (pEndSum * pNow);
      pSum1 += ratio * pNow;
      pSum2 += (1. - ratio) * pNow;
    }
  }

  // If split did not provide an axis then pick random axis to break tie.
  // (Needed for low-mass q-g-qbar with q-qbar perfectly parallel.)
  if (pSum1.mCalc() + pSum2.mCalc() > 0.999999 * mSum) {
    double cthe = 2. * rndmPtr->flat() - 1.;
    double sthe = sqrtpos(1. - cthe * cthe);
    double phi  = 2. * M_PI * rndmPtr->flat();
    Vec4 delta  = 0.5 * min( pSum1.e(), pSum2.e())
        * Vec4( sthe * sin(phi), sthe * cos(phi), cthe, 0.);
    pSum1 += delta;
    pSum2 -= delta;
    loggerPtr->WARNING_MSG("random axis needed to break tie");
  }

  // Set up a string region based on the two effective endpoints.
  StringRegion region;
  region.setUp( pSum1, pSum2, 0, 0);

  // Generate an isotropic decay in the ministring rest frame,
  // suppressed at large pT by a fragmentation pT Gaussian.
  double pAbs2 = 0.25 * ( pow2(m2Sum - mHad1*mHad1 - mHad2*mHad2)
    - pow2(2. * mHad1 * mHad2) ) / m2Sum;
  double pT2 = 0.;
  do {
    double cosTheta = rndmPtr->flat();
    pT2 = (1. - pow2(cosTheta)) * pAbs2;
  } while (pTSelPtr->suppressPT2(pT2) < rndmPtr->flat() );

  // Construct the forward-backward asymmetry of the two particles.
  double mT21 = mHad1*mHad1 + pT2;
  double mT22 = mHad2*mHad2 + pT2;
  double lambda = sqrtpos( pow2(m2Sum  - mT21 - mT22) - 4. * mT21 * mT22 );
  double probReverse = 1. / (1. + exp( min( 50., bLund * lambda) ) );

  // Construct kinematics, as viewed in the transverse rest frame.
  double xpz1 = 0.5 * lambda/ m2Sum;
  if (probReverse > rndmPtr->flat()) xpz1 = -xpz1;
  double xmDiff = (mT21 - mT22) / m2Sum;
  double xe1 = 0.5 * (1. + xmDiff);
  double xe2 = 0.5 * (1. - xmDiff );

  // Distribute pT isotropically in angle.
  double phi = 2. * M_PI * rndmPtr->flat();
  double pT  = sqrt(pT2);
  double px  = pT * cos(phi);
  double py  = pT * sin(phi);

  // Translate this into kinematics in the string frame.
  Vec4 pHad1 = region.pHad( xe1 + xpz1, xe1 - xpz1,  px,  py);
  Vec4 pHad2 = region.pHad( xe2 - xpz1, xe2 + xpz1, -px, -py);

  // Mark hadrons from junction fragmentation with different status.
  int statusHadPos = 82, statusHadNeg = 82;
  if (abs(idHad1) > 1000 && abs(idHad1) < 10000 &&
      abs(idHad2) > 1000 && abs(idHad2) < 10000) {
    if (event[iFront].statusAbs() == 74) statusHadPos = 89;
    if (event[iBack ].statusAbs() == 74) statusHadNeg = 89;
  }
  else if (abs(idHad1) > 1000 && abs(idHad1) < 10000) {
    if ( event[iFront].statusAbs() == 74
      || event[iBack ].statusAbs() == 74) statusHadPos = 89;
  }
  else if (abs(idHad2) > 1000 && abs(idHad2) < 10000) {
    if ( event[iFront].statusAbs() == 74
      || event[iBack ].statusAbs() == 74) statusHadNeg = 89;
  }

  // Save parton vertices. Remove diquark if it is from a junction system.
  Vec4 vProdF = event[iFront].vProd();
  Vec4 vProdL = event[iBack ].vProd();
  if (isJunctionSystem) event.popBack(1);

  // Add produced particles to the event record.
  int iFirst = event.append( idHad1, statusHadPos, iFront, iBack,
    0, 0, 0, 0, pHad1, mHad1);
  int iLast  = event.append( idHad2, statusHadNeg, iFront, iBack,
    0, 0, 0, 0, pHad2, mHad2);

  // Set hadron vertices when they are displaced.
  if (event[iFront].hasVertex()) {
    event[iFirst].vProd( vProdF );
    event[iLast ].vProd( vProdL );
  }

  // Set lifetime of hadrons.
  event[iFirst].tau( event[iFirst].tau0() * rndmPtr->exp() );
  event[iLast ].tau( event[iLast ].tau0() * rndmPtr->exp() );

  // Mark original partons as hadronized and set their daughter range.
  if (!isJunctionSystem) {
    for (int i = 0; i < int(iParton.size()); ++i) {
      event[ iParton[i] ].statusNeg();
      event[ iParton[i] ].daughters(iFirst, iLast);
    }
  }

  // Store breakup vertex information from the fragmentation process.
  if (setVertices) {
    ministringVertices.clear();
    ministringVertices.push_back( StringVertex(true, 0, 0, 1., 0.) );
    ministringVertices.push_back(
      StringVertex(true, 0, 0, 1. - (xe1 + xpz1), xe1 - xpz1) );
    ministringVertices.push_back( StringVertex(true, 0, 0, 0., 1.) );

    // Store hadron production space-time vertices.
    setHadronVertices( event, region, iFirst, iLast);
  }

  // Successfully done.
  return true;

}

//--------------------------------------------------------------------------

// Attempt to produce one particle from a ministring.
// Current algorithm: find the system with largest invariant mass
// relative to the existing one, and boost that system appropriately.
// Optionally force pick lightest hadron possible to increase chances.
// Optionally let existing hadron take recoil instead of system.
// Try more sophisticated alternatives later?? (Z0 mass shifted??)

bool MiniStringFragmentation::ministring2one( int iSub,
  ColConfig& colConfig, Event& event, bool findLowMass, bool systemRecoil) {

  // Cannot handle qq + qbarqbar system.
  if (abs(flav1.id) > 100 && abs(flav2.id) > 100) return false;

  // For closed gluon loop need to pick an initial flavour.
  if (isClosed) do {
    int idStart = flavSelPtr->pickLightQ();
    FlavContainer flavStart(idStart, 1);
    flav1 = flavSelPtr->pick( flavStart);
    flav2 = flav1.anti();
  } while (abs(flav1.id) > 100);

  // Select hadron flavour from available quark flavours,
  // optionally with lowest possible mass.
  int idHad = 0;
  if (findLowMass) idHad = flavSelPtr->combineToLightest( flav1.id, flav2.id);
  else for (int iTryFlav = 0; iTryFlav < NTRYFLAV; ++iTryFlav) {
    idHad = flavSelPtr->combine( flav1, flav2);
    if (idHad != 0) break;
  }
  if (idHad == 0) return false;

  // Find mass.
  double mHad = particleDataPtr->mSel(idHad);

  // Find the untreated parton system, alternatively final hadron,
  // which combines to the largest squared mass above mimimum required.
  int iMax = -1;
  double deltaM2 = mHad*mHad - mSum*mSum;
  double delta2Max = 0.;
  if (systemRecoil) {
    for (int iRec = iSub + 1; iRec < colConfig.size(); ++iRec) {
      double delta2Rec = 2. * (pSum * colConfig[iRec].pSum) - deltaM2
        - 2. * mHad * colConfig[iRec].mass;
      if (delta2Rec > delta2Max) { iMax = iRec; delta2Max = delta2Rec;}
    }
  } else {
    for (int iRec = 0; iRec < event.size(); ++iRec)
    if (event[iRec].isHadron() && event[iRec].status() > 80) {
      double delta2Rec = 2. * (pSum * event[iRec].p()) - deltaM2
        - 2. * mHad * event[iRec].m();
      if (delta2Rec > delta2Max) { iMax = iRec; delta2Max = delta2Rec;}
    }
  }
  if (iMax == -1) return false;

  // Construct kinematics of the hadron and recoiling system (or hadron).
  Vec4   pRec    = (systemRecoil) ? colConfig[iMax].pSum : event[iMax].p();
  double mRec    = (systemRecoil) ? colConfig[iMax].mass : event[iMax].m();
  double vecProd = pSum * pRec;
  double coefOld = mSum*mSum + vecProd;
  double coefNew = mHad*mHad + vecProd;
  double coefRec = mRec*mRec + vecProd;
  double coefSum = coefOld + coefNew;
  double sHat    = coefOld + coefRec;
  double root    = sqrtpos( (pow2(coefSum) - 4. * sHat * mHad*mHad)
    / (pow2(vecProd) - pow2(mSum * mRec)) );
  double k2      = 0.5 * (coefOld * root - coefSum) / sHat;
  double k1      = (coefRec * k2 + 0.5 * deltaM2) / coefOld;
  Vec4 pHad      = (1. + k1) * pSum - k2 * pRec;
  Vec4 pRecNew   = (1. + k2) * pRec - k1 * pSum;

  // Mark hadrons from junction split off with status 89.
  int statusHad = 81;
  if (abs(idHad) > 1000 && abs(idHad) < 10000 &&
      (event[ iParton.front() ].statusAbs() == 74 ||
       event[ iParton.back() ].statusAbs() == 74)) statusHad = 89;

  // Remove diquark if it is from a junction system.
  if (isJunctionSystem) event.popBack(1);

  // Add the produced particle to the event record.
  int iHad = event.append( idHad, statusHad, iParton.front(), iParton.back(),
    0, 0, 0, 0, pHad, mHad);

  // Set production vertex when this is displaced.
  if (event[iParton.front()].hasVertex()) {
    Vec4 vProd = 0.5 * (event[iParton.front()].vProd()
      + event[iParton.back()].vProd());
    event[iHad].vProd( vProd );
  }

  // Set lifetime of hadron.
  event[iHad].tau( event[iHad].tau0() * rndmPtr->exp() );

  // Mark original partons as hadronized and set their daughter range.
  if (!isJunctionSystem) {
    for (int i = 0; i < int(iParton.size()); ++i) {
      event[ iParton[i] ].statusNeg();
      event[ iParton[i] ].daughters(iHad, iHad);
    }
  }

  // Copy down recoiling system, with boosted momentum. Update current partons.
  RotBstMatrix M;
  M.bst(pRec, pRecNew);
  if (systemRecoil) {
    for (int i = 0; i < colConfig[iMax].size(); ++i) {
      int iOld = colConfig[iMax].iParton[i];
      // Do not touch negative iOld = beginning of new junction leg.
      if (iOld >= 0) {
        int iNew;
        // Keep track of 74 throughout the event.
        if (event[iOld].status() == 74) iNew = event.copy(iOld, 74);
        else iNew = event.copy(iOld, 72);
        event[iNew].rotbst(M);
        colConfig[iMax].iParton[i] = iNew;
      }
    }
    colConfig[iMax].pSum = pRecNew;
    colConfig[iMax].isCollected = true;

  // Alternatively copy down modified hadron and boost it (including vertex).
  } else {
    int iNew = event.copy(iMax, event[iMax].status());
    event[iNew].rotbst(M);
  }

  // Calculate hadron production points from breakup vertices
  // using one of the three definitions.
  if (setVertices) {
    Vec4 prodPoint = Vec4( 0., 0., 0., 0.);
    Vec4 pHadron = event[iHad].p();

    // Smearing in transverse space.
    if (smearOn) {

      // Find two spacelike transverse four-vector directions.
      Vec4 eX = Vec4( 1., 0., 0., 0.);
      Vec4 eY = Vec4( 0., 1., 0., 0.);

      // Introduce smearing in transverse space.
      double transX = rndmPtr -> gauss();
      double transY = rndmPtr -> gauss();
      prodPoint = xySmear * (transX * eX + transY * eY) / sqrt(2.);
      // Keep proper or actual time constant when including the smearing.
      // Latter case to be done better when introducing MPI vertices.
      if (constantTau) prodPoint.e( prodPoint.pAbs() );
      else prodPoint = Vec4( 0., 0., 0., 0.);
    }

    // Reduced oscillation period if hadron contains massive quarks.
    int id1 = event[ iParton.front() ].idAbs();
    int id2 = event[ iParton.back() ].idAbs();
    double redOsc = 1.;
    if (id1 == 4 || id1 == 5 || id2 == 4 || id2 == 5) {
      double posMass = (id1 == 4 || id1 == 5) ? particleDataPtr->m0(id1) : 0.;
      double negMass = (id2 == 4 || id2 == 5) ? particleDataPtr->m0(id2) : 0.;
      redOsc = sqrtpos( pow2(pow2(mHad) - pow2(posMass) - pow2(negMass))
        - 4. * pow2(posMass * negMass) ) / pow2(mHad);
    }

    // Find hadron production points according to chosen definition.
    if (hadronVertex == 0) prodPoint += 0.5 * redOsc * pHadron / kappaVtx;
    else if (hadronVertex == 1) prodPoint += redOsc * pHadron / kappaVtx;
    event[iHad].vProd( event[iHad].vProd() + prodPoint * FM2MM );
  }

  // Successfully done.
  return true;

}

//--------------------------------------------------------------------------

// Store two hadron production points in the event record.

void MiniStringFragmentation::setHadronVertices(Event& event,
  StringRegion& region, int iFirst, int iLast) {

  // Initial values.
  vector<Vec4> longitudinal;
  int id1 = event[ iParton.front() ].idAbs();
  int id2 = event[ iParton.back() ].idAbs();

  // Longitudinal space-time location of breakup points.
  for (int i = 0; i < 3; ++i) {
    double xPosIn = ministringVertices[i].xRegPos;
    double xNegIn = ministringVertices[i].xRegNeg;
    Vec4 noOffset = (xPosIn * region.pPos + xNegIn * region.pNeg) / kappaVtx;
    longitudinal.push_back( noOffset );
  }

  // Longitudinal offset of breakup points for massive quarks.
  if (region.massiveOffset( 0, 0, 0, id1, id2, mc, mb)) {
    for (int i = 0; i < 3; ++i) {

      // Endpoint correction separately for each end.
      if (i == 0 && (id1 == 4 || id1 == 5)) {
        Vec4 v1 = longitudinal[i];
        Vec4 v2 = longitudinal[i + 1];
        double mHad =  event[event.size() - 2].m();
        double pPosMass = particleDataPtr->m0(id1);
        longitudinal[i] = v1 + (pPosMass / mHad) * (v2 - v1);
      }
      if (i == 2 && (id2 == 4 || id2== 5)) {
        Vec4 v1 = longitudinal[i];
        Vec4 v2 = longitudinal[i-1] + region.massOffset / kappaVtx;
        double mHad =  event[i - 1 + event.size() - 2].m();
        double pNegMass = particleDataPtr->m0(id2);
        longitudinal[i] = v1 + (pNegMass / mHad) * (v2 - v1);
        if (longitudinal[i].m2Calc()
           < -1e-8 * max(1., pow2(longitudinal[i].e())))
           loggerPtr->WARNING_MSG(
             "negative tau^2 for endpoint massive correction");
      }

      // Add mass offset for all breakup points.
      Vec4 massOffset = region.massOffset / kappaVtx;
      Vec4 position = longitudinal[i] - massOffset;

      // Correction for non-physical situations.
      if (position.m2Calc() < 0.) {
        double cMinus = 0.;
        if (position.m2Calc() > -1e-8 * max(1., pow2(position.e())))
          position.e( position.pAbs() );
        else {
          if(massOffset.m2Calc() > 1e-6)
            cMinus = (longitudinal[i] * massOffset
              - sqrt(pow2(longitudinal[i] * massOffset)
              - longitudinal[i].m2Calc() * massOffset.m2Calc()))
              / massOffset.m2Calc();
          else cMinus = 0.5 * longitudinal[i].m2Calc()
              / (longitudinal[i] * massOffset);
          position = longitudinal[i] - cMinus * massOffset;
        }
      }
      longitudinal[i] = position;
    }
  }

  // Smearing in transverse space.
  vector<Vec4> spaceTime;
  for (int i = 0; i < 3; ++i) {
    Vec4 positionTot =  longitudinal[i];
    if (smearOn) {

      if (!isClosed && (i == 0 || i == 2)) {
        spaceTime.push_back(positionTot);
        continue;
      }
      Vec4 eX = region.eX;
      Vec4 eY = region.eY;

      // Smearing calculated randomly following a gaussian.
      for (int iTry = 0; ; ++iTry) {
        double transX = rndmPtr->gauss();
        double transY = rndmPtr->gauss();
        Vec4 transversePos = xySmear * (transX * eX + transY * eY) / sqrt(2.);
        positionTot = transversePos + longitudinal[i];

        // Keep proper or actual time constant when including the smearing.
        // Latter case to be done better when introducing MPI vertices.
        if (constantTau) {
          double newtime = sqrt(longitudinal[i].m2Calc()
            + positionTot.pAbs2());
          positionTot.e(newtime);
          break;
        } else {
          if (positionTot.m2Calc() >= 0.) break;
          if (iTry == 100) {
            positionTot = longitudinal[i];
            break;
          }
        }
      }
    }
    spaceTime.push_back(positionTot);
  }

  // Find hadron production points according to chosen definition.
  vector<Vec4> prodPoints(2);
  for(int i = 0; i < 2; ++i) {
    Vec4 middlePoint = 0.5 * (spaceTime[i] + spaceTime[i+1]);
    int  iHad = (i == 0) ? iFirst : iLast;
    Vec4 pHad = event[iHad].p();

    // Reduced oscillation period if hadron contains massive quarks.
    double mHad = event[iHad].m();
    int    idQ  = (i == 0) ? id1 : id2;
    double redOsc = (idQ == 4 || idQ == 5)
      ? 1. - pow2(particleDataPtr->m0(idQ) / mHad) : 0.;

    // Set production point according to chosen definition.
    if (hadronVertex == 0) prodPoints[i] = middlePoint;
    else if (hadronVertex == 1)
      prodPoints[i] = middlePoint + 0.5 * redOsc * pHad / kappaVtx;
    else {
      prodPoints[i] = middlePoint - 0.5 * redOsc * pHad / kappaVtx;
      if (prodPoints[i].m2Calc() < 0. || prodPoints[i].e() < 0.) {
        double tau0fac = 2. * (redOsc * middlePoint * pHad
          - sqrt(pow2(middlePoint * redOsc * pHad) - middlePoint.m2Calc()
          * pow2(redOsc * mHad))) / pow2(redOsc * mHad);
        prodPoints[i] = middlePoint - 0.5 * tau0fac * redOsc * pHad / kappaVtx;
      }
    }
    event[iHad].vProd( event[iHad].vProd() + prodPoints[i] * FM2MM );
  }

}

//==========================================================================

} // end namespace Pythia8
