// LowEnergyProcess.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the LowEnergyProcess
// class.

#include "Pythia8/LowEnergyProcess.h"

namespace Pythia8 {

//==========================================================================

// LowEnergyProcess class.
// This class handles low-energy collisions between two hadrons.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to split beam particles before reconnection.
const int LowEnergyProcess::MAXLOOP = 100;

// Gradually reduce assumed quark masses from their constituent values.
const double LowEnergyProcess::MASSREDUCERATE = 0.02;

// Let diffractive mass spectrum begin this much above unexcited mass.
const double LowEnergyProcess::MDIFFMIN = 0.2;

// Pomeron trajectory alpha(t) = 1 + epsilon + alpha' * t
const double LowEnergyProcess::ALPHAPRIME = 0.25;

//--------------------------------------------------------------------------

// Initialize the LowEnergyProcess class as required.

bool LowEnergyProcess::init(StringFragmentation* stringFragPtrIn,
  MiniStringFragmentation* ministringFragPtrIn) {

  // Save pointers.
  stringFragPtr     = stringFragPtrIn;
  ministringFragPtr = ministringFragPtrIn;

  // Relative fraction of s quark production in strin breaks.
  probStoUD       = parm("StringFlav:probStoUD");

  // Mixing for eta and eta'.
  double theta    = parm("StringFlav:thetaPS");
  double alpha    = (theta + 54.7) * M_PI / 180.;
  fracEtass       = pow2(sin(alpha));
  fracEtaPss      = 1. - fracEtass;

  // Longitudinal momentum sharing of valence quarks in hadrons.
  xPowMes         = parm("BeamRemnants:valencePowerMeson");
  xPowBar         = 0.5 * ( parm("BeamRemnants:valencePowerUinP")
                          + parm("BeamRemnants:valencePowerDinP") );
  xDiqEnhance     = parm("BeamRemnants:valenceDiqEnhance");

  // Transverse momentum spread.
  sigmaQ          = parm("StringPT:sigma") / sqrt(2.);

  // Boundary mass between string and ministring handling.
  mStringMin      = parm("HadronLevel:mStringMin");

  // Initialize collision event record.
  leEvent.init( "(low energy event)", particleDataPtr);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Produce outgoing primary hadrons from collision of incoming pair.
// type = 0: mix; = 1: nondiff; = 2 : el; = 3: SD (XB); = 4: SD (AX);
//      = 5: DD; = 6: annihilation.

bool LowEnergyProcess::collide( int i1, int i2, int typeIn, Event& event,
  Vec4 vtx) {

  // Check that incoming hadrons. Store current event size.
  if (!event[i1].isHadron() || !event[i2].isHadron()) return false;
  if (typeIn < 0 || typeIn > 6) return false;
  sizeOld = event.size();

  // Pick event type for typeIn = 0.
  int type = typeIn;
  if (typeIn == 0) {
    // To be done?? Requires access to partial cross sections.
  }

  //  Hadron type and meson/baryon distinction.
  id1       = event[i1].id();
  id2       = event[i2].id();
  isBaryon1 = ( (abs(id1)/1000)%10 > 0 );
  isBaryon2 = ( (abs(id2)/1000)%10 > 0 );

  //  Hadron masses and collision invariant mass.
  m1        = event[i1].m();
  m2        = event[i2].m();
  eCM       = (event[i1].p() + event[i2].p()).mCalc();
  sCM       = eCM * eCM;

  // Reset leEvent event record. Add incoming hadrons as beams in rest frame.
  leEvent.reset();
  leEvent.append( event[i1]);
  leEvent.append( event[i2]);
  leEvent[1].status( -12);
  leEvent[2].status( -12);
  RotBstMatrix MtoCM = toCMframe( leEvent[1].p(), leEvent[2].p());
  leEvent.rotbst( MtoCM);

  // Do inelastic nondiffractive collision.
  if (type == 1 && !nondiff()) return false;

  // Do elastic or diffractive collision.
  if (type > 1 && type < 6 && !eldiff( type)) return false;

  // Do annihilation collision.
  if (type == 6 && !annihilation()) return false;

  // Hadronize new strings and move products to standard event record.
  if (!simpleHadronization( leEvent)) {
    infoPtr->errorMsg( "Error in LowEnergyProcess::collide: "
      "no rescattering since hadronization failed");
    return false;
  }
  for (int i = 3; i < leEvent.size(); ++i) if (leEvent[i].isFinal())
    event.append( leEvent[i]);

  // Boost from collision rest frame to event frame.
  // Set status and mothers. Offset vertex info to collision vertex.
  RotBstMatrix MfromCM = fromCMframe( event[i1].p(), event[i2].p());
  int mother1 = min(i1, i2);
  int mother2 = max(i1, i2);
  for (int i = sizeOld; i < event.size(); ++i) {
    event[i].rotbst( MfromCM);
    event[i].status( 110 + type);
    event[i].mothers( mother2, mother1 );
    event[i].vProdAdd( vtx);
  }

  // Mark incoming colliding hadrons as decayed. Set their daughters.
  event[i1].statusNeg();
  event[i2].statusNeg();
  event[i1].daughters( sizeOld, event.size() - 1);
  event[i2].daughters( sizeOld, event.size() - 1);

  // Add option to do decays?? Not to be done here for space-time tracing.

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do an inelastic nondiffractive scattering.

bool LowEnergyProcess::nondiff() {

  // Check that not stuck in infinite loop. Allow reduced quark masses.
  int    loop = 0;
  double mAbove1, mAbove2;
  Vec4   pc1, pac1, pc2, pac2;
  do {
    do {
      if (++loop == MAXLOOP) {
        infoPtr->errorMsg("Error in LowEnergyProcess::nondiff: "
          " failed to construct valid kinematics");
        return false;
      }
      double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9));

      // Split up hadrons A  and B into q + qbar or q + qq for meson/baryon.
      splitA( redStep);
      splitB( redStep);

      // Assign relative sharing of longitudinal momentum.
      z1     = splitZ( idc1, idac1, mTc1 / eCM, mTac1 / eCM);
      z2     = splitZ( idc2, idac2, mTc2 / eCM, mTac2 / eCM);
      mT1    = sqrt( mTsc1 / z1 + mTsac1 / (1. - z1));
      mT2    = sqrt( mTsc2 / z2 + mTsac2 / (1. - z2));

    // Ensure that hadron beam remnants are not too massive.
    } while (mT1 + mT2 > eCM);

    // Set up kinematics for outgoing beam remnants.
    double e1    = 0.5 * (sCM + mT1 * mT1 - mT2 * mT2) / eCM;
    double pz1   = sqrtpos(e1 * e1 - mT1 * mT1);
    double epz1  = z1 * (e1 + pz1);
    double pzc1  = 0.5 * (epz1 - mTsc1 / epz1 );
    double ec1   = 0.5 * (epz1 + mTsc1 / epz1 );
    pc1.p(   px1,  py1,       pzc1,      ec1 );
    pac1.p( -px1, -py1, pz1 - pzc1, e1 - ec1 );
    double epz2  = z2 * (eCM - e1 + pz1);
    double pzc2  = -0.5 * (epz2 - mTsc2 / epz2 );
    double ec2   =  0.5 * (epz2 + mTsc2 / epz2 );
    pc2.p(   px2,  py2,        pzc2,            ec2 );
    pac2.p( -px2, -py2, -pz1 - pzc2, eCM - e1 - ec2 );

    // Catch reconnected systems with too small masses.
    mAbove1 = (pc1 + pac2).mCalc() - mThreshold( idc1, idac2);
    mAbove2 = (pc2 + pac1).mCalc() - mThreshold( idc2, idac1);
  } while (mAbove1 < 0. || mAbove2 < 0.);

  // Store new reconnected string systems.
  leEvent.append(  idc1, 63, 1, 0, 0, 0, 101,   0,  pc1,  mc1);
  leEvent.append( idac2, 63, 2, 0, 0, 0,   0, 101, pac2, mac2);
  leEvent.append(  idc2, 63, 2, 0, 0, 0, 102,   0,  pc2,  mc2);
  leEvent.append( idac1, 63, 1, 0, 0, 0,   0, 102, pac1, mac1);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do an elastic or diffractive scattering.
// type = 2: elastic; = 3: SD (XB); = 4: SD (AX); = 5: DD.

bool LowEnergyProcess::eldiff( int type) {

  // Classify process type.
  bool excite1 = (type == 3 || type == 5);
  bool excite2 = (type == 4 || type == 5);

  // Find excited mass ranges.
  mA           = m1;
  mB           = m2;
  double mAmin = (excite1) ? m1 + MDIFFMIN : m1;
  double mBmin = (excite2) ? m2 + MDIFFMIN : m2;
  double mAmax = eCM - mBmin;
  double mBmax = eCM - mAmin;
  if (mAmin + mBmin > eCM) {
    infoPtr->errorMsg("Error in LowEnergyProcess::eldiff: "
      " too low invariant mass for diffraction");
    return false;
  }

  // Check that not stuck in infinite loop. Allow reduced quark masses.
  int  loop  = 0;
  bool failM = false;
  do {
    failM = false;
    if (++loop == MAXLOOP) {
      infoPtr->errorMsg("Error in LowEnergyProcess::eldiff: "
        " failed to construct valid kinematics");
      return false;
    }
    double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9));

    // Split up hadron 1 (on side A) and assign excited A mass.
    // Should contain better low-mass description??
    if (excite1) {
      splitA( redStep);
      mA     = mAmin * pow( mAmax / mAmin, rndmPtr->flat() );
      if (mTc1 + mTac1 > mA) failM = true;
    }

    // Split up hadron 2 (on side B) and assign excited B mass.
    // Should contain better low-mass description??
    if (excite2 && !failM) {
      splitB( redStep);
      mB     = mBmin * pow( mBmax / mBmin, rndmPtr->flat() );
      if (mTc2 + mTac2 > mB) failM = true;
    }

  // Ensure that pair of hadron masses not too large.
  if (mA + mB > eCM) failM = true;
  } while (failM);

  // Squared masses, energies and longitudinal momenta of excited hadrons.
  double s1    = m1 * m1;
  double s2    = m2 * m2;
  double sA    = mA * mA;
  double sB    = mB * mB;
  double eA    = 0.5 * (sCM + sA - sB) / eCM;
  double pzA   = sqrtpos(eA * eA - sA);
  Vec4   pA( 0., 0.,  pzA,       eA);
  Vec4   pB( 0., 0., -pzA, eCM - eA);

  // Internal kinematics on side A, boost to CM frame and store constituents.
  if (excite1) {
    double ec1   = 0.5 * (sA + mTsc1 - mTsac1) / mA;
    double pzc1  = sqrtpos(ec1 * ec1 - mTsc1);
    // Diquark always forward. Randomize for meson.
    if ( abs(idac1) > 10 || (abs(idc1) < 10 && abs(idac1) < 10
      && rndmPtr->flat() > 0.5) ) pzc1 = -pzc1;
    Vec4 pc1(   px1,  py1,  pzc1,      ec1);
    Vec4 pac1( -px1, -py1, -pzc1, mA - ec1);
    pc1.bst(pA);
    pac1.bst(pA);
    leEvent.append(  idc1, 63, 1, 0, 0, 0, 101,   0,  pc1,  mc1);
    leEvent.append( idac1, 63, 1, 0, 0, 0,   0, 101, pac1, mac1);

  // Simple copy if not excited, and set momentum as in collision frame.
  } else {
    int iNew = leEvent.copy( 1, 63);
    leEvent[iNew].p( pA);
    leEvent[iNew].vProd( 0., 0., 0., 0.);
  }

  // Internal kinematics on side B, boost to CM frame and store constituents.
  if (excite2) {
    double ec2   = 0.5 * (sB + mTsc2 - mTsac2) / mB;
    double pzc2  = -sqrtpos(ec2 * ec2 - mTsc2);
    // Diquark always forward (on negative side). Randomize for meson.
    if ( abs(idac2) > 10 || (abs(idc2) < 10 && abs(idac2) < 10
      && rndmPtr->flat() > 0.5) ) pzc2 = -pzc2;
    Vec4 pc2(   px2,  py2,  pzc2,      ec2);
    Vec4 pac2( -px2, -py2, -pzc2, mB - ec2);
    pc2.bst(pB);
    pac2.bst(pB);
    leEvent.append(  idc2, 63, 2, 0, 0, 0, 102,   0,  pc2,  mc2);
    leEvent.append( idac2, 63, 2, 0, 0, 0, 0,   102, pac2, mac2);

  // Simple copy if not excited, and set momentum as in collision frame.
  } else {
    int iNew = leEvent.copy( 2, 63);
    leEvent[iNew].p( pB);
    leEvent[iNew].vProd( 0., 0., 0., 0.);
  }

  // Select t value and rotate outgoing particles accordingly.
  double lambda12 = pow2( sCM - s1 - s2) - 4. * s1 * s2;
  double lambdaAB = pow2( sCM - sA - sB) - 4. * sA * sB;
  double tLow     = -0.5 * (sCM - (s1 + s2 + sA + sB) + (s1 - s2)
    * (sA - sB) / sCM + sqrtpos(lambda12 *  lambdaAB) / sCM);
  double tUpp     = ( (sA - s1) * (sB - s2) + (s1 + sB - s2 - sA)
    * (s1 * sB - s2 * sA) / sCM ) / tLow;
  double bNow     = bSlope( type);
  double eBtLow   = exp( bNow * tLow);
  double eBtUpp   = exp( bNow * tUpp);
  double tNow     = log( eBtLow + rndmPtr->flat() * (eBtUpp - eBtLow) ) / bNow;
  double theta    = acos( (2. * tNow - tLow - tUpp) / (tUpp - tLow) );
  double phi      = 2. * M_PI * rndmPtr->flat();
  for (int i = 3; i < leEvent.size(); ++i) leEvent[i].rot( theta, phi);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Do an annihilation collision.
// Unsolved: how handle K0S and K0L, which have alternating flavours??

bool LowEnergyProcess::annihilation() {

  // Vectors of quarks ("ord") and antiquarks ("bar"), and of their pairs.
  vector<int> iqord, iqbar, iord, ibar;

  // Split first and second hadron by flavour content.
  for (int iHad = 0; iHad < 2; ++iHad) {
    pair< int, int>  paircac  = splitFlav( (iHad == 0) ? id1 : id2 );
    int idcAbs   = abs(paircac.first);
    int idacAbs  = abs(paircac.second);
    if (idcAbs == 0 || idacAbs == 0) return false;

    // Store flavour content for later usage.
    if (iHad == 0) {
      idc1  = paircac.first;
      idac1 = paircac.second;
    } else {
      idc2  = paircac.first;
      idac2 = paircac.second;
    }

    // Store flavour content, with a diquark split further.
    if (idcAbs < 10) iqord.push_back( idcAbs );
    else {
      iqbar.push_back( (idcAbs/1000)%10 );
      iqbar.push_back( (idcAbs/100)%10 );
    }
    if (idacAbs < 10) iqbar.push_back( idacAbs );
    else {
      iqord.push_back( (idacAbs/1000)%10 );
      iqord.push_back( (idacAbs/100)%10 );
    }
  }

  // Find potential annihilating quark-antiquark pairs.
  for (int i1 = 0; i1 < int(iqord.size()); ++i1)
  for (int i2 = 0; i2 < int(iqbar.size()); ++i2)
  if (iqbar[i2] == iqord[i1]) {
    iord.push_back(i1);
    ibar.push_back(i2);
  }

  // Return if no annihilation possible.
  if (iord.size() == 0) {
    infoPtr->errorMsg( "Warning in LowEnergyProcess::annihilation: "
      "flavour content does not allow annihilation");
    return false;
  }

  // Annihilate one quark-antiquark pair at random among options.
  int iAnn = max( 0, min( int(iord.size()) - 1,
    int(iord.size() * rndmPtr->flat()) ));
  iqord[iord[iAnn]] = iqord.back();
  iqord.pop_back();
  iqbar[ibar[iAnn]] = iqbar.back();
  iqbar.pop_back();

  // Optionally allow baryon number annihilation - to be refined??
  if (iqord.size() == 2 && iqbar.size() == 2 && rndmPtr->flat() > 0.5) {
    iord.clear();
    ibar.clear();

    // Find potential annihilating second quark-antiquark pairs.
    for (int i1 = 0; i1 < 2; ++i1)
    for (int i2 = 0; i2 < 2; ++i2)
    if (iqbar[i2] == iqord[i1]) {
      iord.push_back(i1);
      ibar.push_back(i2);
    }

    // Annihilate a second quark-antiquark pair if possible.
    if (iord.size() > 0) {
      iAnn = max( 0, min( int(iord.size()) - 1,
        int(iord.size() * rndmPtr->flat()) ));
      iqord[iord[iAnn]] = iqord.back();
      iqord.pop_back();
      iqbar[ibar[iAnn]] = iqbar.back();
      iqbar.pop_back();
    }
  }

  // Optionally allow full annihilation of meson-meson or baryon-antibaryon,
  // and creation of new flavour - to be refined?
  if (iqord.size() == 1 && iqbar.size() == 1 && iqbar[0] == iqord[0]
    && rndmPtr->flat() > 0.5) {
    double rndmq = (2. + probStoUD) * rndmPtr->flat();
    iqord[0] = (rndmq < 1.) ? 1 : ((rndmq < 2.) ? 2 : 3);
    iqbar[0] = iqord[0];
  }

  // Recombine into diquarks where required; for now at random.
  int idcAnn  = 0;
  int idacAnn = 0;
  if (iqord.size() == 1 && iqbar.size() == 1) {
    idcAnn  = iqord[0];
    idacAnn = -iqbar[0];
  } else if (iqord.size() == 3 && iqbar.size() == 0) {
    int iPick = max( 0, min( 2, int(3. * rndmPtr->flat()) ));
    idcAnn  = iqord[iPick];
    int iq4 = iqord[(iPick + 1)%3];
    int iq5 = iqord[(iPick + 2)%3];
    idacAnn = 1000 * max(iq4, iq5) + 100 * min( iq4, iq5)
            + ( (iq5 == iq4) ? 3 : 1 );
  } else if (iqord.size() == 0 && iqbar.size() == 3) {
    int iPick = max( 0, min( 2, int(3. * rndmPtr->flat()) ));
    idacAnn = -iqbar[iPick];
    int iq4 = iqbar[(iPick + 1)%3];
    int iq5 = iqbar[(iPick + 2)%3];
    idcAnn  = -(1000 * max(iq4, iq5) + 100 * min( iq4, iq5)
            + ( (iq5 == iq4) ? 3 : 1 ));
  } else if (iqord.size() == 2 && iqbar.size() == 2) {
    idcAnn  = 1000 * max(iqbar[0], iqbar[1]) + 100 * min( iqbar[0], iqbar[1])
            + ( (iqbar[1] == iqbar[0]) ? 3 : 1 );
    idacAnn = -(1000 * max(iqord[0], iqord[1]) + 100 * min( iqord[0], iqord[1])
            + ( (iqord[1] == iqord[0]) ? 3 : 1 ));
  } else {
    infoPtr->errorMsg( "Error in LowEnergyProcess::: "
      "obtained unphysical flavour content");
    return false;
  }

  // Begin kinematics construction. Allow reduced quark masses.
  int    loop = 0;
  double mcAnn, macAnn, pxAnn, pyAnn, pTsAnn, mTscAnn, mTsacAnn;
  do {
    if (++loop == MAXLOOP) {
      infoPtr->errorMsg("Error in LowEnergyProcess::nondiff: "
        " failed to construct valid kinematics");
      return false;
    }
    double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9));

    // Find constituent masses and scale down if necessary.
    mcAnn  = redStep * particleDataPtr->m0( idcAnn);
    macAnn = redStep * particleDataPtr->m0( idacAnn);

    // Select Gaussian relative transverse momenta for constituents.
    pair<double, double> gauss2 = rndmPtr->gauss2();
    pxAnn  = redStep * sigmaQ * gauss2.first;
    pyAnn  = redStep * sigmaQ * gauss2.second;
    pTsAnn = pxAnn * pxAnn + pyAnn * pyAnn;

    // Construct transverse masses.
    mTscAnn  = pow2(mcAnn)  + pTsAnn;
    mTsacAnn = pow2(macAnn) + pTsAnn;

  // Ensure that hadron beam remnants are not too massive.
  } while (sqrt(mTscAnn) + sqrt(mTsacAnn) > eCM);

  // Decide which side is which.
  double sgnAnn = 1.;
  if (isBaryon1 && isBaryon2) {
    if (abs(idcAnn) > 10) sgnAnn = -1.;
    if (id1 < 0) sgnAnn *= -1.;
  } else if (isBaryon1) {
    if (id1 > 0) sgnAnn = -1.;
  } else if (isBaryon2) {
    if (id2 < 0) sgnAnn = -1.;
  } else {
    bool match12 = idcAnn == idc1 && idacAnn == idac2;
    bool match21 = idcAnn == idc2 && idacAnn == idac1;
    if (match12 && match21) sgnAnn = (rndmPtr->flat() > 0.5) ? 1. : -1.;
    else if (match12)       sgnAnn = 1.;
    else if (match21)       sgnAnn = -1.;
    else                    sgnAnn = (rndmPtr->flat() > 0.5) ? 1. : -1.;
  }

  // Set up kinematics for outgoing beam remnants.
  double ecAnn  = 0.5 * (sCM + mTscAnn - mTsacAnn) / eCM;
  double pzcAnn = sgnAnn * sqrtpos(ecAnn * ecAnn - mTscAnn);
  Vec4 pcAnn(   pxAnn,  pyAnn,  pzcAnn,       ecAnn );
  Vec4 pacAnn( -pxAnn, -pyAnn, -pzcAnn, eCM - ecAnn );

  // Store new single string system.
  leEvent.append(  idcAnn, 63, 1, 2, 0, 0, 101,   0,  pcAnn,  mcAnn);
  leEvent.append( idacAnn, 63, 1, 2, 0, 0,   0, 101, pacAnn, macAnn);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Simple version of hadronization for low-energy hadronic collisions.
// Only accepts simple q-qbar systems and hadrons.

bool LowEnergyProcess::simpleHadronization( Event& event, bool isDiff) {

  // Find the complete colour singlet configuration of the event.
  simpleColConfig.clear();
  for (int i = 0; i < event.size(); ++i)
  if (event[i].isQuark() || event[i].isDiquark()) {
    vector<int> qqPair;
    qqPair.push_back(   i);
    qqPair.push_back( ++i);
    simpleColConfig.simpleInsert( qqPair, event );
  }

  // Process all colour singlet (sub)systems.
  for (int iSub = 0; iSub < simpleColConfig.size(); ++iSub) {

    // String fragmentation of each colour singlet (sub)system.
    if ( simpleColConfig[iSub].massExcess > mStringMin ) {
      if (!stringFragPtr->fragment( iSub, simpleColConfig, event))
        return false;

    // Low-mass string treated separately. Tell if diffractive system.
    } else {
      if (!ministringFragPtr->fragment( iSub, simpleColConfig, event, isDiff))
        return false;
    }
  }

  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Split up hadron A into a colour-anticolour pair, with masses and pT values.

bool LowEnergyProcess::splitA( double redMpT) {

  // Split up flavour of hadron into a colour and an anticolour constituent.
  pair< int, int>  paircac  = splitFlav( id1 );
  idc1   = paircac.first;
  idac1  = paircac.second;
  if (idc1 == 0 || idac1 == 0) return false;

  // Find constituent masses and scale down to less than full mass.
  mc1    = particleDataPtr->m0( idc1);
  mac1   = particleDataPtr->m0( idac1);
  double redNow = redMpT * min( 1., m1 / (mc1 + mac1));
  mc1   *= redNow;
  mac1  *= redNow;

  // Select Gaussian relative transverse momenta for constituents.
  pair<double, double> gauss2 = rndmPtr->gauss2();
  px1    = redMpT * sigmaQ * gauss2.first;
  py1    = redMpT * sigmaQ * gauss2.second;
  pTs1   = px1 * px1 + py1 * py1;

  // Construct transverse masses.
  mTsc1  = pow2(mc1)  + pTs1;
  mTsac1 = pow2(mac1) + pTs1;
  mTc1   = sqrt(mTsc1);
  mTac1  = sqrt(mTsac1);

  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Split up hadron B into a colour-anticolour pair, with masses and pT values.

bool LowEnergyProcess::splitB( double redMpT) {

  // Split up flavour of hadron into a colour and an anticolour constituent.
  pair< int, int>  paircac  = splitFlav( id2 );
  idc2   = paircac.first;
  idac2  = paircac.second;
  if (idc2 == 0 || idac2 == 0) return false;

  // Find constituent masses and scale down to less than full mass.
  mc2    = particleDataPtr->m0( idc2);
  mac2   = particleDataPtr->m0( idac2);
  double redNow = redMpT * min( 1., m2 / (mc2 + mac2));
  mc2   *= redNow;
  mac2  *= redNow;

  // Select Gaussian relative transverse momenta for constituents.
  pair<double, double> gauss2 = rndmPtr->gauss2();
  px2    = redMpT * sigmaQ * gauss2.first;
  py2    = redMpT * sigmaQ * gauss2.second;
  pTs2   = px2 * px2 + py2 * py2;

  // Construct transverse masses.
  mTsc2  = pow2(mc2)  + pTs2;
  mTsac2 = pow2(mac2) + pTs2;
  mTc2   = sqrt(mTsc2);
  mTac2  = sqrt(mTsac2);

  // Done.
  return true;

}

//-------------------------------------------------------------------------

// Split up a hadron into a colour and an anticolour part, of q or qq kinds.

pair< int, int> LowEnergyProcess::splitFlav( int id) {

  // Hadron flavour content.
  int idAbs = abs(id);
  int iq1   = (idAbs/1000)%10;
  int iq2   = (idAbs/100)%10;
  int iq3   = (idAbs/10)%10;
  int iq4, iq5;

  // Nondiagonal mesons.
  if (iq1 == 0 && iq2 != iq3) {
    if (id != 130 && id != 310) {
      if (iq2%2 == 1) swap( iq2, iq3);
      if (id > 0) return make_pair( iq2, -iq3);
      else        return make_pair( iq3, -iq2);

    // K0S and K0L are mixes d sbar and dbar s.
    } else {
      if (rndmPtr->flat() < 0.5) return make_pair( 3, -1);
      else                       return make_pair( 1, -3);
    }

  // Diagonal mesons: assume complete mixing ddbar and uubar.
  } else if (iq1 == 0) {
   if (iq2 < 3 || id == 331) {
     iq4 = (rndmPtr->flat() < 0.5) ? 1 : 2;
     // eta and eta' can also be s sbar.
     if (id == 221 && rndmPtr->flat() < fracEtass) iq4 = 3;
     if (id == 331 && rndmPtr->flat() < fracEtaPss) iq4 = 3;
     return make_pair( iq4, -iq4);
   }

  // Octet baryons.
  } else if (idAbs%10 == 2) {
    // Three identical quarks: emergency in case of higher spin 1/2 multiplet.
    if (iq1 == iq2 && iq2 == iq3) {iq4 = iq1; iq5 = 1100 * iq1 + 3;}
    // Two identical quarks, like normal p or n.
    else if (iq1 == iq2 || iq2 == iq3) {
      double rr6 = 6. * rndmPtr->flat();
      if    (iq1 == iq2 && rr6 < 2.) { iq4 = iq3; iq5 = 1100 * iq1 + 3;}
      else if             (rr6 < 2.) { iq4 = iq1; iq5 = 1100 * iq3 + 3;}
      else if (rr6 < 3.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;}
      else               { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 1;}
    // Three nonidentical quarks, Sigma- or Lambda-like.
    } else {
      int isp = (iq2 > iq3) ? 3 : 1;
      if (iq3 > iq2) swap( iq2, iq3);
      double rr12 = 12. * rndmPtr->flat();
      if      (rr12 < 4.) { iq4 = iq1; iq5 = 1000 * iq2 + 100 * iq3 + isp;}
      else if (rr12 < 5.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + isp;}
      else if (rr12 < 6.) { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + isp;}
      else if (rr12 < 9.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 4 - isp;}
      else                { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + 4 - isp;}
    }
    if (id > 0) return make_pair(  iq4,  iq5);
    else        return make_pair( -iq5, -iq4);

  // Decuplet baryons.
  } else {
    double rr3 = 3. * rndmPtr->flat();
    if (rr3 < 1.)      { iq4 = iq1; iq5 = 1000 * iq2 + 100 * iq3 + 3;}
    else if (rr3 < 2.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;}
    else               { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + 3;}
    if (id > 0) return make_pair(  iq4,  iq5);
    else        return make_pair( -iq5, -iq4);
  }

  // Done. (Fake call to avoid unwarranted compiler warning.)
  return make_pair( 0, 0);

}

//-------------------------------------------------------------------------

// Find relative momentum of colour and anticolour constituents in hadron.

double LowEnergyProcess::splitZ(int iq1, int iq2, double mRat1, double mRat2) {

  // Initial values.
  int iq1Abs = abs(iq1);
  int iq2Abs = abs(iq2);
  if (iq2Abs > 10) swap( mRat1, mRat2);
  double x1, x2, x1a, x1b;

  // Handle mesons.
  if (iq1Abs < 10 && iq2Abs < 10) {
    do x1 = pow2( mRat1 + (1. - mRat1) * rndmPtr->flat() );
    while ( pow(1. - x1, xPowMes) < rndmPtr->flat() );
    do x2 = pow2( mRat2 + (1. - mRat2) * rndmPtr->flat() );
    while ( pow(1. - x2, xPowMes) < rndmPtr->flat() );

  // Handle baryons.
  } else {
    double mRat1ab = 0.5 * mRat1 / xDiqEnhance;
    do x1a = pow2( mRat1ab + (1. - mRat1ab) * rndmPtr->flat() );
    while ( pow(1. - x1a, xPowBar) < rndmPtr->flat() );
    do x1b = pow2( mRat1ab + (1. - mRat1ab) * rndmPtr->flat() );
    while ( pow(1. - x1b, xPowBar) < rndmPtr->flat() );
    x1 = xDiqEnhance * ( x1a + x1b);
    do x2 = pow2( mRat2 + (1. - mRat2) * rndmPtr->flat() );
    while ( pow(1. - x2, xPowBar) < rndmPtr->flat() );
    if (iq2Abs > 10) swap( x1, x2);
  }

  // Return z value.
  return x1 / (x1 + x2);

}

//-------------------------------------------------------------------------

// Overestimate mass of lightest 2-body state for given flavour combination.
// Only account for one c or b in hadron, ond do not consider diquark spin.

double LowEnergyProcess::mThreshold( int iq1, int iq2) {

  // Initial values.
  int iq1Abs = abs(iq1);
  int iq2Abs = abs(iq2);
  if (iq2Abs > 10) swap( iq1Abs, iq2Abs);
  double mThr = 0.14;

  // Mesonic state.
  if (iq1Abs < 10) {
    if (iq2Abs > iq1Abs) swap( iq1Abs, iq2Abs);
    if      (iq1Abs < 3)  mThr += 0.14;
    else if (iq1Abs == 3) mThr += (iq2Abs < 3) ? 0.50 : 1.00;
    else if (iq1Abs == 4) mThr += (iq2Abs < 3) ? 1.90 : 2.00;
    else if (iq1Abs == 5) mThr += (iq2Abs < 3) ? 5.30 : 5.38;

  // Baryonic state.
  } else if (iq2Abs < 10) {
    int iqo1 = (iq1Abs/1000)%10;
    int iqo2 = (iq1Abs/100)%10;
    int iqo3 = iq2Abs;
    if (iqo3 > iqo2) swap( iqo2, iqo3);
    if (iqo2 > iqo1) swap( iqo1, iqo2);
    if      (iqo1 <  3) mThr += 0.95;
    else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.12 : ((iqo3 < 3) ? 1.33 : 1.68);
    else if (iqo1 == 4) mThr += (iqo2 < 3) ? 2.30 : ((iqo3 < 3) ? 2.48 : 2.70);
    else if (iqo1 == 5) mThr += (iqo2 < 3) ? 5.62 : ((iqo3 < 3) ? 5.80 : 6.08);

  // Baryon-antibaryon state.
  } else {
    int iqo1 = (iq1Abs/1000)%10;
    int iqo2 = (iq1Abs/100)%10;
    if      (iqo1 <  3) mThr += 0.95;
    else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.12 : 1.33;
    else if (iqo1 == 4) mThr += (iqo2 < 3) ? 2.30 : 2.48;
    else if (iqo1 == 5) mThr += (iqo2 < 3) ? 5.62 : 5.80;
    iqo1 = (iq2Abs/1000)%10;
    iqo2 = (iq2Abs/100)%10;
    if      (iqo1 <  3) mThr += 0.95;
    else if (iqo1 == 3) mThr += (iqo2 < 3) ? 1.12 : 1.33;
    else if (iqo1 == 4) mThr += (iqo2 < 3) ? 2.30 : 2.48;
    else if (iqo1 == 5) mThr += (iqo2 < 3) ? 5.62 : 5.80;
  }

  // Done.
  return mThr;

}

//-------------------------------------------------------------------------

// Pick slope b of exp(b * t) for elastic and diffractive events.

double LowEnergyProcess::bSlope( int type) {

  // Steeper slope for baryons than mesons.
  // To do: charm and bottom should have smaller slopes.
  double bA = (isBaryon1) ? 2.3 : 1.4;
  double bB = (isBaryon2) ? 2.3 : 1.4;

  // Elastic slope.
  if (type == 2)
    return 2. * bA + 2. * bB + 2. * ALPHAPRIME * log(ALPHAPRIME * sCM);

  // Single diffractive slope for XB and AX, respectively.
  if (type == 3) return 2. * bB + 2. * ALPHAPRIME * log(sCM / (mA * mA));
  if (type == 4) return 2. * bA + 2. * ALPHAPRIME * log(sCM / (mB * mB));

  // Double diffractive slope.
  return 2. * ALPHAPRIME * log(exp(4.) + sCM / (ALPHAPRIME * pow2(mA * mB)) );

}

//==========================================================================

} // end namespace Pythia8
