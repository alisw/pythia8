// main510.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Lisa Carloni
//          Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: Hidden Valley;

// Hidden Valley: omnibus version for various scenarios in 1 TeV e+e-.
// See L. Carloni, J. Rathsman and T. Sjostrand, JHEP 04 (2011) 091
// for details.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

// Avoid division by zero when calculating cosTheta.
const double PABSMIN = 1e-10;

//==========================================================================

int main() {

  // Number of events and CM energy.
  int nEvent = 10000;
  double eCM = 1000.;

  // Superselection of scenarios: 0 = free, 1 = abelian, 2 = non-abelian.
  int supersel = 0;

  // Selection of scenarios. Based on classification in JHEP 04 (2011) 091.
  // 1 = KMA;  e+ e- -> Z' -> qv qvbar -> gammav.
  // 2 = KMNA; e+ e- -> Z' -> qv qvbar -> piv, rhov.
  // 3 = SMA;  e+ e- -> gamma* -> Ev Evbar -> e+ e- qv qvbar -> gammav.
  // 4 = SMNA; e+ e- -> gamma* -> Ev Evbar -> e+ e- qv qvbar -> piv, rhov.
  // 5 = KMA;  e+ e- -> gamma* -> qv qvbar -> gammav.
  // 6 = KMNA; e+ e- -> gamma* -> qv qvbar -> piv, rhov.
  int scenario   = 4;

  // Key mass choices: Ev, qv, gammav = piv.
  // Note: if NA then mqv = mgampiv/2 is forced, i.e. mqv ignored.
  double mEv     = 250.;
  double mqv     = 50.;
  double mgampiv = 10.0;

  // Number of qv replicas in NA models. Production rhov/(piv + rhov).
  int nFlav      = 4;
  double probVec = 0.75;

  // Coupling strength alpha in showers. (Here frozen; can be made running.)
  int alphaHVorder = 0;
  double alphaHV = 0.2;

  // a and b (scaled) parameter in fragmentation function.
  double aLund     = 0.68;    // default 0.68; larger means more particles
  double bLund     = 0.98;    // default 0.98: larger means fewer particles
  double sigmaLund = 0.335;   // default 0.335; larger means broader pT

  // Allow or not ISR of photons off incoming e+e-
  bool allowISR  = true;

  // Max angle for beam pipe. Used to remove much ISR.
  double theBeam = 0.050;

  // Isolation cone for leptons: require no charged inside it,
  // and at most given fraction neutral energy.
  double theIsol  = 0.10;
  double fracIsol = 0.20;

  //------------------------------------------------------------------------

  // Superselection overwrites above values.

  // Abelian (A) scenario.
  if (supersel == 1) {
    scenario = 1;
    mqv      = 20.;
    mgampiv  = 10.;
    alphaHV  = 0.3;
  }

  // Non-Abelian (NA) scenario.
  if (supersel == 2) {
    scenario = 2;
    mqv      = 5.;
    mgampiv  = 10.;
    alphaHV  = 0.15;
    nFlav    = 4;
    // Example of (unrealistic?) changed fragmentation parameters.
    aLund     = 0.12;
    bLund     = 2.0;
    sigmaLund = 0.35;
    // Rescale b and sigma, based on ratio m(rho_HV) / m(rho_QCD).
    double rescale = mgampiv / 0.775;
    bLund         /= pow2(rescale);
    sigmaLund     *= rescale;
  }

  //------------------------------------------------------------------------

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up specifics of each scenario.

  // KMA scenario.
  if (scenario == 1) {
    pythia.readString("HiddenValley:Ngauge = 1");
    pythia.readString("HiddenValley:spinqv = 0");

    // f fbar -> Zv -> qv qvbar.
    pythia.readString("HiddenValley:ffbar2Zv = on");
    pythia.readString("4900023:onMode = off");
    pythia.readString("4900023:onIfAny = 4900101");

  // KMNA scenario.
  } else if (scenario == 2) {
    pythia.readString("HiddenValley:Ngauge = 3");
    pythia.readString("HiddenValley:spinFv = 2");

    // f fbar -> Zv -> qv qvbar.
    pythia.readString("HiddenValley:ffbar2Zv = on");
    pythia.readString("4900023:onMode = off");
    pythia.readString("4900023:onIfAny = 4900101");

  // SMA scenario.
  } else if (scenario == 3) {
    pythia.readString("HiddenValley:Ngauge = 1");
    pythia.readString("HiddenValley:spinFv = 1");
    pythia.readString("HiddenValley:spinqv = 0");
    pythia.readString("HiddenValley:kappa  = 1.");

    // f fbar -> Fv Fvbar -> f qv fbar qvbar
    pythia.readString("HiddenValley:ffbar2EvEvbar = on");

  // SMNA scenario.
  } else if (scenario == 4) {
    pythia.readString("HiddenValley:Ngauge = 3");
    pythia.readString("HiddenValley:spinFv = 1");
    pythia.readString("HiddenValley:spinqv = 0");
    pythia.readString("HiddenValley:kappa  = 1.");

    // f fbar -> Fv Fvbar -> f qv fbar qvbar
    pythia.readString("HiddenValley:ffbar2EvEvbar = on");

  // KMA scenario, with Ev masking as qv.
  } else if (scenario == 5) {
    pythia.readString("HiddenValley:Ngauge = 1");
    pythia.readString("HiddenValley:spinFv = 1");

    // f fbar -> gamma/Z -> qv qvbar.
    pythia.readString("HiddenValley:doKinMix = on");
    pythia.readString("HiddenValley:ffbar2EvEvbar = on");

  // KMNA scenario, with Ev masking as qv.
  } else if (scenario == 6) {
    pythia.readString("HiddenValley:Ngauge = 3");
    pythia.readString("HiddenValley:spinFv = 1");

    // f fbar -> gamma/Z -> qv qvbar.
    pythia.readString("HiddenValley:doKinMix = on");
    pythia.readString("HiddenValley:ffbar2EvEvbar = on");
  }

  // For NA allow HV-fragmentation, with/without off-diagonal mesons.
  if (scenario%2 == 0) {
    pythia.readString("HiddenValley:fragment = on");
    pythia.settings.mode("HiddenValley:nFlav", nFlav);
    pythia.settings.parm("HiddenValley:probVector", probVec);
    if (supersel == 2) {
      pythia.readString("HiddenValley:setabsigma = 2");
      pythia.settings.parm("HiddenValley:aLund", aLund);
      pythia.settings.parm("HiddenValley:bLund", bLund);
      pythia.settings.parm("HiddenValley:sigmaLund", sigmaLund);
    }
  }

  //------------------------------------------------------------------------

  // Set up common (or not used) parameters. Initialize generator.

  // In non-Abelian scenarios force m_qv = m_piv/2 = m_rhov/2.
  if (scenario%2 == 0) mqv = 0.5 * mgampiv;

  // Set Ev, qv, gamma_v, piv, rhov masses.
  pythia.particleData.m0( 4900011, mEv );
  pythia.particleData.m0( 4900101, mqv );
  pythia.particleData.m0( 4900022, mgampiv );
  pythia.particleData.m0( 4900111, mgampiv );
  pythia.particleData.m0( 4900113, mgampiv );
  pythia.particleData.m0( 4900211, mgampiv );
  pythia.particleData.m0( 4900213, mgampiv );

  // In scenarios 5 and 6 Ev masks for qv.
  if (scenario == 5 || scenario == 6) {
    pythia.particleData.m0( 4900011, mqv );
    pythia.particleData.mWidth( 4900011, 0. );
    pythia.readString("4900011:mayDecay = false");
    pythia.readString("4900011:chargeType = 0");
    pythia.readString("4900011:isVisible = false");
  }

  // In scenario 3 and 4 set mEv min
  if (scenario == 3 || scenario == 4) {
    pythia.particleData.mMin( 4900011, 0.1*mEv );
  }

  // Set gamma_v/pi_v/rho_v decay on, except channels close to threshold.
  pythia.readString("4900022:mayDecay = true");
  pythia.readString("4900022:onMode = on");
  pythia.readString("4900111:onMode = on");
  pythia.readString("4900113:onMode = on");
  if (mgampiv < 11.0) {
    pythia.readString("4900022:offIfAny = 5");
    pythia.readString("4900111:offIfAny = 5");
    pythia.readString("4900113:offIfAny = 5");
  }
  if (mgampiv < 4.0) {
    pythia.readString("4900022:offIfAny = 4");
    pythia.readString("4900111:offIfAny = 4");
    pythia.readString("4900113:offIfAny = 4");
  }
  if (mgampiv < 3.6) {
    pythia.readString("4900022:offIfAny = 15");
    pythia.readString("4900111:offIfAny = 15");
    pythia.readString("4900113:offIfAny = 15");
  }

  // Valley shower parameters.
  pythia.readString("HiddenValley:FSR = on");
  pythia.settings.mode("HiddenValley:alphaOrder", alphaHVorder);
  pythia.settings.parm("HiddenValley:alphaFSR", alphaHV);
  pythia.settings.parm("HiddenValley:pTminFSR", max(0.5, 0.1 * mgampiv));

  // Switch off the photons from ISR (from the PDF).
  pythia.settings.flag("PDF:lepton", allowISR);

  // Initialization for e+e- at selected CM energy.
  pythia.readString("Beams:idA = -11");
  pythia.readString("Beams:idB =  11");
  pythia.settings.parm("Beams:eCM", eCM);

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  //------------------------------------------------------------------------

  // Histograms. Initialize event analysis.

  // Basic plots.
  Hist nBackH("number of decays back into SM", 50, -0.5, 49.5);
  Hist nNonDH("number of non-diagonal mesons", 50, -0.5, 49.5);
  Hist eBackH("total energy radiated back", 100, 0., eCM);
  Hist eSpBackH("energy spectrum of decays back", 100, 0., 0.5 * eCM);
  Hist theBackH("theta angle between decays back", 100, 0., M_PI);
  Hist theBackFH("theta angle between decays back, fast", 100, 0., M_PI);
  Hist theBackSH("theta angle between decays back, slow", 100, 0., M_PI);
  Hist pTmissH("missing pT", 100, 0., 0.5 * eCM);
  Hist eISRH("total energy in ISR photons", 100, 0., eCM);
  Hist mISRH("reduced mass after ISR photons", 100, 0., eCM);
  Hist nChargedH("charged multiplicity in detector", 100, -0.5, 99.5);

  // Lepton spectra.
  Hist eLeptonH("lepton energy spectum", 100, 0., 0.5 * eCM);
  Hist isolLeptonH("isolated lepton fraction (E)", 100, 0., 0.5 * eCM);
  Hist neLeptonH("neutral energy fraction to lepton", 100, 0., 1.0);
  Hist eelimLeptonH("eliminated lepton energy spectum", 100, 0., 0.5 * eCM);
  Hist ekeptLeptonH("kept lepton energy spectum", 100, 0., 0.5 * eCM);
  Hist mLeptonPairH("kept lepton pair mass spectrum", 100, 0., 2.5 * mgampiv);

  // Hadron and lepton spectra.
  Hist mHadH("invariant mass of non-leptons", 100, 0.,  2.5 * mgampiv);
  Hist mAllLowH("invariant mass of all, low range", 100, 0.,  2.5 * mgampiv);
  Hist mAllHighH("invariant mass of all, full range", 100, 0., eCM);

  // Jets with Jade algorithm.
  Hist nJetH("number of jets", 100, -0.5, 99.5);
  Hist mJetH("mass distribution of jets", 100, 0., 5. * mgampiv);
  Hist cosTheBeamAllH("cosTheta to beam, all jets", 100, 0., 1.);
  Hist cosTheBeamMassH("cosTheta to beam, massive jets", 100, 0., 1.);
  Hist thetaPairAllH("theta of jet pair, all jets", 100, 0., M_PI);
  Hist thetaPairMassH("theta of jet pair, massive jets", 100, 0., M_PI);
  Hist cosThePairAllH("cosTheta of jet pair, all jets", 100, -1., 1.);
  Hist cosThePairMassH("cosTheta of jet pair, massive jets", 100, -1., 1.);

  // Linearized sphericity and thrust.
  Hist spheriAllH("sphericity almost all events", 100, 0., 1.);
  Hist spheriMassH("sphericity massive events", 100, 0., 1.);
  Hist thrustAllH("thrust almost all events", 100, 0.5, 1.);
  Hist thrustMassH("thrust massive events", 100, 0.5, 1.);

  // Shorthand for some HV codes.
  int gammav = 4900022;
  int piv    = 4900111;
  int rhov   = 4900113;
  int pivND  = 4900211;
  int rhovND = 4900213;

  // Derived quantities for analysis.
  double cosTheBeam = cos(theBeam);
  double cosTheIsol = cos(theIsol);

  // Set up linearized Sphericity and Thrust.
  Sphericity sph(1.);
  Thrust thr;

  // Set up Jade jet finding and use m_{gamma_v/pi_v/rho_v} as scale.
  // Note: easy to switch to Durham instead.
  ClusterJet jade("Jade");
  double mJade   = mgampiv;
  double mJetMin = 0.5 * mgampiv;

  //------------------------------------------------------------------------

  // Begin event generation loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events.
    if (!pythia.next()) continue;

    // First analysis of event, searching for gammav/piv/rhov.
    int    nBack = 0;
    double eBack = 0.;
    vector<int> iBack;
    int    nNonD = 0;
    Vec4 pSumVis = 0.;
    Vec4 pISR    = 0.;
    for (int i = 0; i < event.size(); ++i) {
      int id = event[i].id();
      bool isgampirho  = (id == gammav || id == piv || id == rhov);
      bool noCopyBelow = (event[i].isFinal()
        || (event[i].daughter1() != event[i].daughter2()));

      // Histogram gammav/piv/rhov. Sum energy in ISR photons.
      if (isgampirho && noCopyBelow) {
        ++nBack;
        iBack.push_back(i);
        double eNow = event[i].e();
        eBack += eNow;
        eSpBackH.fill( eNow);
      }
      if (abs(id) == pivND || abs(id) == rhovND) ++nNonD;
      if (event[i].isFinal() && event[i].isVisible())
        pSumVis += event[i].p();
      if (id == 22 && (event[i].status() == 43 || event[i].status() == 44
        || event[i].status() == 62 || event[i].status() == 63))
        pISR += event[i].p();

    // End of particle loop. Histogram properties.
    }
    nBackH.fill( nBack);
    eBackH.fill( eBack);
    nNonDH.fill( nNonD);
    double eSep = 5. * mgampiv;
    for (int j1 = 0; j1 < nBack; ++j1)
    for (int j2 = j1 + 1; j2 < nBack; ++j2) {
      Vec4 p1 = event[iBack[j1]].p();
      Vec4 p2 = event[iBack[j2]].p();
      double the12 = theta( p1, p2);
      theBackH.fill( the12 );
      if (p1.e() > eSep && p2.e() > eSep) theBackFH.fill( the12 );
      if (p1.e() < eSep && p2.e() < eSep) theBackSH.fill( the12 );
    }
    pTmissH.fill( pSumVis.pT());
    eISRH.fill( pISR.e());
    mISRH.fill( (Vec4( 0., 0., 0., eCM) - pISR).mCalc() );

    // Remove particles near beampipe.
    double cosTheNow;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
      cosTheNow = abs(event[i].pz()) / max(PABSMIN, event[i].pAbs());
      if (cosTheNow > cosTheBeam) event[i].statusNeg();
    }

    // Charged multiplicity in detector.
    int nCharged = 0;
    for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) ++nCharged;
    nChargedH.fill( nCharged);

    // Find all final leptons (e+-, mu+-) in an event.
    vector<int>  iLepton;
    vector<bool> isIsol;
    vector<bool> isKept;
    for (int i = 0; i < event.size(); ++i) if ( event[i].isFinal()
    && (event[i].idAbs() == 11 || event[i].idAbs() == 13) ) {
      iLepton.push_back(i);
      isKept.push_back(true);

      // Number of charged and neutral particles in cone around each.
      int nChIsol = 0;
      vector<int> iNeutral;
      double eLepton = event[i].e();
      double eIsol   = 0.;
      for (int j = 0; j < event.size(); ++j)
        if (j != i && event[j].isFinal() && event[i].pAbs2() != 0
          && event[j].pAbs2() != 0
          && costheta( event[i].p(), event[j].p() ) > cosTheIsol) {
        if (event[j].isCharged()) ++nChIsol;
        else {
          iNeutral.push_back(j);
          eIsol += event[j].e();
        }
      }

      // Lepton isolated or not; for isolated resum any photons.
      bool isIsolNow = (nChIsol == 0 && eIsol < fracIsol * eLepton);
      isIsol.push_back( isIsolNow);
      if (isIsolNow)
      for (int j = 0; j < int(iNeutral.size()); ++j) {
        event[i].p() += event[iNeutral[j]].p();
        event[iNeutral[j]].statusNeg();
      }

      // Histogram lepton properties.
      eLeptonH.fill( eLepton);
      if (isIsolNow) isolLeptonH.fill( eLepton);
      neLeptonH.fill( eIsol / eLepton);
    }

    // For SM scenarios search out highest-energy e+ and e-.
    if (scenario == 3 || scenario == 4) {
      int iEpos = 0;
      int iEneg = 0;
      int jEpos = 0;
      int jEneg = 0;
      double eEpos = 0.;
      double eEneg = 0.;
      for (int j = 0; j < int(iLepton.size()); ++j) {
        int i = iLepton[j];
        if (event[i].id() == -11 && event[i].e() > eEpos) {
          iEpos = i;
          jEpos = j;
          eEpos = event[i].e();
        }
        if (event[i].id() == 11 && event[i].e() > eEneg) {
          iEneg = i;
          jEneg = j;
          eEneg = event[i].e();
        }
      }

      // Eliminate them from the event record and mark in lepton list.
      if (iEpos > 0) {
        event[iEpos].statusNeg();
        isKept[jEpos] = false;
        eelimLeptonH.fill( eEpos);
      }
      if (iEneg > 0) {
        event[iEneg].statusNeg();
        isKept[jEneg] = false;
        eelimLeptonH.fill( eEneg);
      }
    }

    // Energy of remaining leptons and mass spectrum of pairs.
    for (int j1 = 0; j1 < int(iLepton.size()); ++j1) if (isKept[j1]) {
      int i1 = iLepton[j1];
      ekeptLeptonH.fill( event[i1].e() );
      for (int j2 = j1 + 1; j2 < int(iLepton.size()); ++j2) if (isKept[j2]) {
        int i2 = iLepton[j2];
        if (event[i1].id() + event[i2].id() == 0) {
          double mPair = (event[i1].p() + event[i2].p()).mCalc();
          mLeptonPairH.fill( mPair );
        }
      }
    }

    // Momentum sum of all visible and of non-leptons.
    Vec4 pSumHad, pSumAll;
    int nSumHad = 0;
    int nSumAll = 0;
    for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal() && event[i].isVisible()) {
      pSumAll += event[i].p();
      ++nSumAll;
      if (event[i].idAbs() != 11 && event[i].idAbs() != 13) {
        pSumHad += event[i].p();
        ++nSumHad;
      }
    }
    double mSumAll = pSumAll.mCalc();
    if (nSumHad > 0) mHadH.fill( pSumHad.mCalc() );
    if (nSumAll > 0) mAllLowH.fill( mSumAll );
    if (nSumAll > 0) mAllHighH.fill( mSumAll );

    // Cluster jets and study angles to beam axis and between jets.
    double mJ, mK, costheBeam, thetaPair, cosThePair;
    if ( mSumAll > mJetMin && jade.analyze(event, 0., mJade, 1, 0) ) {
      nJetH.fill( jade.size());
      for (int j = 0; j < jade.size(); ++j) {
        mJ = jade.p(j).mCalc();
        mJetH.fill( mJ);
        costheBeam = abs( jade.p(j).pz()) / max(PABSMIN, jade.p(j).pAbs());
        cosTheBeamAllH.fill( costheBeam);
        if (mJ > mJetMin) cosTheBeamMassH.fill( costheBeam);
        for (int k = 0; k < j; ++k) {
          mK = jade.p(k).mCalc();
          thetaPair = theta( jade.p(j), jade.p(k));
          cosThePair = cos(thetaPair);
          thetaPairAllH.fill( thetaPair);
          if (mJ > mJetMin && mK > mJetMin) thetaPairMassH.fill( thetaPair);
          cosThePairAllH.fill( cosThePair);
          if (mJ > mJetMin && mK > mJetMin) cosThePairMassH.fill( cosThePair);
        }
      }
    }

    // Linearized Sphericity and Thrust analysis.
    if ( mSumAll > mJetMin) {

      // Bost to rest frame of system and check it worked.
      event.bst( -pSumAll.px()/pSumAll.e(), -pSumAll.py()/pSumAll.e(),
                 -pSumAll.pz()/pSumAll.e() );
      Vec4 pSumCheck;
      for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isVisible())
        pSumCheck += event[i].p();
      if (pSumCheck.pAbs() > 0.1) cout << " Error: failed boost to rest "
        << pSumCheck;

      // Linearized sphericity.
      if (sph.analyze( event)) {
        spheriAllH.fill( sph.sphericity());
        if (mSumAll > 2. * mgampiv) spheriMassH.fill( sph.sphericity());
      }

      // Thrust.
      if (thr.analyze( event)) {
        thrustAllH.fill( thr.thrust());
        if (mSumAll > 2. * mgampiv) thrustMassH.fill( thr.thrust());
      }
    }

  // End of event loop.
  }

  // Statistics. Histograms.
  pythia.stat();
  isolLeptonH /= eLeptonH;
  cout << nNonDH << nBackH << eBackH << eSpBackH << theBackH << theBackFH
       << theBackSH << pTmissH << eISRH << mISRH << nChargedH
       << eLeptonH << isolLeptonH << neLeptonH << eelimLeptonH
       << ekeptLeptonH << mLeptonPairH << mHadH << mAllLowH
       << mAllHighH << nJetH << mJetH << cosTheBeamAllH
       << cosTheBeamMassH << thetaPairAllH << thetaPairMassH
       << cosThePairAllH << cosThePairMassH  << spheriAllH
       << spheriMassH << thrustAllH << thrustMassH;

  // Done.
  return 0;
}
