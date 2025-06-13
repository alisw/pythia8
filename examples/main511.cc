// main511.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: Hidden Valley; hadronization; shower;

// Comparison of a two-flavour Standard Model with a rescaled Hidden Valley.
// Should give almost perfect agreement, modulo minor issues like the eta'.
// Actually, the three HV setups ideally should give identical results.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Set up factor by which the HV masses and energies are rescaled
  // relative to the SM. Operates in a two-flavour scenario for simplicity.
  double rescale = 5.0;

  // Allow simple parton shower or not.
  bool doShower = false;

  // Set number of events and reference CM energy.
  int nEvent    = 1000000;
  double eCMref = 1000.;

  // Book histograms.
  int nPrimMax     = (doShower) ? 100 : 50;
  double pTPrimMax = (doShower) ? 10. : 2.;
  Hist nPrim[4], dndyPrim[4], pTPrim[4], mPrim[4];
  for (int iCase = 0; iCase < 4; ++iCase) {
    nPrim[iCase].book("multiplicity primaries", nPrimMax, -0.5, nPrimMax-0.5);
    dndyPrim[iCase].book("dn/dy primaries",     100, -10., 10.);
    pTPrim[iCase].book("pT primaries",          100, 0., pTPrimMax);
    mPrim[iCase].book("mass primaries",         100, 0., 2.);
  }

  // Loop over four cases: SM and three ways of setting HV parameters.
  for (int iCase = 0; iCase < 4; ++iCase) {
    double rescaleNow = (iCase == 0) ? 1. : rescale;

    // Generator; shorthand for event and particleData.
    Pythia pythia;
    Event& event      = pythia.event;
    ParticleData& pdt = pythia.particleData;
    double temp;

    // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
    pythia.readString("ProcessLevel:all = off");

    // Reduce printout of number of events generated.
    pythia.readString("Next:numberCount = 1000000");

    // Switch off ordinary decays.
    pythia.readString("HadronLevel:Decay = off");

    // Standard model with two light flavours; no strangeness or baryons.
    if (iCase == 0) {
      pythia.readString("StringFlav:probStoUD = 0.");
      pythia.readString("StringFlav:probQQtoQ = 0.");

    // Two-flavour Hidden Valley, like the Standard Model with only u and d.
    } else {
      pythia.readString("HiddenValley:Ngauge = 3");
      pythia.readString("HiddenValley:nFlav = 2");
      pythia.readString("HiddenValley:fragment = on");
      pythia.readString("HiddenValley:separateFlav = on");
      pythia.readString("HiddenValley:probVector = 0.34");
      // Effective parameter, weighing in also etaPrime suppression in SM.
      pythia.readString("HiddenValley:probKeepEta1 = 0.4");

      // Set up Hidden Valley masses and widths, rescaled from Standard Model.
      int idOld[8] = { 1, 2, 111, 211, 221, 113, 213, 223};
      int idNew[8] = { 4900101, 4900102, 4900111, 4900211, 4900221,
        4900113, 4900213, 4900223};
      for (int idLoop = 0; idLoop < 8; ++idLoop) {
        pdt.m0(      idNew[idLoop], rescale * pdt.m0(idOld[idLoop]) );
        if (idLoop > 1) {
          pdt.mWidth(idNew[idLoop], rescale * pdt.mWidth(idOld[idLoop]) );
          pdt.mMin(  idNew[idLoop], rescale * pdt.mMin(idOld[idLoop]) );
          pdt.mMax(  idNew[idLoop], rescale * pdt.mMax(idOld[idLoop]) );
        }
      }
    }

    // Alternatives for HV setup.
    if (iCase == 2) {
       pythia.readString("HiddenValley:setabsigma = 1");
       pythia.settings.parm("HiddenValley:LambdaNPQCD", 0.4 * pdt.m0(113));
       pythia.settings.parm("HiddenValley:LambdaNPHV",  0.4 * pdt.m0(4900113));
    } else if (iCase == 3) {
       pythia.readString("HiddenValley:setabsigma = 2");
       temp = pythia.settings.parm("StringZ:bLund");
       pythia.settings.parm("HiddenValley:bLund", temp / pow2(rescale));
       temp = pythia.settings.parm("StringPT:sigma");
       pythia.settings.parm("HiddenValley:sigmaLund", rescale * temp);
    }

    // Set up SM shower with fixed alphaS for simplicity.
    if (doShower && iCase == 0) {
      pythia.readString("TimeShower:alphaSorder = 0");
      pythia.readString("TimeShower:alphaSvalue = 0.12");
      pythia.readString("TimeShower:pTmin = 0.5");
      pythia.readString("TimeShower:QEDshowerByQ = off");
      pythia.readString("TimeShower:nGluonToQuark = 2");

    // Set up HV shower with fixed alphaS for simplicity.
    } else if (doShower) {
      pythia.readString("HiddenValley:FSR = on");
      pythia.readString("HiddenValley:alphaOrder = 0");
      pythia.readString("HiddenValley:alphaFSR = 0.12");
      pythia.settings.parm("HiddenValley:pTminFSR", rescale * 0.5);
    }

    // Exit with error if Pythia fails to initialize.
    if (!pythia.init()) return 1;

    // Initial flavours and colours.
    int idEnd  = (iCase == 0) ? 2 : 4900102;
    int iCol   = (iCase == 0) ? 101 : 0;
    // Warning: forceTimeShower does not catch already used HV colours,
    // so setting iColHV > 100 would mess up colour assignments in shower.
    int iColHV = (iCase == 0) ? 0 : 100;

    // Initial kinematics.
    double eCM  = rescaleNow * eCMref;
    double mEnd = pdt.m0(idEnd);
    double eEnd = 0.5 * eCM;
    double pEnd = sqrt( eEnd*eEnd - mEnd*mEnd);

    // Begin of event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Set up parton-level configuration.
      event.reset();
      int i1 = event.append(  idEnd, 23, iCol, 0, 0., 0.,  pEnd, eEnd, mEnd);
      int i2 = event.append( -idEnd, 23, 0, iCol, 0., 0., -pEnd, eEnd, mEnd);
      event[i1].colsHV( iColHV, 0);
      event[i2].colsHV( 0, iColHV);

      // Perform shower evolution.
      if (doShower) {
        event[1].scale( eCM);
        event[2].scale( eCM);
        pythia.forceTimeShower( 1, 2, eCM);
      }

      // Generate events. Skip if failure.
      if (!pythia.next()) continue;

      // Loop over all primary hadrons and count them.
      int nPrimary = 0;
      for (int i = 0; i < event.size(); ++i) {
        if (event[i].status() < 0) continue;
        ++nPrimary;

        // Rapidity, pT and mass spectra of primary hadrons.
        dndyPrim[iCase].fill( event[i].y() );
        pTPrim[iCase].fill( event[i].pT() / rescaleNow);
        mPrim[iCase].fill( event[i].m() / rescaleNow);

      // End of hadron loop. Number of primary hadrons.
      }
      nPrim[iCase].fill( nPrimary );

    // End of event loop. Print statistics and rescale histograms.
    }
    pythia.stat();
    nPrim[iCase]    *=   1. / double(nEvent);
    dndyPrim[iCase] *=   5. / double(nEvent);
    pTPrim[iCase]   *= 100. / (pTPrimMax * double(nEvent));
    mPrim[iCase]    *=  50. / double(nEvent);
    cout << nPrim[iCase] << dndyPrim[iCase] << pTPrim[iCase] << mPrim[iCase];

  // End of case loop.
  }

  // Write Python code for plotting.
  HistPlot hpl("plot511");
  hpl.frame("fig511", "Primary hadron multiplicity distribution",
    "$n$", "$\\mathrm{d}P/\\mathrm{d}n$");
  hpl.add( nPrim[0], "h,red", "SM");
  hpl.add( nPrim[1], "h,black", "HV set = 0");
  hpl.add( nPrim[2], "h,blue", "HV set = 1");
  hpl.add( nPrim[3], "h,green", "HV set = 2");
  hpl.plot();
  hpl.frame("", "Primary hadron rapidity distribution",
    "$y$", "$\\mathrm{d}N/\\mathrm{d}y$");
  hpl.add( dndyPrim[0], "-,red", "SM");
  hpl.add( dndyPrim[1], "-,black", "HV set = 0");
  hpl.add( dndyPrim[2], "--,blue", "HV set = 1");
  hpl.add( dndyPrim[3], ".,green", "HV set = 2");
  hpl.plot();
  hpl.frame("", "Primary hadron transverse momentum distribution (rescaled)",
    "$p_{\\perp}$", "$\\mathrm{d}N/\\mathrm{d}p_{\\perp}$");
  hpl.add( pTPrim[0], "-,red", "SM");
  hpl.add( pTPrim[1], "-,black", "HV set = 0");
  hpl.add( pTPrim[2], "--,blue", "HV set = 1");
  hpl.add( pTPrim[3], ".,green", "HV set = 2");
  hpl.plot();
  hpl.frame("", "Primary hadron mass distribution (rescaled)",
    "$m$", "$\\mathrm{d}N/\\mathrm{d}m$");
  hpl.add( mPrim[0], "h,red", "SM");
  hpl.add( mPrim[1], "h,black", "HV set = 0");
  hpl.add( mPrim[2], "h,blue", "HV set = 1");
  hpl.add( mPrim[3], "h,green", "HV set = 2");
  hpl.plot();

  // Done.
  return 0;
}
