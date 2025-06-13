// main216.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; jet finding;

// This is a simple test program, convenient for tutorials.
// It illustrates searches for a 1 TeV Z' decaying to two jets.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events, generated and listed ones.
  int nEvent    = 1000;
  int nListJets = 5;

  // Pythia generator at LHC energy; pp is default beam setup.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");

  // Z' production with 1 TeV mass (+- 200 GeV) and decay only to d/u/s/c/b.
  pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
  pythia.readString("32:m0 = 1000.");
  pythia.readString("PhaseSpace:mHatMin = 800.");
  pythia.readString("PhaseSpace:mHatMax = 1200.");
  pythia.readString("32:onMode = off");
  pythia.readString("32:onIfAny = 1 2 3 4 5");

  // Put other settings here.

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;


  // Required parameters for an anti-kT jet finder.
  int    jetPower  = -1;
  double jetRadius = 0.7;
  double jetpTmin  = 100.;
  // Exclude high-eta and neutrinos (and other invisible) from study.
  double etaMax    = 4.;
  int    nSel      = 2;
  // Set up jet finder, with pion mass assumed for non-photons.
  SlowJet slowJet( jetPower, jetRadius, jetpTmin, etaMax, nSel, 1);

  // Histograms.
  Hist massZp("true Zprime mass", 100, 0., 2000.);
  Hist mass2j("two-jet mass", 100, 0., 2000.);
  Hist deltaM("mass shift", 100, -500., 500.);
  Hist deltaMz("mass shift zoomed", 100, -50., 50.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // A copy of the Z' is always stored in slot 5 of the event record.
    double mZp = pythia.event[5].m();
    massZp.fill( mZp);

    // Analyze Slowet jet properties. List first few.
    slowJet.analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Form invariant mass of two highest-pT jets.
    if ( slowJet.sizeJet() > 1) {
      double m2j = (slowJet.p(0) + slowJet.p(1)).mCalc();
      mass2j.fill( m2j);

      // Shift in mass from true to reconstructed.
      deltaM.fill( m2j - mZp);
      deltaMz.fill( m2j - mZp);
    }

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Normalize histograms to rate per GeV and print them.
  massZp *= 0.05 / double(nEvent);
  mass2j *= 0.05 / double(nEvent);
  deltaM *= 0.10 / double(nEvent);
  deltaMz *= 1.0 / double(nEvent);
  cout << massZp << mass2j << deltaM << deltaMz;

  // Plot histograms for nicer output.
  HistPlot hpl("plot216");
  hpl.frame( "fig216", "Zprime mass distribution", "$m$ (GeV)",
    "Probability");
  hpl.add( massZp, "h,blue", "true mass");
  hpl.add( mass2j, "h,red", "reconstructed mass");
  hpl.plot();
  hpl.plotFrame("", deltaM, "Zprime mass shift from true to reconstructed",
    "$\\Delta m$ (GeV)", "Probability", "h,magenta", "reconstructed$-$true");
  hpl.plotFrame("", deltaMz,
    "Zprime mass shift from true to reconstructed (zoomed-in)",
    "$\\Delta m$ (GeV)", "Probability", "h,magenta", "reconstructed$-$true");

  // Done.
  return 0;
}
