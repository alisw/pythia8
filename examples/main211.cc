// main211.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; jet finding; slowjet;

// This is a simple test program.
// It studies jet production at the LHC, using SlowJet,
// here as a Pythia-centric wrapper for FastJet.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events, generated and listed ones.
  int nEvent    = 1000;
  int nListJets = 5;

  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 200.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Standard parameters for a jet finder.
  double radius   = 0.7;
  double pTjetMin = 10.;
  double etaMax   = 4.;
  // Exclude neutrinos (and other invisible) from study.
  int    nSel     = 2;

  // Set up SlowJet jet finder, with anti-kT clustering (-1)
  // and pion mass assumed for non-photons.
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

  // Histograms.
  Hist nJets("number of jets", 50, -0.5, 49.5);
  Hist pTjets("pT for jets", 100, 0., 500.);
  Hist yJets("y for jets", 100, -5., 5.);
  Hist phiJets("phi for jets", 100, -M_PI, M_PI);
  Hist distJets("R distance between jets", 100, 0., 10.);
  Hist pTdiff("pT difference between jets", 100, -100., 400.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Analyze Slowet jet properties. List first few.
    slowJet.analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Fill SlowJet inclusive jet distributions.
    nJets.fill( slowJet.sizeJet() );
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      pTjets.fill( slowJet.pT(i) );
      yJets.fill( slowJet.y(i) );
      phiJets.fill( slowJet.phi(i) );
    }

    // Fill SlowJet distance between jets.
    for (int i = 0; i < slowJet.sizeJet() - 1; ++i)
    for (int j = i + 1; j < slowJet.sizeJet(); ++j) {
      double dY = slowJet.y(i) - slowJet.y(j);
      double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dY) + pow2(dPhi) );
      distJets.fill( dR );
    }

    // Fill SlowJet pT-difference between jets (to check ordering of list).
    for (int i = 1; i < slowJet.sizeJet(); ++i)
      pTdiff.fill( slowJet.pT(i-1) - slowJet.pT(i) );

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();
  cout << nJets << pTjets << yJets << phiJets << distJets << pTdiff;

  // Done.
  return 0;
}
