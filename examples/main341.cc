// main341.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; DIS

// Basic setup for Deeply Inelastic Scattering at HERA.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Beam energies, minimal Q2, number of events to generate.
  double eProton   = 920.;
  double eElectron = 27.5;
  double Q2min     = 25.;
  int    nEvent    = 10000;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up incoming beams, for frame with unequal beam energies.
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton.
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA", eProton);
  // BeamB = electron.
  pythia.readString("Beams:idB = 11");
  pythia.settings.parm("Beams:eB", eElectron);

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // Uncomment to allow charged current.
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Histograms.
  double Wmax = sqrt(4.* eProton * eElectron);
  Hist Qhist("Q [GeV]", 100, 0., 50.);
  Hist Whist("W [GeV]", 100, 0., Wmax);
  Hist xhist("x", 100, 0., 1.);
  Hist yhist("y", 100, 0., 1.);
  Hist pTehist("pT of scattered electron [GeV]", 100, 0., 50.);
  Hist pTrhist("pT of radiated parton [GeV]", 100, 0., 50.);
  Hist pTdhist("ratio pT_parton/pT_electron", 100, 0., 5.);

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Fill kinematics histograms.
    Qhist.fill( sqrt(pythia.info.Q2DIS()) );
    Whist.fill( pythia.info.WDIS() );
    xhist.fill( pythia.info.xDIS() );
    yhist.fill( pythia.info.yDIS() );
    pTehist.fill( event[6].pT() );

    // pT spectrum of partons being radiated in shower.
    for (int i = 0; i < event.size(); ++i) if (event[i].statusAbs() == 43) {
      pTrhist.fill( event[i].pT() );
      pTdhist.fill( event[i].pT() / event[6].pT() );
    }

  // End of event loop. Statistics and histograms.
  }
  pythia.stat();
  cout << Qhist << Whist << xhist << yhist << pTehist << pTrhist << pTdhist;

  // Done.
  return 0;
}
