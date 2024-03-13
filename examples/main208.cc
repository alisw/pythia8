// main208.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This example illustrates Vincia (or Pythia) showers for
// double-dissociative photon-initiated gamma-gamma -> mu+mu- at LHC.

#include "Pythia8/Pythia.h"

using namespace Pythia8;
int main() {

  // Construct Pythia.
  Pythia pythia;

  // Set up the process.
  pythia.readString("PhotonCollision:gmgm2mumu = on");
  pythia.readString("PhaseSpace:pTHatMin = 5.");

  // Use Vincia (or Pythia) showers.
  pythia.readString("PartonShowers:model = 2"); // 1: Pythia, 2: Vincia.

  // Simple shower settings.
  // pythia.readString("PartonShowers:model = 1");
  // pythia.readString("SpaceShower:dipoleRecoil = on");
  // pythia.readString("SpaceShower:pTmaxMatch = 2");
  // pythia.readString("SpaceShower:pTdampMatch = 1");

  // Vincia shower settings.
  pythia.readString("Vincia:EWMode = 2"); // vincia Multipole QED shower
  pythia.readString("Vincia:pTmaxMatch = 2"); // 2: power showers
  pythia.readString("Vincia:pTdampMatch = 1"); // 1: dampening.

  // Switch off hadronisation and MPI for simplicity and speed.
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");

  // Initialise
  pythia.init();

  // Book a histogram.
  // Second-to-last argument => use logarithmic x axis.
  // Last argument = true => compute stat uncertainties.
  Hist hisPTgam("log(pT(gamma)/GeV)", 100, 1.e-6, 10000., true, true);
  Hist hisPTele("log(pT(e)/GeV)", 100, 1.e-6, 10000., true, true);
  Hist hisPTmu("log(pT(mu)/GeV)", 100, 1.e-6, 10000., true, true);
  Hist hisPTtau("log(pT(tau)/GeV)", 100, 1.e-6, 10000., true, true);
  Hist hisPToth("log(pT(other)/GeV)", 100, 1.e-6, 10000., true, true);

  // Number of events.
  int nEvent = 1000;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate next event.
    pythia.next();

    // Look for final-state photons.
    for ( int i = 1; i < pythia.event.size(); ++i ) {
      if ( !pythia.event[i].isFinal() ) continue;
      if ( pythia.event[i].id() == 22) {
        hisPTgam.fill(pythia.event[i].pT());
      } else if ( pythia.event[i].idAbs() == 11) {
        hisPTele.fill(pythia.event[i].pT());
      } else if ( pythia.event[i].idAbs() == 13) {
        hisPTmu.fill(pythia.event[i].pT());
      } else if ( pythia.event[i].idAbs() == 15) {
        hisPTtau.fill(pythia.event[i].pT());
      } else {
        hisPToth.fill(pythia.event[i].pT());
      }
    }

  }

  // Give statistics.
  pythia.stat();

  // Stat output
  cout << hisPTgam << hisPTele << hisPTmu << hisPTtau << hisPToth;

  // Done.
  return 0;
}
