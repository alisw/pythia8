// main151.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian T. Preuss <christian.preuss@uni-goettingen.de>

// Keywords: MC@NLO; aMC@NLO; MadGraph5_aMC@NLO

// This is a simple program to demonstrate MC@NLO matching
// with events MadGraph5_aMC@NLO.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;
  pythia.readFile("main151.cmnd");
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Number of events.
  int nEvent = pythia.mode("Main:numberOfEvents");

  // Histogram.
  Hist pTZ("pT Z", 50, 0., 200.);

  // Event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find the final copy of the W and its pT.
    int iZ = 0;
    for (int i = pythia.event.size() - 1; i > 0; --i)
      if (pythia.event[i].idAbs() == 23) {iZ = i; break;}
    pTZ.fill(pythia.event[iZ].pT());
  }

  // Print statistics and histogram.
  pythia.stat();
  cout << pTZ;

  // Done.
  return 0;
}
