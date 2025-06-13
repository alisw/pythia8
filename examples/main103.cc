// main103.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: electron-positron; basic usage; charged multiplicity

// This is a simple test program. It fits on one slide in a talk.
// It studies the charged multiplicity distribution at LEP 1.

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Incoming e+e- beams ar LEP1 energies. No photons from e+-.
  Pythia pythia;
  pythia.readString("Beams:idA = -11");
  pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:eCM = 91.2");
  pythia.readString("PDF:lepton = off");
  // Process: hadronic decays of Z0. Initialize. Histogram.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;
  Hist mult("charged multiplicity", 100, -0.5, 99.5);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << mult;
  return 0;
}
