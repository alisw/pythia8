// main111.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: electron-positron; basic usage; command file

// This is a simple test program, equivalent with main103.cc,
// but with data-base settings collected in a separate .cmnd file.
// It studies the charged multiplicity distribution at LEP 1.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator. Read in command file. Initialize.
  Pythia pythia;
  pythia.readFile("main111.cmnd");
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Read number of events to generate. Book histogram.
  int nEvent = pythia.mode("Main:numberOfEvents");
  Hist mult("charged multiplicity", 100, -0.5, 99.5);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

  // End of event loop. Statistics. Normalized histogram. Done.
  }
  pythia.stat();
  mult /= nEvent;
  cout << mult;
  return 0;
}
