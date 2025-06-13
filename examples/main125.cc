// main125.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; LHE file

// This is a simple test program. It shows how PYTHIA 8 can write
// a Les Houches Event File v. 3.0 with hadron-level events, based on
// its process-level input events. Event weights are carried along.
// Note: HepMC is more commonly used for hadron-level output than LHEF.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Set up for external input of LHEF 3.0 process-level events.
  Pythia pythia;
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = wbj_lhef3.lhe");
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Create and open file for LHEF 3.0 hadron-level output.
  LHEF3FromPythia8 myLHEF3(&pythia.event, &pythia.info);
  myLHEF3.openLHEF("weakbosons.lhe");

  // Write out initialization info on the file.
  myLHEF3.setInit();

  // Event generation loop.
  for (int iEvent = 0; iEvent < 10; ++iEvent) {

    // Generate next event.
    if (!pythia.next()) {
      if( pythia.info.atEndOfFile() ) break;
      else continue;
    }

    // Store and write event info.
    myLHEF3.setEvent();

  } // End loop over events to generate.

  // Statistics: full printout.
  pythia.stat();

  // Write endtag. Overwrite initialization info with new cross sections.
  myLHEF3.closeLHEF(true);

  // Done.
  return 0;
}
