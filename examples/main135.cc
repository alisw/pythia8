// main135.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Mikhail Kirsanov <Mikhail.Kirsanov@cern.ch>
//          Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: basic usage; hepmc

// This program illustrates how HepMC files can be written by Pythia8,
// but reduced in size by only including final-state particles.
// HepMC events are output to the main135.hepmc file.

// WARNING: typically one needs 10 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "Pythia8/Pythia.h"

// Preferably use HepMC3, but alternatively HepMC2.
#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif

using namespace Pythia8;

//==========================================================================

int main() {

  // Interface for conversion from Pythia8::Event to HepMC event.
  // Specify file where HepMC events will be stored.
  Pythia8ToHepMC toHepMC("main135.hepmc");

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Process selection. LHC initialization. Histogram.
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Compress event record so that only final particles remain.
    // Need lines 0 - 2 to remain as "vertex" for the rest.
    int sizeOld = event.size();
    int sizeNew = 3;
    for (int i = 3; i < sizeOld; ++i) if (event[i].isFinal()) {
      event[sizeNew] = event[i];
      // Mother information not relevant, so blanked out.
      event[sizeNew].mothers( 1, 2);
      ++sizeNew;
    }
    event.popBack( sizeOld - sizeNew);
    event[1].daughters( 3, sizeNew - 1);
    event[2].daughters( 3, sizeNew - 1);

    // List first compressed event.
    if (iEvent == 0) event.list();

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 1; i < event.size(); ++i)
      if (event[i].isCharged()) ++nCharged;
    mult.fill( nCharged );

    // Construct new empty HepMC event, fill it and write it out.
    toHepMC.writeNextEvent( pythia );

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << mult;

  // Done.
  return 0;
}
