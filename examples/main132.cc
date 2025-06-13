// main132.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Mikhail Kirsanov <Mikhail.Kirsanov@cern.ch>

// Keywords: basic usage; command file; command line option; hepmc; LHE file

// This program illustrates how HepMC files can be written by Pythia8.
// Input and output files are specified on the command line, e.g. like
//     ./main132 -c main132.cmnd -o main132.hepmc > main132.log
// Either internal Pythia processes or Les Houches Event Files can be used.
// The main program contains no analysis; this is intended to happen later.
// It therefore "never" has to be recompiled to handle different tasks.

// WARNING: typically one needs 25 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"
#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif

using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Set up command line options.
  InputParser ip("This program illustrates how HepMC files can be written by"
    " Pythia8.", {"./main132 -c main132.cmnd -o main132.hepmc"});
  ip.require("c", "Use this user-written command file.", {"-cmnd"});
  ip.require("o", "Specify HepMC output filename.", {"-out"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Confirm that external files will be used for input and output.
  string cmnd(ip.get<string>("c")), out(ip.get<string>("o"));
  cout << "\n >>> PYTHIA settings will be read from file '" << cmnd
       << "' <<< \n >>> HepMC events will be written to file '"
       << out << "' <<< \n";

  // Interface for conversion from Pythia8::Event to HepMC event.
  // Specify file where HepMC events will be stored.
  Pythia8ToHepMC toHepMC(out);

  // Generator.
  Pythia pythia;

  // Read in commands from external file.
  pythia.readFile(cmnd);

  // Extract settings to be used in the main program.
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate event.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) {
        cout << " Aborted since reached end of Les Houches Event File\n";
        break;
      }

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Construct new empty HepMC event, fill it and write it out.
    toHepMC.writeNextEvent( pythia );

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Done.
  return 0;
}
