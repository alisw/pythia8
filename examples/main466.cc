// main466.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim

// Keywords: hadron widths

// Code to create parameterization tables for hadron widths.
// Useful if resonances are added or particle properties are changed.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"
using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Set up command line options.
  InputParser ip("Create parameterization tables for hadron widths.",
    {"./main465 -a 2212 -b 2212"});
  ip.add("p", "50", "Precision, provided as integer.", {"-precision"});
  ip.add("o", "main466.dat", "Output file for the widths.", {"-out"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Initialize Pythia.
  Pythia pythia;
  pythia.readString("ProcessLevel:all = off");

  if (!pythia.init()) {
    cout << endl << " Pythia failed to initialize. \n"
     " If this happened because hadron widths are unavailable or invalid,\n"
     " particle data should still be loaded. In this case, this code should\n"
     " still be used to generate hadron widths. Therefore, execution will\n"
     " continue." << endl << endl;
  }

  // Perform parameterization.
  HadronWidths& hadronWidths = pythia.hadronWidths;
  hadronWidths.parameterizeAll(ip.get<int>("p"));
  hadronWidths.save(ip.get<string>("o"));

  // Done.
  return 0;

}
