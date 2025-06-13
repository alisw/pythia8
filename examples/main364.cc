// main364.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Philip Ilten <philten@cern.ch>

// Keywords: B decays; external decays; EvtGen

// An example where decays are performed with EvtGen. See the
// documentation in Pythia8Plugins/EvtGen.h for further details on the
// EvtGenDecays class. In this example the decay B0 -> nu_e e+ D*-[->
// gamma D-[-> nu_ebar e- pi0]] is forced. The invariant mass spectrum
// of the final state electron and pion is then plotted.

// The syntax to run this example is:
//     ./main364 -d <EvtGen decay file> -p <EvtGen particle data>
//               -x <PYTHIA8DATA>
//               <flag to use EvtGen>

// This example has only been tested with EvtGen version 1.3.0. The
// EvtGen package is designed to use Pythia 8 for any decays that it
// cannot perform. For this to be possible, EvtGen must be linked
// against the Pythia 8 shared library. To modify how this example
// program is compiled (i.e to remove linking against the
// EvtGenExternal library) change the main364 rule in the Makefile of
// this directory.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"
#include "Pythia8Plugins/EvtGen.h"
using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Set up command line options.
  InputParser ip("Optionally perform decays with EvtGen.",
    {"./main364 -d <EvtGen decay file> -p <EvtGen particle data>"
        "\n\t          -x <PYTHIA8DATA> -e"});
  ip.require("d", "EvtGen decay file (e.g. DECAY_2010.DEC).", {"-dec"});
  ip.require("p", "EvtGen particle data (e.g. evt.pdl).", {"-pdl"});
  ip.require("x", "PYTHIA8DATA path.", {"-xml"});
  ip.add("e", "false", "Flag to use EvtGen.");

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Get the options.
  string dec = ip.get<string>("d");
  string pdl = ip.get<string>("p");
  string xml = ip.get<string>("x");
  bool use = ip.get<bool>("e");

  // Intialize Pythia.
  Pythia pythia;
  pythia.readString("Print:quiet = on");
  pythia.readString("HardQCD:hardbbbar = on");
  if (!use) {
    cout << "Not using EvtGen." << endl;
    pythia.readString("511:onMode = off");
    pythia.readString("511:onIfMatch = 12 -11 -413");
    pythia.readString("413:onMode = off");
    pythia.readString("413:onIfMatch = 411 22");
    pythia.readString("411:onMode = off");
    pythia.readString("411:onIfMatch = -11 12 111");
  } else cout << "Using EvtGen." << endl;

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Initialize EvtGen.
  EvtGenDecays *evtgen = 0;
  if (use) {
    setenv("PYTHIA8DATA", xml.c_str(), 1);
    evtgen = new EvtGenDecays(&pythia, dec, pdl);
    evtgen->readDecayFile("main364.dec");
  }

  // The event loop.
  Hist mass("m(e, pi0) [GeV]", 100, 0., 2.);
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    // Generate the event.
    if (!pythia.next()) continue;
    // Perform the decays with EvtGen.
    if (evtgen) evtgen->decay();
    // Analyze the event.
    Event &event = pythia.event;
    for (int iPrt = 0; iPrt < (int)event.size(); ++iPrt) {
      if (event[iPrt].idAbs() != 511) continue;
      int iB0(event[iPrt].iBotCopyId()), iDsm(-1), iDm(-1), iE(-1), iPi0(-1);
      for (int iDtr = event[iB0].daughter1(); iDtr <= event[iB0].daughter2();
           ++ iDtr) {
        if (event[iDtr].idAbs() == 413) {
          iDsm = event[iDtr].iBotCopyId();
          continue;
        }
      }
      if (iDsm == -1) continue;
      for (int iDtr = event[iDsm].daughter1(); iDtr <= event[iDsm].daughter2();
           ++ iDtr) {
        if (event[iDtr].idAbs() == 411) {
          iDm = event[iDtr].iBotCopyId();
          continue;
        }
      }
      if (iDm == -1) continue;
      for (int iDtr = event[iDm].daughter1(); iDtr <= event[iDm].daughter2();
           ++ iDtr) {
        if (event[iDtr].idAbs() == 11)  iE   = event[iDtr].iBotCopyId();
        if (event[iDtr].idAbs() == 111) iPi0 = event[iDtr].iBotCopyId();
      }
      if (iE == -1 || iPi0 == -1) continue;
      mass.fill((event[iE].p() + event[iPi0].p()).mCalc());
    }
  }

  // Print the statistics and histogram.
  pythia.stat();
  mass /= mass.getEntries();
  cout << mass;
  if (evtgen) delete evtgen;
  return 0;
}
