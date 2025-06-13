// main143.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Rene Brun, Axel Naumann and Bernhard Meirose

// Keywords: analysis; root

// This is a simple test program, based on main101.cc,
// but modified to use ROOT for histogramming.
// It studies the charged multiplicity distribution at the LHC.

// WARNING: for currently unknown reasons it may hang
// with an empty canvas on a Mac.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Beams:eCM = 14000.");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("hist.root", "RECREATE");

  // Book histogram.
  TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);

  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;

    // Fill charged multiplicity in histogram. End event loop.
    mult->Fill( nCharged );
  }

  // Statistics on event generation.
  pythia.stat();

  // Show histogram. Possibility to close it.
  mult->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();

  // Save histogram on file and close file.
  mult->Write();
  delete outFile;

  // Done.
  return 0;
}
