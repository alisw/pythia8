// main144.C is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// ROOT macro demonstrating how one can read in ROOT files generated with
// main144.

// The following steps describe how this macro can be run.
//
// (1) Configure Pythia with ROOT.
//         ./configure --with--root
// (2) Compile main144.
//         make main144
//     This build command triggers the following steps, which are
//     outline here for clarity.
//     (a) Generate source for ROOT dictionary (main144Dct.cc) with CLING.
//             rootcling -f main144Dct.cc -c main144Dct.h
//     (b) Compile the ROOT dictionary library, main144Dct.so.
//     (c) Compile main144, linked against main144Dct.so.
// (3) Run main144, with the ROOT command options. This generates
//     the output main144.root.
//         ./main144 -c main144.cmnd -c main144Root.cmnd
// (4) Run this macro with ROOT.
//         root main144.C
//     This will genererate two plots, one of the track pT and the other of
//     the pion pseudo-rapidity
//
// Alternatively, step (4) can be replaced with an interactive ROOT session.
//     root [0] .L main144.C
//     root [1] read({"main144.root"}) 
//
// Finally, one can build a shared reader library from the macro.
//     root main144.C+
// The + symbol here directs ROOT to compile a shared library from the
// macro, main144Reader_C.so, which can then be linked to compiled
// code. In this code, the read method must be declared (but not
// defined).

// Generally ROOT headers do not need to be included in macros.
// These are included so that the macro can be compiled as a library.
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"

// The ROOT event.
#include "main144Dct.h"

// Load the Pythia ROOT event library.
R__LOAD_LIBRARY(main144Dct.so)

// Read one or more files, loop over events, draw histograms.
void read(std::vector<std::string> filenames = {"main144.root"}) {

  // Create a TChain and add the files (one TTree per file).
  TChain* tree = new TChain("t");
  for (auto &name : filenames) tree->Add(name.c_str());
  // Map the TChain to a RootEvent pointer.
  RootEvent* evt = nullptr;
  tree->SetBranchAddress("events", &evt);
  
  // Create some histograms to fill.
  TH1D* hPt = new TH1D("hpT", "all particle p_{T} [GeV]", 100, 0.0, 10.0);
  TH1D* hEtaPi = new TH1D("hEtaPi", "charged #pi #eta", 20, -10, 10.0);

  // Loop over the events.
  int nEntries = tree->GetEntries();
  for (int iEntry = 0; iEntry < nEntries; iEntry++) {

    // Retrieve the event from the TChain. Now evt will contain the event.
    tree->GetEntry(iEntry);

    // Grab the event weight.
    double w = evt->weight;

    // Loop over the particles.
    for (RootParticle &prt : evt->particles) {
      hPt->Fill(prt.pT, w);
      if (abs(prt.pid) == 211) hEtaPi->Fill(prt.eta, w); 
    }
  }
  
  // (4) Draw the histograms.
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  hPt->Draw();
  TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  hEtaPi->Draw();
  
}

// Define the main execution for this ROOT macro.
void main144Reader() {

  // Read the output of main144, which is produced with the following.
  //     ./main144 -c main144.cmnd -c main144Root.cmnd
  read({"main144.root"});

}
