// main122.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; LHE file; matplotlib

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in PYTHIA.
// It uses two LHE files, ttbar.lhe and ttbar2.lhe, which are combined
// using Beams:newLHEFsameInit = on to skip new initialization second time.
// Then the second file is viewed as a simple continuation of the first,
// just split for practical reasons, rather than as a separate new run
// with a new set of processes.
// In the first file top decays have been performed, in the second not,
// and are instead handled by the internal PYTHIA resonance-decay machinery.
// Furthermore the internal top production processes are switched on and
// mixed in, giving an unrealistic "double up" total top cross section.
// Much of this of course is not intended to be realistic,
// but rather illustrates several tricks that can be useful.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  //  Number of listed events. Allow for possibility of a few faulty events.
  int nPrintLHA  = 1;
  int nPrintRest = 0;
  int nAbort     = 10;

  // Generator
  Pythia pythia;

  // Switch on internal ttbar production.
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");

  // Use same top mass as in LHE files to simplify comparison.
  pythia.readString("6:m0 = 175.");

  // No automatic event listings - do it manually below.
  pythia.readString("Next:numberShowLHA = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize Les Houches Event File run.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = ttbar.lhe");
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Book histogram.
  Hist nChargedL("nch, LHEF events",100,-0.5,399.5);
  Hist nChargedI("nch, internal events",100,-0.5,399.5);

  // Set counters.
  int iPrintLHA  = 0;
  int iPrintRest = 0;
  int iAbort     = 0;
  int iFile      = 1;

  // Begin event loop
  while (iAbort < nAbort) {

    // Generate until none left in input file.
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) {

        // First time open next file, second time stop event loop.
        if (iFile == 1) {
          pythia.readString("Beams:newLHEFsameInit = on");
          pythia.readString("Beams:LHEF = ttbar2.lhe");

          // If Pythia fails to initialize, exit with error.
          if (!pythia.init()) {
            cout << "Failed init() at file " << iFile << endl;
            return 1;
          }

          ++iFile;
          continue;
        } else break;
      }
      ++iAbort;
      continue;
    }

    // List first few Les Houches and other events.
    if (pythia.info.isLHA() && iPrintLHA < nPrintLHA) {
      pythia.LHAeventList();
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
      ++iPrintLHA;
    } else if (!pythia.info.isLHA() && iPrintRest < nPrintRest) {
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
      ++iPrintRest;
    }

    // Sum up final charged multiplicity and fill in histogram.
    int nChg = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
      ++nChg;
    if (pythia.info.isLHA()) nChargedL.fill(nChg);
    else                     nChargedI.fill(nChg);

  // End of event loop.
  }

  // Give statistics. Print histograms.
  pythia.stat();
  cout << nChargedL << nChargedI;

  // Plot histograms.
  HistPlot hpl("plot122");
  hpl.frame("fig122", "Charged particle multiplicities",
     "$n_{\\mathrm{ch}}$", "$\\mathrm{d}P/\\mathrm{d}n_{\\mathrm{ch}}$");
  hpl.add(nChargedL);
  hpl.add(nChargedI);
  hpl.plot();

  // Done.
  return 0;
}
