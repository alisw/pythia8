// main114.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// Keywords: basic usage; plotting; yoda;

// This is a simple test program.
// It illustrates how YODA2 histograms can be used through the Pythia8Yoda
// simplified interface. This can be useful for more advanced statistical
// treatments than the built-in histogram types allow, for example using
// profile histograms or custom types, as shown below.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/Pythia8Yoda.h"

using namespace Pythia8;

//==========================================================================

int main() {
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("SoftQCD:all = on");
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;
  Pythia8Yoda p8y("MAIN114", "main114.yoda");
  // Booking of normal types, as defined in Pythia8Yoda.
  auto h1 = p8y.bookHisto1D(15, 0.5, 191.5, "charged-mult");
  auto h2 = p8y.bookProfile1D(15, 0.5, 191.5, "mean-pt");
  // Booking of custom YODA types. Two examples:
  // 1) a 1D histogram with discrete string binning, used
  // to count the number of pi, K and p per event, and
  // 2) a 1D Profile histogram with discrete string binning,
  // used to count <pT> of said particles.
  // First the types:
  using IntegerHisto1D = YODA::BinnedHisto1D<string>;
  using IntegerProfile1D = YODA::BinnedProfile1D<string>;
  // Then the desired axis binning.
  const auto pidAxis = YODA::Axis<string>({"211", "321", "2212"});
  // And finally the two histograms.
  auto h3 = p8y.book<IntegerHisto1D>(pidAxis,
    "/MAIN114/pid-mult", "pid-mult");
  auto h4 = p8y.book<IntegerProfile1D>(pidAxis,
    "/MAIN114/pid-mean-pt", "pid-mean-pt");

  // Begin event loop. Generate event. Skip if error. List first one.
  int nEvents = 10000;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;
    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    double sumpT = 0.;
    for (auto p : pythia.event) {
      if (p.isFinal() && p.isHadron() && p.isCharged() &&
        abs(p.eta()) < 2.5) {
        ++nCharged;
        h3->fill(to_string(abs(p.id())));
        sumpT += p.pT();
        h4->fill(to_string(abs(p.id())), p.pT());
      }
    }
    if (nCharged > 0) {
      h1->fill(nCharged);
      h2->fill(nCharged, sumpT/double(nCharged));
    }
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  // Normalize to number of events generated.
  h1->scaleW(1./double(nEvents));
  h3->scaleW(1./double(nEvents));

  // Write histograms to file.
  p8y.write();
  cout << "Histograms have been written. "
    "You can now plot them with e.g. rivet-mkhtml or Python.\n";

  // Done.
  return 0;
}
