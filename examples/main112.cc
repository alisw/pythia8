// main112.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; Z production; Tevatron; matplotlib

// This is a simple test program, equivalent with main102,
// for the pT_Z spectrum at the Tevatron, but with plotting of histogram.
// Assuming that you have Python3 installed on your platform, do as follows.
// After the program has run, type "python3 plot112.py" (without the " ")
// in a terminal window, and open "fig112.pdf" in a PDF viewer.


#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Generator. Process selection.
  Pythia pythia;
  pythia.readString("Beams:idB = -2212");
  pythia.readString("Beams:eCM = 1960.");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("PhaseSpace:mHatMin = 80.");
  pythia.readString("PhaseSpace:mHatMax = 120.");

  // Tevatron initialization. Histogram.
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;
  Hist pTZ("dN/dpTZ", 100, 0., 100.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    // Loop over particles in event. Find last Z0 copy. Fill its pT.
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 23) iZ = i;
    pTZ.fill( pythia.event[iZ].pT() );

  // End of event loop. Statistics. Output normalized histogram.
  }
  pythia.stat();
  pTZ /= 1000;
  cout << pTZ;

  // Write Python code that can generate a PDF file with the spectrum.
  // Use plotFrame when a frame should only contain one histogram.
  // Colours and other choices can be omitted, but are shown as illustration.
  HistPlot hpl("plot112");
  hpl.plotFrame("fig112", pTZ, "Z transverse momentum", "$p_{\\perp}$ (GeV)",
    "$\\mathrm{d}N/\\mathrm{d}p_{\\perp}$ (1/GeV)", "-");

  // Done
  return 0;
}
