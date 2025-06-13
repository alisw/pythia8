// main465.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim

// Keywords: rescattering; low energy; cross sections; resonances

// Calculate and plot resonance cross sections for the specified process.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"
using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Set up command line options.
  InputParser ip("Calculate resonance cross section for specified process.",
    {"./main465 -a 2212 -b 2212"});
  ip.add("a", "0", "ID for beam A.", {"-ida"});
  ip.add("b", "0", "ID for beam B.", {"-idb"});
  ip.add("n", "300", "Number of bins.", {"-nbins"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Initialize Pythia.
  Pythia pythia;
  pythia.readFile("main465.cmnd");
  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl;
    return 1;
  }

  int idA = ip.get<int>("a");
  int idB = ip.get<int>("b");
  if (idA == 0) idA = pythia.mode("Main:spareMode1");
  if (idB == 0) idB = pythia.mode("Main:spareMode2");
  double eMin = pythia.parm("Main:spareParm1");
  double eMax = pythia.parm("Main:spareParm2");
  if (eMin < pythia.particleData.m0(idA) + pythia.particleData.m0(idB)) {
    eMin = pythia.particleData.m0(idA) + pythia.particleData.m0(idB);
    cout << "Warning, setting eMin to nominal mass sum of " << eMin << ".\n";
  }
  int nBin = ip.get<int>("n");

  // Get possible resonances.
  set<int> resonances = pythia.hadronWidths.getResonances(idA, idB);

  if (resonances.size() == 0) {
    cout << "No resonances for input particles " << idA << " " << idB << endl;
    return -1;
  }

  // Define plot.
  HistPlot plt("plot465");
  plt.frame("fig465", "", "$\\sqrt{s}$ [GeV]", "$\\sigma$ [mb]");

  // Plot all resonances.
  for (int res : resonances) {
    Hist sigRes = Hist::plotFunc( [&](double eCM) {
      return pythia.getSigmaPartial(idA, idB, eCM, res);
    }, pythia.particleData.name(res), nBin, eMin, eMax);

    plt.add(sigRes, "-");
  }

  // Add total and miscellaneous cross sections at the end.
  Hist sigTot = Hist::plotFunc(
    [&](double eCM) {
      // type == 0 is equivalent to getSigmaTotal.
      return pythia.getSigmaPartial(idA, idB, eCM, 0); },
    "Total", nBin, eMin, eMax);
  Hist sigTotRes = Hist::plotFunc(
    [&](double eCM) { return pythia.getSigmaPartial(idA, idB, eCM, 9); },
    "Sum of resonances", nBin, eMin, eMax);
  Hist sigOther = sigTot - sigTotRes;

  plt.add(sigTotRes, "-");
  plt.add(sigOther, "-", "Other");
  plt.add(sigTot, "k-");

  // Plot.
  plt.plot();

  // Done.
  return 0;
}
