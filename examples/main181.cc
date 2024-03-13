// main181.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.m.utheim@jyu.fi>

// Keywords: PDFs; arXiv:2108.03481 [hep-ph]

// Plots parton distribution functions for the specified hadrons and partons.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// This struct defines which partons and hadrons should be shown in a plot.
struct PlotCase {
  vector<int> partonIds, hadronIds;
};

// These cases are the ones that will be plotted.
vector<PlotCase> plotCases = {
 { { 2, 1, -1, -2 }, { 2212, 2112 } },
 { { 2, 1 }, { 1114, 2114, 2214, 2224 } },
 { { 3, 2, 1 }, { 3112, 3212, 3222 } },
 { { 1, 2, -1, -2 }, { 211, 111, -211 } },
 { { 1, 2, 3 }, { 221, 331, 130, 310 } },
 { { 1, -1, 3, -3 }, { 311, 321 } },
 { { -1, -2, -3, 4 }, { 421, 431 } },
 { { 2, 3, 4, 5 }, { 521, 531, 541 } },
 { { 3, 4, 5 }, { 333, 443, 553 } },
};

// The Q^2 at which the PDFs are plotted.
constexpr double Q2 = 100.;

// Each curve is shifted up by a small amount to help separate the curves.
// Set this to 0 to plot the "correct" PDFs.
constexpr double dOffset = 0.002;
double offset;

// Flags that determine whether to show full/valence/sea content.
constexpr bool showFull = true, showValence = true, showSea = true;


int main() {

  // Initialization.
  Pythia pythia;
  pythia.readString("ProcessLevel:all = off");

  // Use the following settings to use SU21 PDF sets for p/pi:
  //pythia.readString("PDF:pSet = 24");
  //pythia.readString("PDF:piSet = 3");
  pythia.init();

  // Loop over and plot the cases.
  HistPlot plt("main181plot");
  for (auto plotCase : plotCases) {

    // Loop over each hadron and plot full contents.
    if (showFull) {
      plt.frame("out181plot", "PDF comparison", "x", "x f(x)");
      offset = 0.;
      for (int id : plotCase.hadronIds) {
        PDFPtr pdf = pythia.getPDFPtr(id);
        if (!pdf) {
          cout << "Failed to get PDF pointer for " << id << endl;
          return -1;
        }
        for (int q : plotCase.partonIds) {
          Hist h = Hist::plotFunc([&](double x) {
            return offset + pdf->xf(q, x, Q2);
            }, to_string(q) + "/" + to_string(id), 100, 0., 1.);
          plt.add(h, "-");
          offset += dOffset;
        }
      }
      plt.plot();
    }

    // Loop over each hadron and plot valence contents.
    if (showValence) {
      plt.frame("out181plot", "Valence PDF comparison", "x", "x v(x)");
      offset = 0.;
      for (int id : plotCase.hadronIds) {
        PDFPtr pdf = pythia.getPDFPtr(id);
        if (!pdf) {
          cout << "Failed to get PDF pointer for " << id << endl;
          return -1;
        }
        for (int q : plotCase.partonIds) {
          Hist h = Hist::plotFunc([&](double x) {
            return offset + pdf->xfVal(q, x, Q2);
            }, to_string(q) + "/" + to_string(id), 100, 0., 1.);
          plt.add(h, "-");
          offset += dOffset;
        }
      }
      plt.plot();
    }

    // Loop over each hadron and plot sea contents.
    if (showSea) {
      plt.frame("out181plot", "Sea PDF comparison", "x", "x s(x)");
      offset = 0.;
      for (int id : plotCase.hadronIds) {
        PDFPtr pdf = pythia.getPDFPtr(id);
        if (!pdf) {
          cout << "Failed to get PDF pointer for " << id << endl;
          return -1;
        }
        for (int q : plotCase.partonIds) {
          plt.add(Hist::plotFunc([&](double x) {
            return offset + pdf->xfSea(q, x, Q2);
          }, to_string(q) + "/" + to_string(id), 100, 0., 1.), "-");
          offset += dOffset;
        }
      }
      plt.plot();
    }

  }

  return 0;
}
