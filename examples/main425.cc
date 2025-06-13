// main425.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim
//          Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: heavy ions, cross sections.

// This example calculates proton-oxygen cross sections at varying
// energies using the Angantyr module in Pythia. Illustrates the
// difference between the "generated" cross section and the cross
// section calculated from the Glauber modelling.

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;
//==========================================================================

int main() {

  // Input in fixed-target or CM frame.
  bool useFixed   = true;
  // Reuse precalculated initialization data or not.
  bool usePrecalc = true;

  // Set up momentum grid for fixed-target option.
  double pMin = 1e2, pMax = 1e11;
  int nPts = 4;
  vector<double> pLabs = logSpace(nPts, pMin, pMax);
  double dr = pow(pMax / pMin, 1. / (nPts - 1));

  // Histograms.
  Hist sigTotAn("Total", nPts, pMin / sqrt(dr), pMax * sqrt(dr), true);
  Hist sigInelAn("Inelastic", nPts, pMin / sqrt(dr), pMax * sqrt(dr), true);
  Hist sigTotPy("Total", nPts, pMin / sqrt(dr), pMax * sqrt(dr), true);
  Hist sigInelPy("Inelastic", nPts, pMin / sqrt(dr), pMax * sqrt(dr), true);

  // Iterate over momenta. Initialize for p 16O(xygen).
  for (double pNow : pLabs) {
    Pythia pythia;
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 1000080160");
    pythia.readString("HadronLevel:all = off");

    // Set up in fixed-target frame.
    if (useFixed) {
      pythia.readString("Beams:frameType = 3");
      pythia.settings.parm("Beams:pzA", pNow);
      pythia.settings.parm("Beams:pzB", 0.);

    // Alternatively initialize for equivalent proton-nucleon CM energy.
    } else {
      pythia.readString("Beams:frameType = 1");
      double eCMNow = ( Vec4(0., 0., pNow, pNow * sqrt(1 + pow2(0.938 / pNow)))
                      + Vec4(0., 0., 0., 0.938) ).mCalc();
      pythia.settings.parm("Beams:eCM", eCMNow);
    }

    // Optionally reuse initialization information (if it exists, see main424).
    if (usePrecalc) {
      pythia.readFile("main424.cmnd");
    }

    // Initialize.
    if (!pythia.init()) {
      cout << "Pythia failed to initialize." << endl;
      return -1;
    }

    // Generate events to get statistics.
    for (int iEvent = 0; iEvent < 10000; ++iEvent) {
      pythia.next();
      if (iEvent == 0) pythia.event.list();
    }

    // Read out total and inelastic cross section two ways:

    // First we have the full total and inelastic non-diffractive
    // cross for p-O, as obtained from the Glauber caclulation.
    double hiTot  = pythia.info.hiInfo->glauberTot();
    double hiInel = pythia.info.hiInfo->glauberND();

    // Then we have the cross section of the generated events
    // categorised according to the type primary p-nucleon
    // scattering. Note that since Angantyr does not actually generate
    // proper inelastic events, not even the total cross section is
    // the same as the total p-O cross section, but also the inelastic
    // non-diffractive may differ, since also a diffractive primary
    // p-nucleon scattering may correspond to a non-diffractive p-O
    // scattering.
    double pyTot  = pythia.info.sigmaGen();
    double pyInel = pythia.info.sigmaGen(101);

    // Fill histograms with cross sections above.
    sigTotAn.fill(pNow, hiTot);
    sigInelAn.fill(pNow, hiInel);
    sigTotPy.fill(pNow, pyTot);
    sigInelPy.fill(pNow, pyInel);

    // Print statistics.
    pythia.stat();
  }

  // Print histograms.
  cout << sigTotAn << sigInelAn << sigTotPy << sigInelPy;

  // Plot histograms.
  HistPlot plt("plot425");
  plt.frame("fig425", "p ${}^{16}$O cross sections",
    "$p_{Lab}$ (GeV)", "$\\sigma$ (mb)", 6.4, 4.8);
  plt.add(sigTotAn, "-", "Total p-O from Glauber");
  plt.add(sigInelAn, "-", "Inelastic non-diffactive p-O from Glauber");
  plt.add(sigTotPy, "-", "Total generated");
  plt.add(sigInelPy, "--", "Generated with non-diffractive primary p-N");
  plt.plot(0.5 * pMin, 2. * pMax, 0., 1200., false, true);

  return 0;
}
