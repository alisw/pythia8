// main486.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: total cross section; partial cross sections

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Study total and inelastic cross section for various beam combinations,
// using the public methods in the Pythia class. These methods are intended
// for fast switching, and only provide the SaS/DL ansats at high energies.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Histogram binning. Beam combinations.
  double eCMmin = 2.;
  double eCMmax = 100000.;
  int nECM = 80;
  int idA[5] = { 2212,  2212,  211,  321,  211};
  int idB[5] = { 2212, -2212, 2212, 2212, -211};

  // Histograms.
  Hist sigTot[5], sigInel[5];
  for (int i = 0; i < 5; ++i) {
    sigTot[i].book("Total cross section",       nECM, eCMmin, eCMmax, true);
    sigInel[i].book("Inelastic cross section",  nECM, eCMmin, eCMmax, true);
  }
  double eCM, sigT, sigI;

  // Generator. Need minimal setup to get going with cross sections only.
  Pythia pythia;
  pythia.readString("SoftQCD:elastic = on");
  pythia.readString("PartonLevel:all = off");
  pythia.readString("HadronLevel:all = off");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Loop over beam combinations.
  for (int iSig = 0; iSig < 5; ++iSig) {

    // Loop through logarithmic grid of collision energies.
    for (int iECM = 0; iECM < nECM; ++iECM) {
      eCM = eCMmin * pow( eCMmax / eCMmin, (iECM + 0.5) / double(nECM) );

      // Calculate cross sections. Inelastic = total - elastic.
      sigT = pythia.getSigmaTotal( idA[iSig], idB[iSig], eCM);
      sigI = sigT - pythia.getSigmaPartial( idA[iSig], idB[iSig], eCM, 2);

      // Fill histograms.
      sigTot[iSig].fill(  eCM, sigT);
      sigInel[iSig].fill( eCM, sigI);

    // End loops over energies and beams.
    }
  }

  // Plot histograms.
  HistPlot hpl("plot486");
  hpl.frame("fig486", "Total cross sections",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigTot[0], "-,black", "pp");
  hpl.add( sigTot[1], "--,black", "p$\\overline{\\mathrm{p}}$");
  hpl.add( sigTot[2], "-,blue", "$\\pi^+$p");
  hpl.add( sigTot[3], "-,green", "K$^+$p");
  hpl.add( sigTot[4], "-,red", "$\\pi^+\\pi^-$");
  hpl.plot();
  hpl.frame("", "Inelastic cross sections",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigInel[0], "-,black", "pp");
  hpl.add( sigInel[1], "--,black", "p$\\overline{\\mathrm{p}}$");
  hpl.add( sigInel[2], "-,blue", "$\\pi^+$p");
  hpl.add( sigInel[3], "-,green", "K$^+$p");
  hpl.add( sigInel[4], "-,red", "$\\pi^+\\pi^-$");
  hpl.plot();

  // Done.
  return 0;
}
