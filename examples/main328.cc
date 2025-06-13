// main328.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: total cross section; partial cross sections

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Study total, elastic and diffractive cross sections in several models.
// Most only available for pp and p pbar; and RPP only for total and elastic.
// ABMST has been modified for diffraction, see html manual for details.

// Note: This main program circumvents the need to generate loads of events
// to just output total cross sections, by creating a free standing object
// of type SigmaTotal, and copy relevant information from the Pythia object.
// See main486.cc for example how to use the public cross section methods.

#include "Pythia8/Pythia.h"
#include "Pythia8/SigmaTotal.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Basic parameters and histogram binning.
  int idA = 2212;
  int idB = 2212;
  double eCMmin = 10.;
  double eCMmax = 100000.;
  int nECM = 40;

  // Histograms.
  Hist sigTot[5], sigEl[5], sigND[5], sigSD[5], sigDD[5], sigCD[5];
  for (int i = 0; i < 5; ++i) {
    sigTot[i].book("Total cross section",       nECM, eCMmin, eCMmax, true);
    sigEl[i].book("Elastic cross section",      nECM, eCMmin, eCMmax, true);
    sigND[i].book("Non-diff cross section",     nECM, eCMmin, eCMmax, true);
    sigSD[i].book("Single diff cross section",  nECM, eCMmin, eCMmax, true);
    sigDD[i].book("Double diff cross section",  nECM, eCMmin, eCMmax, true);
    sigCD[i].book("Central diff cross section", nECM, eCMmin, eCMmax, true);
  }
  double eCM, sigT,sigE, sigN, sigS, sigD, sigC;

  // Loop over parametrizations.
  for (int iSig = 0; iSig < 4; ++iSig) {
    int modeSig = iSig + 1;

    // A SigmaTotal object, which will be used to calculate the cross
    // sections.
    SigmaTotal st;
    // Generator. Cross section model choice.
    Pythia pythia;
    pythia.settings.mode("Beams:idA", idA);
    pythia.settings.mode("Beams:idB", idB);
    pythia.settings.mode("SigmaTotal:mode", modeSig);
    pythia.settings.mode("SigmaDiffractive:mode", modeSig);

    // Switch off most of event generation and initialize.
    pythia.readString("SoftQCD:elastic = on");
    pythia.readString("PartonLevel:all = off");
    pythia.readString("HadronLevel:all = off");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Make a copy of the Pythia info object to initialize
    // pointers in SigmaTotal.
    Info pythiaInfo(pythia.info);
    // Initialize the pointers.
    st.initInfoPtr(pythiaInfo);
    // Initialize the SigmaTotal instance.
    st.init();

    // Loop through logarithmic grid of collision energies.
    for (int iECM = 0; iECM < nECM; ++iECM) {
      eCM = eCMmin * pow( eCMmax / eCMmin, (iECM + 0.5) / double(nECM) );

      // Calculate cross sections.
      st.calc( idA, idB, eCM);
      sigT = st.sigmaTot();
      sigE = st.sigmaEl();
      if (iSig < 3) {
        sigN = st.sigmaND();
        sigS = st.sigmaXB() + st.sigmaAX();
        sigD = st.sigmaXX();
        sigC = st.sigmaAXB();
      }

      // Fill histograms.
      sigTot[iSig].fill( eCM, sigT);
      sigEl[iSig].fill(  eCM, sigE);
      if (iSig < 3) {
        sigND[iSig].fill(  eCM, sigN);
        sigSD[iSig].fill(  eCM, sigS);
        sigDD[iSig].fill(  eCM, sigD);
        sigCD[iSig].fill(  eCM, sigC);
      }
    }

  // End loop over parametrizations.
  }

  // Plot histograms.
  HistPlot hpl("plot328");
  hpl.frame("fig328", "Total cross section in pp collisions",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigTot[0], "-,black", "SaS/DL");
  hpl.add( sigTot[1], "-,red", "MBR");
  hpl.add( sigTot[2], "-,blue", "ABMST");
  hpl.add( sigTot[3], "-,green", "RPP2016");
  hpl.plot();
  hpl.frame("", "Elastic cross section in pp collisions",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigEl[0], "-,black", "SaS/DL");
  hpl.add( sigEl[1], "-,red", "MBR");
  hpl.add( sigEl[2], "-,blue", "ABMST");
  hpl.add( sigEl[3], "-,green", "RPP2016");
  hpl.plot();
  hpl.frame("", "Nondiffractive cross section in pp collisions",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigND[0], "-,black", "SaS/DL");
  hpl.add( sigND[1], "-,red", "MBR");
  hpl.add( sigND[2], "-,blue", "ABMST mod");
  hpl.plot();
  hpl.frame("", "Single diffractive cross section in pp collisions",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigSD[0], "-,black", "SaS/DL");
  hpl.add( sigSD[1], "-,red", "MBR");
  hpl.add( sigSD[2], "-,blue", "ABMST mod");
  hpl.plot();
  hpl.frame("", "Double diffractive cross section in pp collisions",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigDD[0], "-,black", "SaS/DL");
  hpl.add( sigDD[1], "-,red", "MBR");
  hpl.add( sigDD[2], "-,blue", "ABMST mod");
  hpl.plot();
  hpl.frame("", "Central diffractive cross section in pp collisions",
    "$E_{\\mathrm{CM}}$ (GeV)", "$\\sigma$ (mb)", 8.0, 5.4);
  hpl.add( sigCD[0], "-,black", "SaS/DL");
  hpl.add( sigCD[1], "-,red", "MBR");
  hpl.add( sigCD[2], "-,blue", "ABMST mod");
  hpl.plot();

  // Done.
  return 0;
}
