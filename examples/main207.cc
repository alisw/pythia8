// main207.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Peter Skands <peter.skands@monash.edu>
// Thanks to: M. di Mauro, for the original template for this example.

// Keywords: Vincia; weak showers; LHEF; dark matter;

// Example showing how to run Vincia's electroweak shower with LHEF input.
// This requires the LHEF files to contain helicities for the hard partons.
// In this example, the LHEF file contains dark-matter particles annihilating
// to electron-positron pairs.

// Note: emitted weak bosons decay inclusively; it would be up to the user
// themselves to filter events with decays to specific channels if desired.

#include "Pythia8/Pythia.h"
#include "math.h"
#include <iostream>
#include <string>

using namespace std;
using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;

  // Read Pythia settings from command file.
  pythia.readFile("main207.cmnd");

  // Initialization.
  pythia.init();

  // Allow for possibility of a few faulty events.
  int nError      = pythia.settings.mode("Main:timesAllowErrors");

  double logxMin = -9.;
  double logxMax = 0.;
  double nBins = 100;
  double DeltaBin = (logxMax-logxMin)/nBins;

  // Histogram particle spectra.
  Hist gamma("gamma spectrum", nBins, logxMin, logxMax);
  Hist electron("e+- spectrum", nBins, logxMin, logxMax);
  Hist proton("p spectrum", nBins, logxMin, logxMax);
  Hist nue("nu_e spectrum", nBins, logxMin, logxMax);
  Hist numu("nu_mu spectrum", nBins, logxMin, logxMax);
  Hist nutau("nu_tau spectrum", nBins, logxMin, logxMax);
  Hist rest("remaining particle spectrum", nBins, logxMin, logxMax);

  // Begin event loop.
  int nEvent = 0;
  int iError = 0;
  while (true) {

    // Generate the next event.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // Otherwise count event failure and continue/exit as necessary.
      cout << "Warning: event " << nEvent << " failed" << endl;
      if (++iError == nError) {
        cout << "Error: too many event failures... exiting" << endl;
        break;
      }

      continue;
    }

    /*
     * Process dependent checks and analysis go here ...
     */

    // Add total number of events attempted.
    ++nEvent;

    // Get DM mass.
    double mDM = pythia.event[1].m();

    // Loop over all particles and select final-state ones for histrograms.
    for (int i = 0; i < pythia.event.size(); ++i) {
      if ( !pythia.event[i].isFinal() ) continue;
      int id       = pythia.event[i].id();
      int idAbs    = pythia.event[i].idAbs();
      double e     = pythia.event[i].e();
      if (e <= 0) continue;
      double logx  = log10(e/mDM);
      // Select photons.
      if (idAbs == 22) gamma.fill(logx);
      // Select electrons (positron equivalent).
      else if (id == -11)  electron.fill(logx);
      // Select protons.
      else if (id == -2212 or idAbs == 2112) proton.fill(logx);
      // Select various neutrinos.
      else if (id == 12) nue.fill(logx);
      else if (id == 14) numu.fill(logx);
      else if (id == 16) nutau.fill(logx);
      else rest.fill(logx);
    }
  }


  //Statistic and histrograms.
  pythia.stat();

  gamma.operator*=(1./nEvent/DeltaBin);
  electron.operator*=(1./nEvent/DeltaBin);
  proton.operator*=(1./nEvent/DeltaBin);
  nue.operator*=(1./nEvent/DeltaBin);
  numu.operator*=(1./nEvent/DeltaBin);
  nutau.operator*=(1./nEvent/DeltaBin);

  // Make tables, including statistical errors (last argument = true).
  gamma.table("main207-gamma.dat", false, true, true);
  electron.table("main207-positron.dat", false, true, true);
  proton.table("main207-antiproton.dat", false, true, true);
  nue.table("main207-nue.dat", false, true, true);
  numu.table("main207-numu.dat", false, true, true);
  nutau.table("main207-nutau.dat", false, true, true);
  rest.table("main207-rest.dat", false, true, true);

  cout << gamma << electron << proton << nue << numu << nutau << rest;

  //Done
  return 0;
}
