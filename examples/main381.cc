// main381.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: Higgs; electron-positron;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Simple example of Higgs pruduction at future e+e- colliders.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent = 25000;

  // Generator. Incoming beams. (Switch off iniial-state photon radiation.)
  Pythia pythia;
  pythia.readString("Beams:idA = -11");
  pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:eCM = 500.");
  //pythia.readString("PDF:lepton = off");

  // All Higgs production channels.
  pythia.readString("HiggsSM:all = on");
  //pythia.readString("23:onMode = off");
  //pythia.readString("23:onIfAny = 1 2 3 4 5");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  Hist mult("charged multiplicity all", 50, -1., 99.);
  Hist multS("charged multiplicity Higgsstrahlung", 50, -1., 99.);
  Hist multG("charged multiplicity gammagamma fusion", 50, -1., 99.);
  Hist multZ("charged multiplicity ZZ fusion", 50, -1., 99.);
  Hist multW("charged multiplicity WW fusion", 50, -1., 99.);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
        ++nCharged;
    mult.fill( nCharged );

    // Split by process type.
    int code = pythia.info.code();
    if      (code == 904) multS.fill( nCharged );
    else if (code == 903) multG.fill( nCharged );
    else if (code == 906) multZ.fill( nCharged );
    else if (code == 907) multW.fill( nCharged );
    else cout << " Warning: unexpected process code " << code << endl;

  // End of event loop. Statistics. Normalize and print histograms.
  }
  pythia.stat();
  mult  /= double(nEvent);
  multS /= double(nEvent);
  multG /= double(nEvent);
  multZ /= double(nEvent);
  multW /= double(nEvent);
  cout << mult << multS << multG << multZ << multW;

  //Plot histograms.
  HistPlot hpl("plot381");
  hpl.frame( "fig381", "Charged multiplicity for Higgs production "
    "in a 500 GeV e$^+$e$^-$ collider", "$n_{\\mathrm{charged}}$",
    "Probability");
  hpl.add( mult, "h,black", "all events");
  hpl.add( multS, "h,red", "Higgsstrahlung");
  hpl.add( multW, "h,blue", "WW fusion");
  hpl.add( multZ, "h,olive", "ZZ fusion");
  hpl.add( multG, "h,magenta", "$\\gamma\\gamma$ fusion");
  hpl.plot();

  // Done
  return 0;
}
