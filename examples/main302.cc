// main302.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: colour reconnection; e+e- events;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Plot the colour reconnection rate in W+W- hadronic events
// as a function of collision energy.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent = 4000;

  // Histograms.
  Hist recrateSK1("reconnection rate SK1", 40, 0., 400.);
  Hist recrateSK2("reconnection rate SK2", 40, 0., 400.);

  // Loop over reconnection model and CM energy.
  for (int iRec = 1; iRec < 3; ++iRec) {
    Hist& recrate = (iRec == 1) ? recrateSK1 : recrateSK2;
    for (int iEcm = 14; iEcm < 40; ++iEcm) {
      double eCM = 10. * iEcm - 5.;

      // Quiet generator creation. Shorthand for the event record.
      Pythia pythia ("../share/Pythia8/xmldoc", false);
      Event& event = pythia.event;

      // Incoming beams and energy. Optionally no ISR.
      pythia.readString("Beams:idA = -11");
      pythia.readString("Beams:idB = 11");
      pythia.settings.parm("Beams:eCM", eCM);
      pythia.readString("PDF:lepton = off");

      // Hard process, with hadronic W decays.
      pythia.readString("WeakDoubleBoson:ffbar2WW = on");
      pythia.readString("24:onMode = off");
      pythia.readString("24:onIfAny = 1 2 3 4 5");

      // Switch on CR.
      pythia.readString("ColourReconnection:reconnect = on");
      pythia.settings.mode("ColourReconnection:mode", iRec + 2);
      pythia.readString("ColourReconnection:forceResonance = on");

      // Reduce printout. Switch off hadronization and decays.
      pythia.readString("Print:quiet = on");
      pythia.readString("HadronLevel:Hadronize = off");
      pythia.readString("HadronLevel:Decay = off");

      // If Pythia fails to initialize, exit with error.
      if (!pythia.init()) return 1;

      int nRecon = 0;

      // Begin event loop.
      for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;

        // Check if colour reconnection has occured.
        bool hasRec = false;
        for (int i = 6; i < event.size(); ++i)
          if (event[i].statusAbs() == 79) hasRec = true;
        if (hasRec) ++nRecon;

      // End of loops. Fill reconnection rate.
      }
      recrate.fill( eCM, double(nRecon) / double(nEvent));
    }
  }

  // Print and plot histograms.
  cout << recrateSK1 << recrateSK2;
  HistPlot hpl("plot302");
  hpl.frame( "fig302", "Reconnection rate as a function of CM energy for "
    "e$^+$e$^-$ $\\to$ W$^+$W$^-$", "$E_{\\mathrm{CM}}$ (GeV)", "Probability");
  hpl.add( recrateSK1, "h,blue", "SK I");
  hpl.add( recrateSK2, "h,red", "SK II");
  hpl.plot();

  // Done.
  return 0;
}
