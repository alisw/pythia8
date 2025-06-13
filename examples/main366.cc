// main366.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: charm; fixed target;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// This program shows how to calculate charm hadron asymmetries in
// fixed-target pi- p collisions, and compare results with data.
// The data are at different energies, but model is only for one.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Switch on the QCDCR model of Jesper&Peter.
  bool qcdCR = false;

  // Number of events. Number of cases. Ecm from fixed-target beam energy.
  int nEvent   = 100000;
  int nc       = 2;
  double eBeam = 500.;
  double eCM   = sqrt(2.* eBeam * 0.938);

  // Histograms.
  double sigmaT[4];
  Hist DnegxF[4], DposxF[4], DasymxF[4], cbarxF[4], cqrkxF[4], casymxF[4],
    corigxF[4];
  for (int ic = 0; ic < 3; ++ic) {
    DnegxF[ic].book( "xF(D-) "                 , 64, -0.8, 0.8);
    DposxF[ic].book( "xF(D+) "                 , 64, -0.8, 0.8);
    DasymxF[ic].book("(D- - D+)/(D- + D+)(xF) ", 64, -0.8, 0.8);
    cbarxF[ic].book( "xF(cbar)"                , 64, -0.8, 0.8);
    cqrkxF[ic].book( "xF(c)"                   , 64, -0.8, 0.8);
    casymxF[ic].book("(cbar - c) / (cbar + c)" , 64, -0.8, 0.8);
    corigxF[ic].book("cbar + c" ,                64, -0.8, 0.8);
  }
  Hist zeroxF( "", 64, -0.8, 0.8);

  // Loop over cases.
  for (int ic = 0; ic < nc; ++ic) {

    // Generator. Process selection.
    Pythia pythia;
    pythia.readString("Beams:idA = -211");

    // Fixed target energy 500 GeV converted.
    pythia.settings.parm("Beams:eCM", eCM);
    if (ic == 0)      pythia.readString("HardQCD:qqbar2ccbar = on");
    else if (ic == 1) pythia.readString("HardQCD:gg2ccbar = on");
    else              pythia.readString("HardQCD:hardccbar = on");
    pythia.readString("Next:numberShowProcess = 1");
    pythia.readString("Next:numberShowEvent = 1");
    pythia.readString("Next:numberCount = 100000");

    // Alternative setup with Jesper&Peter's colour reconnection.
    if (qcdCR) {
      pythia.readString("ColourReconnection:mode = 1");
      pythia.readString("BeamRemnants:remnantMode = 1");
    }

    // Simplify event generation.
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("PartonLevel:FSR = off");
    pythia.readString("111:mayDecay = off");
    pythia.readString("BeamRemnants:primordialKTsoft = 0.5");
    pythia.readString("BeamRemnants:primordialKThard = 0.5");

    // PDF selection.
    //pythia.readString("PDF:pSet = 1");
    //pythia.readString("PDF:piSet = 1");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Event shorthand.
    Event& event = pythia.event;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Loop over particles in the event. Get basic properties.
      for (int i = 1; i < event.size(); ++i) {
        int idNow    = event[i].id();
        double xFnow = 2. * event[i].pz() / eCM;

        // Original c quarks before selection.
        if (abs(idNow) == 4 && event[event[i].daughter1()].statusAbs() > 80)
          corigxF[ic].fill( xFnow);

        // Analyze D+- hadrons. Trace back to mother c quark.
        if (abs(idNow) == 411) {
          int cType = (idNow == 411) ? 4 : -4;
          if (cType == -4) DnegxF[ic].fill( xFnow);
          else             DposxF[ic].fill( xFnow);
          int iMot1 = event[i].mother1();
          int iMot2 = max(iMot1, event[i].mother2());
          if (event[iMot1].statusAbs() > 80) {
            int iUp = iMot1;
            iMot1   = event[iUp].mother1();
            iMot2   = max(iMot1, event[iUp].mother2());
          }
          int nWork = 0;
          for (int iM = iMot1; iM <= iMot2; ++iM)
          if (event[iM].id() == cType
            || (event[iM].id()/1000 == cType && (event[iM].id()/10)%10 == 0)) {
            double xFmot = 2. * event[iM].pz() / eCM;
            if (cType == -4) cbarxF[ic].fill( xFmot);
            else             cqrkxF[ic].fill( xFmot);
            ++nWork;
          }
          if (nWork != 1) {
            cout << " Missed:  i = " << i << " with iMot1 = " << iMot1
                 << " and iMot2 = " << iMot2 << endl;
            event.list();
          }
        }

      // End of loop over particles in event.
      }

    // End of event loop. Statistics. Histogram handling.
    }
    pythia.stat();
    sigmaT[ic] = pythia.info.sigmaGen();
    DnegxF[ic] *= 10. / nEvent;
    DposxF[ic] *= 10. / nEvent;
    DasymxF[ic] = (DnegxF[ic] - DposxF[ic]) / (DnegxF[ic] + DposxF[ic]);
    cbarxF[ic] *= 10. / nEvent;
    cqrkxF[ic] *= 10. / nEvent;
    casymxF[ic] = (cbarxF[ic] - cqrkxF[ic]) / (cbarxF[ic] + cqrkxF[ic]);
    corigxF[ic] *= 10. / nEvent;
    cout << DnegxF[ic] << DposxF[ic] << DasymxF[ic]
         << cbarxF[ic] << cqrkxF[ic] << casymxF[ic] << corigxF[ic];

  // End of collision type loop. Combined statistics.
  }
  sigmaT[2] = sigmaT[0] + sigmaT[1];
  DnegxF[2] = ( sigmaT[0] * DnegxF[0] + sigmaT[1] * DnegxF[1] ) / sigmaT[2];
  DposxF[2] = ( sigmaT[0] * DposxF[0] + sigmaT[1] * DposxF[1] ) / sigmaT[2];
  DasymxF[2] = (DnegxF[2] - DposxF[2]) / (DnegxF[2] + DposxF[2]);
  cbarxF[2] = ( sigmaT[0] * cbarxF[0] + sigmaT[1] * cbarxF[1] ) / sigmaT[2];
  cqrkxF[2] = ( sigmaT[0] * cqrkxF[0] + sigmaT[1] * cqrkxF[1] ) / sigmaT[2];
  casymxF[2] = (cbarxF[2] - cqrkxF[2]) / (cbarxF[2] + cqrkxF[2]);
  corigxF[2] = ( sigmaT[0] * corigxF[0] + sigmaT[1] * corigxF[1] ) / sigmaT[2];

  // Plot histograms.
  HistPlot hpl("plot366");
  hpl.frame("fig366", "$c + \\overline{c}$ production as function of $x_F$",
    "$x_F$", "(1/$N$) d$N(c + \\overline{c})$/d$x_F$", 8.0, 5.4);
  hpl.add( corigxF[0], "-,red",  "$q\\overline{q} \\to c\\overline{c}$");
  hpl.add( corigxF[1], "-,blue", "$gg \\to c\\overline{c}$");
  hpl.add( corigxF[2], "-,black", "combined");
  hpl.plot();
  hpl.frame("", "$D^-$ production as function of $x_F$",
    "$x_F$", "(1/$N$) d$N(D+)$/d$x_F$", 8.0, 5.4);
  hpl.add( DnegxF[0], "-,red",   "$q\\overline{q} \\to c\\overline{c}$");
  hpl.add( DnegxF[1], "-,blue",  "$gg \\to c\\overline{c}$");
  hpl.add( DnegxF[2], "-,black",  "combined");
  hpl.add( cbarxF[0], "--,orange",
    "$\\overline{c}$ in $q\\overline{q} \\to c\\overline{c}$");
  hpl.add( cbarxF[1], "--,cyan",
    "$\\overline{c}$ in $gg \\to c\\overline{c}$");
  hpl.add( cbarxF[2], "--,grey", "combined");
  hpl.plot();
  hpl.frame("", "$D^+$ production as function of $x_F$",
    "$x_F$", "(1/$N$) d$N(D^+)$/d$x_F$", 8.0, 5.4);
  hpl.add( DposxF[0], "-,red",   "$q\\overline{q} \\to c\\overline{c}$");
  hpl.add( DposxF[1], "-,blue",  "$gg \\to c\\overline{c}$");
  hpl.add( DposxF[2], "-,black",  "combined");
  hpl.add( cqrkxF[0], "-,orange",
    "$c$ in $q\\overline{q} \\to c\\overline{c}$");
  hpl.add( cqrkxF[1], "--,cyan",
    "$c$ in $gg \\to c\\overline{c}$");
  hpl.add( cqrkxF[2], "--,grey", "combined");
  hpl.plot();
  hpl.frame("", "Asymmetry $A(x_F) = (D^- - D^+)/(D^- + D^+)$",
    "$x_F$", "A($x_F$)", 8.0, 5.4);
  hpl.add( DasymxF[0], "-,red",
    "$q\\overline{q} \\to c\\overline{c}$ @ 500 GeV");
  hpl.add( DasymxF[1], "-,blue", "$gg \\to c\\overline{c}$ @ 500 GeV");
  hpl.add( DasymxF[2], "-,orange", "combined");
  hpl.addFile( "main366E769asym.tab", "x,black", "WA82 @ 340 GeV", "y");
  hpl.addFile( "main366E791asym.tab", "*,black", "E769 @ 250 GeV", "y");
  hpl.addFile( "main366WA82asym.tab", "o,black", "E791 @ 500 GeV", "y");
  hpl.add( zeroxF, "-,black");
  hpl.plot( -0.8, 0.8, -1., 1.);

  // Done.
  return 0;
}
