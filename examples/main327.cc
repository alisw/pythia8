// main327.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: forward physics; hadron production;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Study forward p production with diffraction + nondiffractive.
// Compare default, ditto without popcorn for remnant diquark, and
// the QCDCR model.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events, beam energy and more.
  int nEvent = 5000;
  double eBeam = 3500.;
  double etaForw = 10.76;

  // Include single diffraction on the elastic side of the event or not.
  bool includeSDel = false;

  // Histograms.
  int nSides[4] = { 0 };
  Hist yp[4], xFp[4], pTp[4], yn[4], xFn[4], pTn[4], xFpDel[4], xFnDel[4],
    xFpEta[4], xFnEta[4], pTpEta[4], pTnEta[4], ypi[4], xFpi[4], ygam[4],
    xFgam[4], xFdi[4], xFpn[4], xFpi2[4], xFdiq[4];
  for (int i = 0; i < 4; ++i) {
    yp[i].book("y_proton", 100, 0., 10.);
    xFp[i].book("xF_proton", 95, 0.05, 1.);
    pTp[i].book("pT_proton, xF > 0.4", 100, 0., 2.);
    yn[i].book("y_neutron", 100, 0., 10.);
    xFn[i].book("xF_neutron", 95, 0.05, 1.);
    pTn[i].book("pT_neutron, xF > 0.4", 100, 0., 2.);
    xFpDel[i].book("xF_proton Delta", 100, 0., 1.);
    xFnDel[i].book("xF_neutron Delta", 100, 0., 1.);
    xFpEta[i].book("xF_proton, eta > 10.76", 100, 0., 1.);
    xFnEta[i].book("xF_neutron eta > 10.76", 100, 0., 1.);
    pTpEta[i].book("pT_proton, eta > 10.76", 100, 0., 2.);
    pTnEta[i].book("pT_neutron eta > 10.76", 100, 0., 2.);
    ypi[i].book("y_pi+-", 100, 0., 10.);
    xFpi[i].book("xF_pi+-", 95, 0.05, 1.);
    ygam[i].book("y_gamma", 100, 0., 10.);
    xFgam[i].book("xF_gamma", 95, 0.05, 1.);
    xFdi[i].book("xF_diquark", 95, 0.05, 1.);
    xFpn[i].book("xF_proton+neutron", 90, 0.1, 1.);
    xFpi2[i].book("xF_pi+-", 90, 0.1, 1.);
    xFdiq[i].book("xF_diquark", 100, 0., 1.);
  }

  // Loop over cases
  for (int colType = 0; colType < 3; ++colType) {

    // Generator. Process selection.
    Pythia pythia;
    pythia.settings.parm("Beams:eCM", 2. * eBeam);
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("Next:numberShowProcess = 2");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString(" Next:numberCount = 100000");

    // Options with no remnant popcorn and QCDCR.
    if (colType == 1) {
      pythia.readString("BeamRemnants:dampPopcorn = 0.0");
    } else if (colType == 2) {
      pythia.readString("ColourReconnection:mode = 1");
      pythia.readString("BeamRemnants:remnantMode = 1");
    }

    // Map of remnant constituents.
    map<int, int> remnantTypes;

    // Reuse MPI initialization.
    pythia.readString("MultipartonInteractions:reuseInit = 3");
    pythia.readString("MultipartonInteractions:initFile = main327.mpi");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Event shorthand.
    Event& event = pythia.event;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Classify whether both sides should be bookkept.
      int procCode = pythia.info.code();
      bool doPos   = includeSDel || (procCode != 104);
      bool doNeg   = includeSDel || (procCode != 103);
      if (doPos) ++nSides[colType];
      if (doNeg) ++nSides[colType];

      // Loop over particles in the event. Get basic properties.
      for (int i = 1; i < event.size(); ++i) {
        bool isPos    = (event[i].pz() > 0.);
        bool isNeg    = (event[i].pz() < 0.);
        bool doFill   = (doPos && isPos) || (doNeg && isNeg);
        if (!doFill) continue;
        double yAbs   = abs(event[i].y());
        double etaAbs = abs(event[i].eta());
        double xFAbs  = abs(event[i].pz()) / eBeam;
        double pTAbs  = event[i].pT();

        // Proton spectra in y, xF, pT. Also with eta cut.
        if (event[i].isFinal() && event[i].id() == 2212) {
          yp[colType].fill( yAbs );
          xFp[colType].fill( xFAbs );
          xFpn[colType].fill( xFAbs );
          if (xFAbs > 0.4) pTp[colType].fill( pTAbs );
          xFpDel[colType].fill( xFAbs );
          if (etaAbs > etaForw) {
            xFpEta[colType].fill( xFAbs );
            pTpEta[colType].fill( pTAbs );
          }
        // Antiproton spectrum subtraction in xF.
        } else if (event[i].isFinal() && event[i].id() == -2212) {
          xFpDel[colType].fill( xFAbs, -1.);

        // Neutron spectra in y, xF, pT. Also with eta cut.
        } else if (event[i].isFinal() && event[i].id() == 2112) {
          yn[colType].fill( yAbs );
          xFn[colType].fill( xFAbs );
          xFpn[colType].fill( xFAbs );
          if (xFAbs > 0.4) pTn[colType].fill( pTAbs );
          xFnDel[colType].fill( xFAbs );
          if (etaAbs > etaForw) {
            xFnEta[colType].fill( xFAbs );
            pTnEta[colType].fill( pTAbs );
          }
        // Antineutron spectrum subtraction in xF.
        } else if (event[i].isFinal() && event[i].id() == -2112) {
          xFnDel[colType].fill( xFAbs, -1.);

        // Remnant diquark spectra in xF.
        } else if (colType == 0 && event[i].isDiquark()
          && event[i].status() == -63) {
            if (procCode == 101) xFdi[0].fill( xFAbs );
            else                 xFdi[1].fill( xFAbs );

        // pi+- spectrum in y and xF.
        } else if (event[i].isFinal() && event[i].idAbs() == 211) {
          ypi[colType].fill( yAbs );
          xFpi[colType].fill( xFAbs );
          xFpi2[colType].fill( xFAbs );

        // Photon spectrum in y and xF.
        } else if (event[i].isFinal() && event[i].id() == 22) {
             ygam[colType].fill( event[i].y() );
             xFgam[colType].fill( event[i].pz() / eBeam);
        }

      // End of loop over particles in event.
      }

      // Fill map with remnant composition and histogram with diquark xF.
      for (int i = 1; i < event.size(); ++i) if (event[i].statusAbs() == 63) {
        ++remnantTypes[ event[i].id()];
        if (event[i].id() > 1000)
          xFdiq[colType].fill( abs(event[i].pz()) / eBeam );
      }

    // End of event loop. Statistics.
    }
    pythia.stat();

    // Print remnant composition.
    for (auto iter = remnantTypes.begin(); iter != remnantTypes.end(); ++iter)
      cout << " id = " << setw(6) << iter->first << " count = " << setw(6)
      << iter->second << endl;

  // End of collision type loop.
  }

  // Normalize histograms.
  for (int i = 0; i < 4; ++i) {
    yp[i]     *=  10. / max( 1, nSides[i]);
    xFp[i]    *= 100. / max( 1, nSides[i]);
    pTp[i]    *=  50. / max( 1, nSides[i]);
    yn[i]     *=  10. / max( 1, nSides[i]);
    xFn[i]    *= 100. / max( 1, nSides[i]);
    pTn[i]    *=  50. / max( 1, nSides[i]);
    xFdi[i]   *= 100. / max( 1, nSides[i]);
    xFpDel[i] *= 100. / max( 1, nSides[i]);
    xFnDel[i] *= 100. / max( 1, nSides[i]);
    xFpEta[i] *= 100. / max( 1, nSides[i]);
    xFnEta[i] *= 100. / max( 1, nSides[i]);
    pTpEta[i] *=  50. / max( 1, nSides[i]);
    pTnEta[i] *=  50. / max( 1, nSides[i]);
    ypi[i]    *=  10. / max( 1, nSides[i]);
    xFpi[i]   *= 100. / max( 1, nSides[i]);
    ygam[i]   *=  10. / max( 1, nSides[i]);
    xFgam[i]  *= 100. / max( 1, nSides[i]);
    xFpn[i]   *= 100. / max( 1, nSides[i]);
    xFpi2[i]  *= 100. / max( 1, nSides[i]);
    xFdiq[i]  *= 100. / max( 1, nSides[i]);
  }

  // Plot histograms.
  HistPlot hpl("plot327");
  hpl.frame("fig327", "Nucleon rapidity distribution",
    "$y$", "d$n_{\\mathrm{p}}$/d$y$", 8.0, 5.4);
  hpl.add( yp[0], "-,red", "p inelastic old default");
  hpl.add( yn[0], "--,red", "n inelastic old default");
  hpl.add( yp[1], "-,blue", "p inelastic no popcorn");
  hpl.add( yn[1], "--,blue", "n inelastic no popcorn");
  hpl.add( yp[2], "-,green", "p inelastic QCDCR");
  hpl.add( yn[2], "--,green", "n inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "Nucleon Feynman-$x$ distribution", "$x_{\\mathrm{F}}$",
    "d$n_{\\mathrm{p}}$/d$x_{\\mathrm{F}}$", 8.0, 5.4);
  hpl.add( xFp[0], "-,red", "p inelastic old default");
  hpl.add( xFn[0], "--,red", "n inelastic old default");
  hpl.add( xFp[1], "-,blue", "p inelastic no popcorn");
  hpl.add( xFn[1], "--,blue", "n inelastic no popcorn");
  hpl.add( xFp[2], "-,green", "p inelastic QCDCR");
  hpl.add( xFn[2], "--,green", "n inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "Nucleon - antinucleon Feynman-$x$ distribution",
    "$x_{\\mathrm{F}}$", "d$n_{\\mathrm{p}}$/d$x_{\\mathrm{F}}$", 8.0, 5.4);
  hpl.add( xFpDel[0], "-,red", "p inelastic old default");
  hpl.add( xFnDel[0], "--,red", "n inelastic old default");
  hpl.add( xFpDel[1], "-,blue", "p inelastic no popcorn");
  hpl.add( xFnDel[1], "--,blue", "n inelastic no popcorn");
  hpl.add( xFpDel[2], "-,green", "p inelastic QCDCR");
  hpl.add( xFnDel[2], "--,green", "n inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "Nucleon transverse momentum distribution, $x_F > 0.4$",
    "$p_{\\perp}$", "d$n_{\\mathrm{p}}$/d$p_{\\perp}$", 8.0, 5.4);
  hpl.add( pTp[0], "-,red", "p inelastic old default");
  hpl.add( pTn[0], "--,red", "n inelastic old default");
  hpl.add( pTp[1], "-,blue", "p inelastic no popcorn");
  hpl.add( pTn[1], "--,blue", "n inelastic no popcorn");
  hpl.add( pTp[2], "-,green", "p inelastic QCDCR");
  hpl.add( pTn[2], "--,green", "n inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "Nucleon Feynman-$x$ distribution, $\\eta > 10.76$",
    "$x_{\\mathrm{F}}$", "d$n_{\\mathrm{p}}$/d$x_{\\mathrm{F}}$", 8.0, 5.4);
  hpl.add( xFpEta[0], "-,red", "p inelastic old default");
  hpl.add( xFnEta[0], "--,red", "n inelastic old default");
  hpl.add( xFpEta[1], "-,blue", "p inelastic no popcorn");
  hpl.add( xFnEta[1], "--,blue", "n inelastic no popcorn");
  hpl.add( xFpEta[2], "-,green", "p inelastic QCDCR");
  hpl.add( xFnEta[2], "--,green", "n inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "Nucleon transverse momentum distribution, $\\eta > 10.76$",
    "$p_{\\perp}$", "d$n_{\\mathrm{p}}$/d$p_{\\perp}$", 8.0, 5.4);
  hpl.add( pTpEta[0], "-,red", "p inelastic old default");
  hpl.add( pTnEta[0], "--,red", "n inelastic old default");
  hpl.add( pTpEta[1], "-,blue", "p inelastic no popcorn");
  hpl.add( pTnEta[1], "--,blue", "n inelastic no popcorn");
  hpl.add( pTpEta[2], "-,green", "p inelastic QCDCR");
  hpl.add( pTnEta[2], "--,green", "n inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "$\\pi^{\\pm}$ and $\\gamma$ rapidity distribution",
    "$y$", "d$n$/d$y$", 8.0, 5.4);
  hpl.add( ypi[0], "-,red", "pi, inelastic old default");
  hpl.add( ypi[1], "-,blue", "pi, inelastic no popcorn");
  hpl.add( ypi[2], "-,green", "pi, inelastic QCDCR");
  hpl.add( ygam[0], "--,red", "gamma, inelastic old default");
  hpl.add( ygam[1], "--,blue", "gamma, inelastic no popcorn");
  hpl.add( ygam[2], "--,green", "gamma, inelastic QCDCR");
  hpl.plot();
  hpl.frame("", "$\\pi^{\\pm}$ and $\\gamma$ Feynman-$x$ distribution",
    "$x_{\\mathrm{F}}$", "d$n$/d$x_{\\mathrm{F}}$", 8.0, 5.4);
  hpl.add( xFpi[0], "-,red", "pi, inelastic old default");
  hpl.add( xFpi[1], "-,blue", "pi, inelastic no popcorn");
  hpl.add( xFpi[2], "-,green", "pi, inelastic QCDCR");
  hpl.add( xFgam[0], "--,red", "gamma, inelastic old default");
  hpl.add( xFgam[1], "--,blue", "gamma, inelastic no popcorn");
  hpl.add( xFgam[2], "--,green", "gamma, inelastic QCDCR");
  hpl.plot(0.05, 1.0, 1e-5, 1e2, true);
  hpl.frame("", "Diquark Feynman-$x$ distribution", "$x_{\\mathrm{F}}$",
    "d$n_{\\mathrm{p}}$/d$x_{\\mathrm{F}}$", 8.0, 5.4);
  hpl.add( xFdi[0], "-,red", "nondiffractive");
  hpl.add( xFdi[1], "-,blue", "diffractive");
  hpl.plot();
  hpl.frame("", "Nucleon Feynman-$x$ distribution",
    "$x_{\\mathrm{F}}$", "d$n$/d$x_{\\mathrm{F}}$", 6.4, 4.8);
  hpl.add( xFpn[0], "-,red", "default");
  hpl.add( xFpn[1], "-,blue", "no popcorn");
  hpl.add( xFpn[2], "-,green", "QCDCR");
  hpl.plot(0.1, 1.0, 0., 3.0);
  hpl.frame("", "$\\pi^{\\pm}$ Feynman-$x$ distribution",
    "$x_{\\mathrm{F}}$", "d$n$/d$x_{\\mathrm{F}}$", 6.4, 4.8);
  hpl.add( xFpi[0], "-,red", "default");
  hpl.add( xFpi[1], "-,blue", "no popcorn");
  hpl.add( xFpi[2], "-,green", "QCDCR");
  hpl.plot(0.1, 1.0, 1e-3, 20., true);
  hpl.frame("", "diquark Feynman-$x$ distribution",
    "$x_{\\mathrm{F}}$", "d$n$/d$x_{\\mathrm{F}}$", 8.0, 5.4);
  hpl.add( xFdiq[0], "-,red", "default");
  hpl.add( xFdiq[1], "-,blue", "no popcorn");
  hpl.add( xFdiq[2], "-,green", "QCDCR");
  hpl.plot();

  // Done.
  return 0;
}
