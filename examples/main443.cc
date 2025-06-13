// main443.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: thermal model;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Compare pT spectra and particle composition in thermal vs ordinary
// hadronization models, and also ina "straw man" alternative ("mT2").
// Note that the models are out-of-the-box, without any parameter tuning.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Histograms in line printer mode or as pdf file.
  int nEvent = 2000;
  bool histAsPDF = true;

  // Histogram multiplicities and pT spectra.
  Hist nCh[3], pTpi[3], pTK[3], pTp[3];
  for (int iModel = 0; iModel < 3; ++iModel) {
    nCh[iModel].book("charged multiplicity", 100, 1., 401.);
    pTpi[iModel].book("pT for pi+-", 100, 0., 5.);
    pTK[iModel].book("pT for K+-", 100, 0., 5.);
    pTp[iModel].book("pT for p/pbar", 100, 0., 5.);
  }

  // Tabulate particle composition.
  string nameHad[18] = { "pi+-", "pi0", "K+-", "K0S,L", "eta,eta'", "rho+-",
    "rho0,omega", "K*+-0", "phi0", "p(bar)", "n(bar)", "Lambda(bar)",
    "Sigma(bar)", "Xi(bar)", "Delta(bar)", "Sigma*(bar)", "Xi*(bar)",
    "Omega(bar)"};
  int rates[3][18] = {{0}};

  // Loop over normal and thermal generation.
  for (int iModel = 0; iModel < 3; ++iModel) {

    // Generator. Common process selection. Reduce printout. pi0 stable.
    Pythia pythia;
    pythia.readString("Beams:eCM = 8000.");
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("111:mayDecay = off");

    // Model-specific setup. Initialization.
    if (iModel == 1) pythia.readString("StringPT:thermalModel = on");
    if (iModel == 2) pythia.readString("StringPT:mT2suppression = on");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Event shorthand.
    Event& event = pythia.event;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Find charged multiplicity and fill histograms.
      int nCharged = 0;
      for (int i = 0; i < event.size(); ++i)
        if (event[i].isFinal() && event[i].isCharged()) ++nCharged;
      nCh[iModel].fill( nCharged );

      // Find pT spectra of some hadron species.
      for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
        int idAbs = event[i].idAbs();
        double pT = event[i].pT();
        if (idAbs ==  211) pTpi[iModel].fill( pT);
        if (idAbs ==  321) pTK[iModel].fill( pT);
        if (idAbs == 2212) pTp[iModel].fill( pT);
      }

      // Tabulate hadron composition, including unstable intermediaries,
      // but excluding charm and bottom.
      for (int i = 0; i < event.size(); ++i) {
        int idAbs = event[i].idAbs();
        if      (idAbs == 211) ++rates[iModel][0];
        else if (idAbs == 111) ++rates[iModel][1];
        else if (idAbs == 321) ++rates[iModel][2];
        else if (idAbs == 130 || idAbs == 310) ++rates[iModel][3];
        else if (idAbs == 221 || idAbs == 331) ++rates[iModel][4];
        else if (idAbs == 213) ++rates[iModel][5];
        else if (idAbs == 113 || idAbs == 223) ++rates[iModel][6];
        else if (idAbs == 313 || idAbs == 323) ++rates[iModel][7];
        else if (idAbs == 333) ++rates[iModel][8];
        else if (idAbs == 2212) ++rates[iModel][9];
        else if (idAbs == 2112) ++rates[iModel][10];
        else if (idAbs == 3122) ++rates[iModel][11];
        else if (idAbs == 3112 || idAbs == 3212 || idAbs == 3222)
          ++rates[iModel][12];
        else if (idAbs == 3312 || idAbs == 3322) ++rates[iModel][13];
        else if (idAbs == 1114 || idAbs == 2114 || idAbs == 2214
          || idAbs == 2224) ++rates[iModel][14];
        else if (idAbs == 3114 || idAbs == 3214 || idAbs == 3224)
          ++rates[iModel][15];
        else if (idAbs == 3314 || idAbs == 3324) ++rates[iModel][16];
        else if (idAbs == 3334) ++rates[iModel][17];
     }

    // End of event loop. Statistics. Normalize histograms.
    }
    pythia.stat();
    nCh[iModel]  *= 0.5 / double(nEvent);
    pTpi[iModel] *= 20. / double(nEvent);
    pTK[iModel]  *= 20. / double(nEvent);
    pTp[iModel]  *= 20. / double(nEvent);

  // End of model loop.
  }

  // Print historams.
  if (!histAsPDF) {
    for (int iModel = 0; iModel < 3; ++iModel) cout << nCh[iModel];
    for (int iModel = 0; iModel < 3; ++iModel) cout << pTpi[iModel];
    for (int iModel = 0; iModel < 3; ++iModel) cout << pTK[iModel];
    for (int iModel = 0; iModel < 3; ++iModel) cout << pTp[iModel];

  // Alternatively plot histograms.
  } else {
    HistPlot hpl("plot443");
    hpl.frame("fig443", "charged multiplicity",
      "$n_{\\mathrm{charged}}$", "Probability", 8.0, 5.4);
    hpl.add( nCh[0], "-,black",  "default");
    hpl.add( nCh[1], "-,red",    "thermal");
    hpl.add( nCh[2], "-,blue",   "$m_{\\perp}^2$-suppressed");
    hpl.plot();
    hpl.frame("", "$\\pi^{\\pm}$ transverse momentum spectrum",
      "$p_{\\perp}$", "$\\mathrm{d}n/\\mathrm{d}p_{\\perp}$", 8.0, 5.4);
    hpl.add( pTpi[0], "-,black",  "default");
    hpl.add( pTpi[1], "-,red",    "thermal");
    hpl.add( pTpi[2], "-,blue",   "$m_{\\perp}^2$-suppressed");
    hpl.plot();
    hpl.frame("", "K$^{\\pm}$ transverse momentum spectrum",
      "$p_{\\perp}$", "$\\mathrm{d}n/\\mathrm{d}p_{\\perp}$", 8.0, 5.4);
    hpl.add( pTK[0], "-,black",  "default");
    hpl.add( pTK[1], "-,red",    "thermal");
    hpl.add( pTK[2], "-,blue",   "$m_{\\perp}^2$-suppressed");
    hpl.plot();
    hpl.frame("", "p,$\\overline{\\mathrm{p}}$ transverse momentum spectrum",
      "$p_{\\perp}$", "$\\mathrm{d}n/\\mathrm{d}p_{\\perp}$", 8.0, 5.4);
    hpl.add( pTp[0], "-,black",  "default");
    hpl.add( pTp[1], "-,red",    "thermal");
    hpl.add( pTp[2], "-,blue",   "$m_{\\perp}^2$-suppressed");
    hpl.plot();
  }

  // Print table.
  double norm = 1. / double(nEvent);
  cout << "\n\n Particle composition per event, including unstable"
       << "\n    Particle     default     thermal     mT2-suppressed" << endl
       << fixed << setprecision(3);
  for (int i = 0; i < 18; ++i) cout << setw(12) << nameHad[i]
       << setw(12) << norm * rates[0][i] << setw(12) << norm * rates[1][i]
       << setw(12) << norm * rates[2][i] << endl;

  // Done.
  return 0;
}
