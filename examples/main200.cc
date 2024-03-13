// main200.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Peter Skands <peter.skands@monash.edu>

// Keywords: Vincia; Dire; electron-positron; histograms

// Simple example of the VINCIA (or DIRE) shower model(s), on Z decays at
// LEP I, with some basic event shapes, spectra, and multiplicity counts.
// Also useful as a basic test of the respective final-state showers.
// Also: how to book and get statistical uncertainties for Pythia histograms.

// Include Pythia8 header(s) and namespace.
#include "Pythia8/Pythia.h"
using namespace Pythia8;

// Main Program
int main() {

  //************************************************************************
  // Define Pythia 8 generator

  Pythia pythia;

  // Read user settings from file
  pythia.readFile("main200.cmnd");

  //************************************************************************

  // Shorthands
  Event& event = pythia.event;
  Settings& settings = pythia.settings;

  // Extract settings to be used in the main program.
  int    nEvent     = settings.mode("Main:numberOfEvents");
  int    nAbort     = settings.mode("Main:timesAllowErrors");
  bool   vinciaOn   = settings.mode("PartonShowers:model") == 2;
  bool   helicityOn = vinciaOn && settings.flag("Vincia:helicityShower");
  int    iEventPri  = -1;
  double eCM        = settings.parm("Beams:eCM");

  //************************************************************************

  // Initialize
  if(!pythia.init()) { return EXIT_FAILURE; }

  //************************************************************************

  // Define a few PYTHIA utilities
  Thrust Thr(1);
  Sphericity SphLin(1, 2);

  //************************************************************************

  // Define PYTHIA histograms. The two last (boolean) arguments are optional
  // and specify whether to use a logarithmic x axis and wether to compute
  // statistical uncertainties when filling the histograms, respectively.

  // Book Histograms. Last arg = true => include statistical uncertainties.
  Hist histNQuarks("nQuarkPairs (Not counting Born)",
    50, -1.0, 99.0, false, true);
  Hist histNPartons("nPartons",  50, -0.5, 49.5, false, true);
  Hist histNCharged("nCharged",  50, -1.0, 99.0, false, true);
  Hist histNBaryons("nBaryons",  25, -1.0, 49.0, false, true);
  Hist histNGamma("nPhotons", 50, -0.5, 49.5, false, true);
  // Examples of explicit log(x) scales (second-to-last arg = false).
  Hist histX1("Ln(x) for 1st branching (QCD)", 100, -5.0, 0.0, false, true);
  Hist histX1Gamma("Ln(x) for 1st branching (QED)", 100, -5.0, 0.0,
    false, true);
  Hist histXUD("Ln(x) for up and down quarks", 25, -5.0, 0.0, false, true);
  Hist histXStrange("Ln(x) for strange quarks", 25, -5.0, 0.0, false, true);
  Hist histXCharm("Ln(x) for charm quarks", 25, -5.0, 0.0, false, true);
  Hist histXBottom("Ln(x) for bottom quarks", 25, -5.0, 0.0, false, true);
  Hist histMUD("Invariant mass of u-ubar and d-dbar pairs",
    100, 0., 50., false, true);
  Hist histMSS("Invariant mass of s-sbar pairs", 100, 0., 50.,
    false, true);
  Hist histMCC("Invariant mass of c-cbar pairs", 100, 0., 50.,
    false, true);
  Hist histMBB("Invariant mass of b-bbar pairs", 100, 0., 50.,
    false, true);

  // Thrust, C, and D parameters.
  // Example of implicit log(x) scales (second-to-last arg = true).
  int nBinsShapes = 100;
  Hist histT("1-T", nBinsShapes, 0.001, 0.5, true, true);
  Hist histC("C", nBinsShapes, 0.001, 1.0, true, true);
  Hist histD("D", nBinsShapes, 0.001, 1.0, true, true);
  double wHistT = nBinsShapes/0.5;
  double wHistCD = nBinsShapes/1.0;


  //************************************************************************

  // EVENT GENERATION LOOP.
  // Generation, event-by-event printout, analysis, and histogramming.

  // Counter for negative-weight events
  double weight = 1.0;
  double sumWeights = 0.0;

  // Begin event loop
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    bool aborted = !pythia.next();
    if(aborted){
      event.list();
      if (++iAbort < nAbort){
        continue;
      }
      cout << " Event generation aborted prematurely, owing to error!\n";
      cout << " Event number was : " << iEvent << endl;
      break;
    }

    // Check for weights
    weight = pythia.info.weight();
    sumWeights += weight;

    // Print event with helicities
    if (iEvent == 0 || iEvent == iEventPri)
      if (helicityOn) event.list(helicityOn);

    // Count FS charged hadrons, partons, and quarks
    int nCh = 0;
    int nBaryons = 0;
    int nPartons = 0;
    int nQuarks  = 0;
    int nGam     = 0;
    // Vectors to store indices of final copies of quark pairs
    // that do not come from the Born.
    vector<int> iQuarkNew;
    for (int i = 1;i<event.size();i++) {
      // Plot x distributions of first branching
      if (i >= 8 && i <= 10) {
        if (event[i].id() == 21)
          histX1.fill(log(2*event[i].pAbs()/eCM), weight);
        else if (event[i].id() == 22)
          histX1Gamma.fill(log(2*event[i].pAbs()/eCM), weight);
      }
      // Find last parton-level partons
      int iDau1 = event[i].daughter1();
      if (iDau1 == 0 || abs(event[iDau1].status()) > 80) {
        // Count up partons and quarks + antiquarks.
        if (event[i].isQuark() || event[i].isGluon()) nPartons++;
        if (event[i].id() == 22) nGam++;
        if (event[i].isQuark()) {
          nQuarks++;
          // If this was not a (recoiling copy of a) Born quark.
          if (event[i].iTopCopyId() >= 8) iQuarkNew.push_back(i);
        }
      }
      if (event[i].isHadron() && event[i].isFinal()) {
        int idAbs = abs(event[i].id());
        if (idAbs > 1000 && idAbs < 10000) nBaryons++;
        if (event[i].isCharged()) nCh++;
      }
    }

    // Don't include original Born pair in count of quark pairs.
    histNQuarks.fill( 0.5*nQuarks - 1., weight);
    histNPartons.fill( nPartons, weight);
    histNCharged.fill( nCh , weight);
    histNBaryons.fill( nBaryons, weight);
    histNGamma.fill( nGam, weight);

    // Histogram quark x fractions and invariant masses
    for (unsigned int iQ = 0; iQ < iQuarkNew.size(); ++iQ) {
      int i = iQuarkNew[iQ];
      int idAbs = event[i].idAbs();
      double x = 2*event[i].pAbs()/eCM;
      if (idAbs <= 2) histXUD.fill(log(x), weight);
      else if (idAbs == 3) histXStrange.fill(log(x), weight);
      else if (idAbs == 4)  histXCharm.fill(log(x), weight);
      else if (idAbs == 5)  histXBottom.fill(log(x), weight);
      // Inner loop to histogram invariant masses of same-flavour quark pairs.
      for (unsigned int jQ = iQ+1; jQ < iQuarkNew.size(); ++jQ) {
        if ( event[iQuarkNew[jQ]].idAbs() != idAbs) continue;
        double mQQ = sqrt(m2(event[i], event[iQuarkNew[jQ]]));
        if (idAbs <= 2) histMUD.fill( mQQ, weight);
        else if (idAbs == 3) histMSS.fill( mQQ, weight);
        else if (idAbs == 4) histMCC.fill( mQQ, weight);
        else if (idAbs == 5) histMBB.fill( mQQ, weight);
      }
    }

    // Histogram thrust
    Thr.analyze( event );
    histT.fill(1.0-Thr.thrust(), wHistT*weight);

    // Histogram Linear Sphericity values
    if (nPartons >= 2.0) {
      SphLin.analyze( event );
      double evC = 3*(SphLin.eigenValue(1)*SphLin.eigenValue(2)
        + SphLin.eigenValue(2)*SphLin.eigenValue(3)
        + SphLin.eigenValue(3)*SphLin.eigenValue(1));
      double evD = 27*SphLin.eigenValue(1)*SphLin.eigenValue(2)
        *SphLin.eigenValue(3);
      histC.fill(evC, wHistCD*weight);
      histD.fill(evD, wHistCD*weight);
    }

  }

  //************************************************************************

  // POST-RUN FINALIZATION
  // Normalization, Statistics, Output.

  //Normalize histograms to effective number of positive-weight events.
  double normFac = 1.0/sumWeights;
  histT          *= normFac;
  histC          *= normFac;
  histD          *= normFac;
  histNQuarks    *= normFac;
  histNPartons   *= normFac;
  histNCharged   *= normFac;
  histNBaryons   *= normFac;
  histNGamma     *= normFac;
  histXUD        *= normFac;
  histXStrange   *= normFac;
  histXCharm     *= normFac;
  histXBottom    *= normFac;
  histX1         *= normFac;
  histX1Gamma    *= normFac;

  // Print a few histograms.
  cout << histNPartons << endl;
  cout << histNGamma << endl;
  cout << histNCharged << endl;
  cout << histNBaryons << endl;
  cout << histT << endl;
  cout << histC << endl;
  cout << histD << endl;
  cout << histX1 << endl;
  cout << histX1Gamma << endl;
  cout << histXStrange << endl;
  cout << histXCharm << endl;
  cout << histXBottom << endl;
  cout << histMUD << endl;
  cout << histMSS << endl;
  cout << histMCC << endl;
  cout << histMBB << endl;

  // Print out end-of-run information.
  pythia.stat();

  cout << endl;
  cout << fixed;
  cout << " <nPartons>    = " << num2str(histNPartons.getXMean())
       << "  +/- " << num2str(histNPartons.getXMeanErr()) << endl;
  cout << " <nGluonSplit> = " << num2str(histNQuarks.getXMean())
       << "  +/- " << num2str(histNQuarks.getXMeanErr())
       << " Flavours( <ud>, s, c, b ) = "
       << num2str(histXUD.getWeightSum()/4.) << " , "
       << num2str(histXStrange.getWeightSum()/2.) << " , "
       << num2str(histXCharm.getWeightSum()/2.) << " , "
       << num2str(histXBottom.getWeightSum()/2.) <<endl;
  cout << " <nPhotons>    = " << num2str(histNGamma.getXMean())
       << "  +/- " << num2str(histNGamma.getXMeanErr()) << endl;
  cout << " <nCharged>    = " << num2str(histNCharged.getXMean())
       << "  +/- " << num2str(histNCharged.getXMeanErr()) << endl;
  cout << " <nBaryons>    = " << num2str(histNBaryons.getXMean())
       << "  +/- " << num2str(histNBaryons.getXMeanErr()) << endl;
  cout<<endl;

  // Done.
  return 0;
}
