// main171.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: Hidden Valley

// Test of Hidden Valley production in a few different channels.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Some main settings. Many other options exist; see examples further down.

  // Process type: 1 = Fv pair, 2 = Zv portal to qv, 3 = Higgs portal to gv.
  int  procType = 1;

  // Details of simulation: QCD/QED, gammav massive or not, simplifications.
  bool doHVQCD     = true;
  bool brokenQED   = false;
  bool onlyFSR     = false;
  bool onlyHVinFSR = false;

  // Number of events.
  int nEvent = 1000;

  // Begin to set up generator: beams.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");

  // Production process via Fv Fvbar.
  if (procType == 1) {
    pythia.readString("HiddenValley:gg2DvDvbar = on");
    pythia.readString("4900001:m0 = 1000.");

  // Production process via Zv.
  } else if (procType == 2) {
    pythia.readString("HiddenValley:ffbar2Zv = on");
    pythia.readString("4900023:m0   = 1000.");
    pythia.readString("4900023:mMin =  500.");
    pythia.readString("4900023:mMax = 1500.");
    pythia.readString("4900023:onMode = off");
    pythia.readString("4900023:onIfAny = 4900101 4900102 4900103");

  // Production process via H, and decay to gv gv or gammav gammav.
  } else {
    pythia.readString("HiggsSM:all = on");
    pythia.readString("25:onMode = 0");
    if (doHVQCD)
      pythia.readString("25:addChannel = 1 0.1 100 4900021 4900021");
    else
      pythia.readString("25:addChannel = 1 0.1 100 4900022 4900022");
  }

  // Hidden-Valley parton shower.
  pythia.readString("Hiddenvalley:FSR = on");
  if (doHVQCD) pythia.readString("Hiddenvalley:alphaOrder = 1");
  pythia.readString("Hiddenvalley:Lambda = 4.");
  pythia.readString("HiddenValley:pTminFSR = 6.");

  // Fragmentation process of qv to HV-mesons.
  if (doHVQCD) pythia.readString("HiddenValley:Ngauge = 3");
  else         pythia.readString("HiddenValley:Ngauge = 1");
  pythia.readString("HiddenValley:nFlav = 3");
  pythia.readString("HiddenValley:fragment = on");

  // Mass of a Hidden-Valley photon.
  if (brokenQED) pythia.readString("4900022:m0 = 1.");
  else           pythia.readString("4900022:mayDecay = off");

  // Switch off unwanted parts for HV-only simulation.
  if (onlyFSR) {
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("PartonLevel:MPI = off");
    pythia.readString("HadronLevel:all = off");
  }
  if (onlyHVinFSR) {
    pythia.readString("TimeShower:QCDshower = off");
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("TimeShower:QEDshowerByL = off");
    pythia.readString("TimeShower:QEDshowerByGamma = off");
  }

  // Restrict output. Initialize.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.init();

  // Book histograms.
  Hist nGluonv( "number of HV gluons",  100, -0.5, 99.5);
  Hist nGammav( "number of HV gammas",  100, -0.5, 99.5);
  Hist nHadronv("number of HV hadrons", 100, -0.5, 99.5);

   // Begin event loop. Generate event. Extra HV colour output.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent == 0 && (pythia.process.hasHVcols() || event.hasHVcols())) {
      pythia.process.listHVcols();
      event.listHVcols();
    }

    // Number of "final" gv, gammav and hadronv in current event.
    int nGluonvNow  = 0;
    int nGammavNow  = 0;
    int nHadronvNow = 0;
    for (int i = 0; i < event.size(); ++i) {
      int idNow = event[i].idAbs();
      int idDau = event[ event[i].daughter1() ].idAbs();
      if      (idNow == 4900021 && idDau != 4900021) ++nGluonvNow;
      else if (idNow == 4900022 && idDau != 4900022) ++nGammavNow;
      else if (idNow >  4900110) ++nHadronvNow;
    }
    nGluonv.fill( nGluonvNow);
    nGammav.fill( nGammavNow);
    nHadronv.fill( nHadronvNow);

  // End of event loop. Print statistics and histograms.
  }
  pythia.stat();
  cout << nGluonv << nGammav << nHadronv;

  // Done.
  return 0;
}
