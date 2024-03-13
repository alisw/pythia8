// main40.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Philip Ilten <philten@cern.ch>
//          Naomi Cooke <naomi.cooke@cern.ch>
//          Leif Lonnblad <leif.lonnblad@fysik.lu.se>
//          Steve Mrenna <mrenna@fnal.gov>

// Keywords: onia

// This calculates the inclusive branching fractions for the Standard Model
// Higgs into quarkonia using the LETO parton shower.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Single-particle gun. The particle must be a colour singlet.

void fillParticle(int id, Event& event, ParticleData& pdt, Rndm& rndm) {

  // Reset event record to allow for new event.
  event.reset();

  // Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(id);
  double ee = mm;
  double pp = 0.;

  // Angles as input or uniform in solid angle.
  double cThe = 2. * rndm.flat() - 1.;
  double sThe = sqrtpos(1. - cThe * cThe);
  double phi  = 2. * M_PI * rndm.flat();

  // Store the particle in the event record.
  event.append( id, 1, 0, 0, pp * sThe * cos(phi),
    pp * sThe * sin(phi), pp * cThe, ee, mm);

}

//==========================================================================

// Print a table.

void printTable(ParticleData &pdt, string title, vector<int> states,
  vector< pair<string, string> > labels,
  map<string, map<int, int> > &counters, double scale = 1.,
  double enhance = 1.) {
  cout << "----------------------------------------------------------------\n"
       << title + "\n"
       << "----------------------------------------------------------------\n";
  cout << setw(7) << " " << setw(14) << left << " state" << right;
  for (auto label : labels) cout << setw(10) << label.second;
  cout << "\n";
  for (int state : states) {
    cout << setw(7) << state << setw(14) << left
         <<  " (" + pdt.name(state) + ")" << right;
    for (auto label : labels) {
      double scaleNow = label.first == "nFDH" ? scale : enhance*scale;
      if (scale*enhance == 1) cout << fixed << setprecision(0);
      else cout << scientific << setprecision(2);
      cout << setw(10) << counters[label.first][state]/scaleNow;
    }
    cout << "\n";
  }

}

//==========================================================================

int main() {

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Configure Pythia and initialize.
  pythia.readFile("main40.cmnd");
  pythia.readString("ProcessLevel:all = off");
  pythia.init();

  // Map of counters for onium production.
  map<string, map<int, int> > counters = {
    {"nShower", {}}, {"nOctet", {}}, {"nSinglet", {}}, {"nGluon", {}},
    {"nQuark", {}}, {"nFDO", {}}, {"nFDH", {}}, {"nGG", {}}, {"nCC", {}},
    {"nBB", {}}, {"nVV", {}}};

  // Begin of event loop.
  int nEvent = pythia.settings.mode("Main:numberOfEvents");
  int nAcc = 0;
  bool hasPrintedOnium = false;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Set up single particle, with random direction in solid angle.
    fillParticle(25, event, pdt, pythia.rndm);

    // Generate events. Continue if failure.
    if (!pythia.next()) continue;
    else nAcc++;

    // Loop over the particles.
    int nOnia = 0;
    for (int iPrt = 0; iPrt < event.size(); ++iPrt) {

      // Find the onium (only color singlet states).
      Particle *oPrt = &event[iPrt];
      ParticleDataEntry *oPDE = &oPrt->particleDataEntry();
      if (!oPDE->isOnium() || oPDE->isOctetHadron()) continue;

      // Find the mother and ignore onia from hadronization.
      int oTop = oPrt->iTopCopyId();
      if (event[oTop].mother2() != 0) continue;
      Particle *mPrt = &event[event[oTop].mother1()];
      ParticleDataEntry *mPDE = &mPrt->particleDataEntry();

      // Consider first hadron as mother.
      while (mPrt->mother1() > 0 && event[mPrt->mother1()].isHadron()) {
        mPrt = &event[mPrt->mother1()];
        mPDE = &mPrt->particleDataEntry();
      }

      // Fill the counters.
      int oID(oPrt->idAbs());
      if (!mPDE->isHadron()) nOnia++;

      // Fill if octet or singlet.
      if (mPDE->isOctetHadron())  counters["nOctet"][oID]++;
      else if (!mPDE->isHadron()) counters["nSinglet"][oID]++;

      // Fill the production mechanism.
      if      (mPDE->isGluon())  counters["nGluon"][oID]++;
      else if (mPDE->isQuark())  counters["nQuark"][oID]++;
      else if (mPDE->isOnium())  counters["nFDO"][oID]++;
      else if (mPDE->isHadron()) counters["nFDH"][oID]++;

      // Fill the Higgs branching; only consider shower production.
      if (!mPDE->isHadron()) {
        counters["nShower"][oID]++;
        int pID(event[2].idAbs());
        if      (pID == 21)              counters["nGG"][oID]++;
        else if (pID == 4)               counters["nCC"][oID]++;
        else if (pID == 5)               counters["nBB"][oID]++;
        else if (pID == 23 || pID == 24) counters["nVV"][oID]++;

        // Make a printout of first event that contains a shower onium.
        if (!hasPrintedOnium) {
          cout << "\nFirst event with an onium from shower (i = "
               << iPrt << ", id = " << oPrt->id() << "):"
               << " iEvent = " << iEvent << endl;
          event.list();
          hasPrintedOnium = true;
        }

      }
    }
    if (nOnia > 1) cout << "WARNING: more than one onia found\n";

    // End of event loop.
  }

  // Define counters and order of states to print.
  vector<int> states = {
    441, 443, 10441, 20443, 445, 551, 553, 100553, 200553, 10551, 20553, 555};
  vector< pair<string, string> > labels = {
    {"nFDO", "FD-onium"}, {"nFDH", "FD-hadron"}, {"nShower", "shower"},
    {"nOctet", "octet"}, {"nSinglet", "singlet"}, {"nGluon", "g->X1"},
    {"nQuark", "Q->X1"}, {"nGG", "H->gg"}, {"nCC", "H->gg"}, {"nBB", "H->bb"},
    {"nVV", "H->VV"}};

  // Print the tables.
  printTable(pdt, "raw event counts", states, labels, counters);
  printTable(pdt, "brancing fractions", states, labels, counters, nAcc,
    pythia.settings.parm("OniaShower:ldmeFac"));
  return 0;

}
