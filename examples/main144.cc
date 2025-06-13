// main144.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// Keywords: analysis; hepmc; command file; command line option; root; rivet;
//           tuning

// Streamlined event generation with possibility to output ROOT files,
// output HepMC files and run RIVET analyses, all by specifying output modes
// in a cmnd file, where also the event generator settings are specified.
// The example is run with command line options, run ./main144 -h to see a
// full list. See ROOT Usage for information about ROOT output, RIVET Usage
// for information about RIVET and HepMC Interface for information about HepMC.

#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/InputParser.h"
#include <chrono>
#ifdef RIVET
#include "Pythia8Plugins/Pythia8Rivet.h"
#endif
#ifdef PY8ROOT
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#endif

// Use the Pythia namespace.
using namespace Pythia8;

//==========================================================================

// Define filling ROOT particle and events if using ROOT.

#ifdef PY8ROOT
#include "main144Dct.h"

// Fill a ROOT particle.
RootParticle::RootParticle(Pythia8::Particle &prt) {
  phi = prt.phi();
  eta = prt.eta();
  y   = prt.y();
  pT  = prt.pT();
  pid = prt.id();
}

// Fill a ROOT event to a TTree.
void RootEvent::fill(const Pythia8::Info &infoIn, vector<RootParticle> &prtsIn,
  TTree *treeIn) {
  weight = infoIn.weight();
  particles = prtsIn;
  treeIn->Fill();
}
#endif

//==========================================================================

// Implement a user supplied UserHooks derived class inside this
// wrapper, which will allow you to give settings that can be supplied
// in command files.

class UserHooksWrapper : public UserHooks {

public:

  // Add the settings you want available in the run card in this method.
  void additionalSettings(Settings* settingsIn) {
    settings = settingsIn;
    settings->addFlag("UserHooks:doMPICut", false);
    settings->addMode("UserHooks:nMPICut", 0, true, false, 0, 0);
  }

  // Check if parton level can be vetoed.
  bool canVetoPartonLevel() final {
    return settings->flag("UserHooks:doMPICut");}

  // Check if parton level should be vetoed.
  bool doVetoPartonLevel(const Event&) final {
    return infoPtr->nMPI() < settings->mode("UserHooks:nMPICut");}

private:

  Settings* settings{};

};

//==========================================================================

// Main example execution.

int main(int argc, char* argv[]) {

  // Parser object for command line input.
  InputParser ip("Run Pythia with cmnd file input, and get Rivet, HepMC or"
    " standard Pythia output.",
    {"./main144 [options]", "./main144 -c main144.cmnd -n 1000 -o myoutput"},
    "Additional options in cmnd file:\n"
    "\tMain:writeLog = on\n\t\tRedirect output to <-o prefix>.log.\n"
    "\tMain:writeHepMC = on \n\t\tWrite HepMC output, requires HepMC linked.\n"
    "\tMain:writeRoot = on \n\t\tWrite a ROOT tree declared in "
    "RootEvent.h, requires ROOT linked.\n"
    "\tMain:runRivet = on \n\t\tRun Rivet analyses, requires Rivet linked.\n"
    "\tMain:rivetAnalyses = {ANALYSIS1,ANALYSIS2,...}\n "
    "\t\tComma separated list of Rivet analyses to run.\n"
    "\t\tAnalysis names can be post-fixed with analysis parameters.\n"
    "\t\tANALYSIS:parm=value:parm2=value2:...\n"
    "\tMain:rivetRunName = STRING \n\t\tAdd an optional run name to "
    "the Rivet analysis.\n"
    "\tMain:rivetIgnoreBeams = on\n\t\tIgnore beams in Rivet. \n"
    "\tMain:rivetDumpPeriod = NUMBER\n\t\tDump Rivet histograms "
    "to file evert NUMBER of events.\n"
    "\tMain:rivetDumpFile = STRING\n\t\t Specify alternative "
    "name for Rivet dump file. Default = OUT.\n");

  // Set up command line options.
  ip.require("c", "User-written command file, can use multiple times.",
    {"-cmnd"});
  ip.add("s", "-1", "Specify seed for the random number generator.",
    {"-seed"});
  ip.add("o", "main144", "Output prefix for log file, Rivet, HepMC, and ROOT.",
    {"-out"});
  ip.add("n", "-1", "Number of events. Overrides the command files.",
    {"-nevents"});
  ip.add("l", "false", "Silence the splash screen.");
  ip.add("t", "false", "Time event generation.", {"-time"});
  ip.add("v", "false", "Print Pythia version number and exit.", {"-version"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Print version number and exit.
  if (ip.get<bool>("v")) {
    cout << "PYTHIA version: " << PYTHIA_VERSION << endl;
    return 0;
  }

  // Get the command files.
  vector<string> cmnds = ip.getVector<string>("c");
  if (cmnds.size() == 0) {
    cout << "Please provide one or more command files with the -c option."
         << endl;
    return 1;
  }

  // Random number seed.
  string seed = ip.get<string>("s");
  // Output filename.
  string out = ip.get<string>("o");
  // Time event generation.
  bool writeTime = ip.get<bool>("t");
  // Command line number of event, overrides the one set in input .cmnd file.
  int nev = ip.get<int>("n");

  // Catch the splash screen in a buffer.
  stringstream splashBuf;
  std::streambuf* sBuf = cout.rdbuf();
  cout.rdbuf(splashBuf.rdbuf());
  // The Pythia object.
  Pythia pythia;
  // Direct cout back.
  cout.rdbuf(sBuf);

  // UserHooks wrapper.
  shared_ptr<UserHooksWrapper> userHooksWrapper =
    make_shared<UserHooksWrapper>();
  userHooksWrapper->additionalSettings(&pythia.settings);
  pythia.setUserHooksPtr(userHooksWrapper);

  // Some extra parameters.
  pythia.settings.addFlag("Main:writeLog", false);
  pythia.settings.addFlag("Main:writeHepMC", false);
  pythia.settings.addFlag("Main:writeRoot", false);
  pythia.settings.addFlag("Main:runRivet", false);
  pythia.settings.addFlag("Main:rivetIgnoreBeams", false);
  pythia.settings.addMode("Main:rivetDumpPeriod", -1, true, false, -1, 0);
  pythia.settings.addWord("Main:rivetDumpFile", "");
  pythia.settings.addWord("Main:rivetRunName", "");
  pythia.settings.addWVec("Main:rivetAnalyses", {});
  pythia.settings.addWVec("Main:rivetPreload", {});

  // Read the command files.
  for (int iCmnd = 0; iCmnd < (int)cmnds.size(); ++iCmnd)
    if (!cmnds[iCmnd].empty()) pythia.readFile(cmnds[iCmnd]);

  // Set seed after reading input.
  if(seed != "-1") {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = "+seed);
  }

  // Read the extra parameters.
  if (nev > -1) pythia.settings.mode("Main:numberOfEvents", nev);
  int nEvent                   = pythia.mode("Main:numberOfEvents");;
  int nError                   = pythia.mode("Main:timesAllowErrors");
  bool writeLog                = pythia.flag("Main:writeLog");
  bool writeHepmc              = pythia.flag("Main:writeHepMC");
  bool writeRoot               = pythia.flag("Main:writeRoot");
  bool runRivet                = pythia.flag("Main:runRivet");
  bool countErrors             = nError > 0;

  // Check if Rivet, HepMC, and ROOT are requested and available.
  bool valid = true;
#ifndef RIVET
  valid = valid && !runRivet && !writeHepmc;
  if (runRivet)
    cout << "Option Main::runRivet = on requires the Rivet library.\n";
  if (writeHepmc)
    cout << "Option Main::writeHepMC = on requires the HepMC library.\n";
#endif
#ifndef PY8ROOT
  valid = valid && !writeRoot;
  if (writeRoot)
    cout << "Option Main::writeRoot = on requires the ROOT library.\n";
#endif
  if (!valid) return 1;

  // Rivet and HepMC initialization.
#ifdef RIVET
  // Initialize HepMC.
  Pythia8ToHepMC hepmc;
  if (writeHepmc) hepmc.setNewFile(out + ".hepmc");

  // Initialize Rivet.
  Pythia8Rivet rivet(pythia, out + ".yoda");
  rivet.ignoreBeams(pythia.flag("Main:rivetIgnoreBeams"));
  rivet.dump(pythia.settings.mode("Main:rivetDumpPeriod"),
    pythia.settings.word("Main:rivetDumpFile"));

  // Load the analyses.
  vector<string> rivetAnalyses = pythia.settings.wvec("Main:rivetAnalyses");
  for (int iAna = 0; iAna < (int)rivetAnalyses.size(); ++iAna)
    rivet.addAnalysis(rivetAnalyses[iAna]);

  // Pre-load the YODA histograms.
  vector<string> rivetPreload = pythia.settings.wvec("Main:rivetPreload");
  for (int iYoda = 0; iYoda < (int)rivetPreload.size(); ++iYoda)
    rivet.addPreload(rivetPreload[iYoda]);

  // Add the run name.
  rivet.addRunName(pythia.settings.word("Main:rivetRunName"));
#endif

  // ROOT initialization.
#ifdef PY8ROOT
  // Create the ROOT TFile and TTree.
  TFile *file;
  TTree *tree;
  RootEvent *evt;
  if (writeRoot) {

    // Open the ROOT file.
    file = TFile::Open((out + ".root").c_str(), "recreate" );
    tree = new TTree("t", "Pythia8 event tree");
    evt  = new RootEvent();

    // Set the TTree branch to the ROOT event.
    tree->Branch("events", &evt);
  }
#endif

  // Logfile initialization.
  ofstream logBuf;
  streambuf *oldCout;
  if (writeLog) {
    oldCout = cout.rdbuf(logBuf.rdbuf());
    logBuf.open(out + ".log");
  }

  // Remove splash screen, if requested.
  ostream cnull(NULL);
  if (ip.get<bool>("l")) cnull << splashBuf.str();
  else cout << splashBuf.str();

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Make a sanity check of initialized Rivet analyses.
#ifdef RIVET
  if (!runRivet && rivetAnalyses.size() > 0 )
    cout << "Rivet analyses are set with Main:rivetAnalyses, "
         << "but Main:runRivet = off.\n";
#endif

  // Loop over events.
  auto startAllEvents = std::chrono::high_resolution_clock::now();
  for ( int iEvent = 0; iEvent < nEvent; ++iEvent ) {
    auto startThisEvent = std::chrono::high_resolution_clock::now();

    // Exit if too many failures.
    if (!pythia.next()) {
      if (countErrors && --nError < 0) {
        pythia.stat();
        cout << " \n *-------  PYTHIA STOPPED!  -----------------------*\n"
             << " | Event generation failed due to too many errors. |\n"
             << " *-------------------------------------------------*\n";
        return 1;
      }
      continue;
    }

    // Calculate the event time.
    auto stopThisEvent = std::chrono::high_resolution_clock::now();
    auto eventTime = std::chrono::duration_cast<std::chrono::milliseconds>
      (stopThisEvent - startThisEvent);
    double tt = eventTime.count();

    // Run the Rivet analyses.
#ifdef RIVET
    if (runRivet) {
      if (writeTime) rivet.addAttribute("EventTime", tt);
      rivet();
    }
    if (writeHepmc) hepmc.writeNextEvent(pythia);
#endif

    // Write to ROOT file output.
#ifdef PY8ROOT
    if (writeRoot) {
      vector<RootParticle> prts;
      for (int iPrt = 0; iPrt < pythia.event.size(); ++iPrt) {
        Particle& prt = pythia.event[iPrt];

        // Any particle cuts can be placed here. Here, only final
        // state particles are kept.
        if (!prt.isFinal()) continue;

        // Push back the ROOT particle.
        prts.push_back(RootParticle(prt));
      }
      // Fill the ROOT event and tree.
      evt->fill(pythia.info, prts, tree);
    }
#endif
  }

  // Finalize.
  pythia.stat();
#ifdef PY8ROOT
  if (writeRoot) {
    tree->Print();
    tree->Write();
    delete file, tree, evt;
  }
#endif

  // Print timing.
  auto stopAllEvents = std::chrono::high_resolution_clock::now();
  auto durationAll = std::chrono::duration_cast<std::chrono::milliseconds>
    (stopAllEvents - startAllEvents);
  if (writeTime) {
    cout << " \n *-------  Generation time  -----------------------*\n"
         << " | Event generation, analysis and writing to files  |\n"
         << " | took: " << double(durationAll.count()) << " ms or "
         << double(durationAll.count())/double(nEvent)
         << " ms per event     |\n"
         << " *-------------------------------------------------*\n";
  }

  // Put cout back in its place.
  if (writeLog) cout.rdbuf(oldCout);
  return 0;

}
