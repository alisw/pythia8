// PythiaParallel.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Marius Utheim, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef PYTHIA_PARALLEL_H
#define PYTHIA_PARALLEL_H

#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Class for doing Pythia runs in parallel.

class PythiaParallel {

public:

  // Constructor.
  PythiaParallel(string xmlDir = "../share/Pythia8/xmldoc",
    bool printBanner = true);

  // Read in one update for a setting or particle data from a single line.
  bool readString(string setting, bool warn = true) {
    return pythiaHelper.readString(setting, warn); }

  // Read in updates for settings or particle data from user-defined file.
  bool readFile(string fileName, bool warn = true,
    int subrun = SUBRUNDEFAULT);
  bool readFile(string fileName, int subrun) {
    return readFile(fileName, true, subrun); }
  bool readFile(istream& is = cin, bool warn = true,
    int subrun = SUBRUNDEFAULT);
  bool readFile(istream& is, int subrun) {
    return readFile(is, true, subrun);}

  // Initialize all Pythia objects.
  bool init();
  bool init(function<bool(Pythia*)> additionalSetup);

  // Perform the specified action for each Pythia instance.
  void foreach(function<void(Pythia*)> action);

  // Perform the specified action for each instance in parallel.
  void foreachAsync(function<void(Pythia*)> action);

  // Write final statistics, combining errors from each Pythia instance.
  void stat() { pythiaHelper.stat(); }

  // Generate events in parallel.
  vector<long> run(long nEvents, function<void(Pythia*)> callback);
  vector<long> run(function<void(Pythia*)> callback) {
    return run(settings.mode("Main:numberOfEvents"), callback); }

  // Pythia object used for loading data.
  Pythia pythiaHelper;

  // Weighted average of the generated cross section for each Pythia instance.
  double sigmaGen() const { return sigmaGenSave; }

  // Sum of weights from all Pythia instances.
  double weightSum() const { return weightSumSave; }

  // The settings that will be used to initialize Pythia instances.
  Settings& settings;

  // The particle database that will be used to initialize Pythia instances.
  ParticleData& particleData;

private:

  // Object used for logging.
  Logger& logger;

  // Flag if initialized.
  bool isInit = false;

  // Sum of weights and weighted cross section.
  double weightSumSave, sigmaGenSave;

  // Configuration flags.
  int numThreads;
  bool processAsync;
  bool balanceLoad;
  bool doNext;

  // Internal Pythia objects.
  vector<unique_ptr<Pythia> > pythiaObjects;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PythiaParallel_H
