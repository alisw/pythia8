// PythiaParallel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Marius Utheim, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// PythiaParallel class.

#include "Pythia8/PythiaParallel.h"

namespace Pythia8 {

//==========================================================================

// The PythiaParallel class.

//--------------------------------------------------------------------------

// Contructor.

PythiaParallel::PythiaParallel(string xmlDir, bool printBanner)
  : pythiaHelper(xmlDir, printBanner), settings(pythiaHelper.settings),
    particleData(pythiaHelper.particleData), logger(pythiaHelper.logger)
{ }

//--------------------------------------------------------------------------

// Read in updates for settings or particle data from user-defined file.

bool PythiaParallel::readFile(string fileName, bool warn, int subrun) {
  ifstream is(fileName);
  if (!is.good()) {
    logger.ERROR_MSG("did not find file", fileName);
    return false;
  }
  // Hand over real work to next method.
  return readFile( is, warn, subrun);
}

bool PythiaParallel::readFile(istream& is, bool warn, int subrun) {
  if (isInit) {
    logger.ERROR_MSG("cannot change further settings after constructing");
    return false;
  }
  return pythiaHelper.readFile(is, warn, subrun);
}

//--------------------------------------------------------------------------

// Initialize all Pythia objects.

bool PythiaParallel::init() {
  return init(function<bool(Pythia*)>());
}

bool PythiaParallel::init(function<bool(Pythia*)> customInit) {

  // Initialize error printing.
  logger.init(settings);

  // Read settings.
  int hardwareThreads = thread::hardware_concurrency();
  numThreads = settings.mode("Parallelism:numThreads");
  if (numThreads == 1)
    logger.WARNING_MSG("running on single thread");
  else if (numThreads == 0) {
    if (hardwareThreads == 0) {
      logger.ABORT_MSG(
        "cannot get hardware_concurrency, numThreads must be set manually");
      return false;
    }
    numThreads = hardwareThreads;
    settings.mode("Parallelism:numThreads", hardwareThreads);
    logger.INFO_MSG("detected number of hardware threads",
      to_string(hardwareThreads));
  }
  else if (numThreads > hardwareThreads) {
    logger.WARNING_MSG(
      "requested numThreads is larger than hardware_concurrency",
      to_string(hardwareThreads));
  }
  processAsync = settings.flag("Parallelism:processAsync");
  balanceLoad  = settings.flag("Parallelism:balanceLoad");
  doNext       = settings.flag("Parallelism:doNext");

  if (!doNext && !processAsync) {
    logger.WARNING_MSG(
      "setting both doNext and processAsync to off prevents parallelism");
  }

  // Set seeds.
  vector<int> seeds = settings.mvec("Parallelism:seeds");
  if (seeds.size() == 0) {
    seeds = vector<int>(numThreads);
    int seed0 = settings.mode("Random:seed");
    if (seed0 < 0)
      seed0 = 19780503;
    else if (seed0 == 0)
      seed0 = int(time(0));
    for (int i = 0; i < numThreads; ++i)
      seeds[i] = seed0 + i;
    settings.mvec("Parallelism:seeds", seeds);
  }

  // Create instances in parallel.
  pythiaObjects = vector<unique_ptr<Pythia>>(numThreads);

  vector<thread> initThreads;
  bool initSuccess = true;

  for (int iPythia = 0; iPythia < numThreads; iPythia += 1) {
    initThreads.emplace_back([=, &seeds, &initSuccess]() {
      Pythia* pythiaPtr = new Pythia(settings, particleData, false);
      pythiaObjects[iPythia] = unique_ptr<Pythia>(pythiaPtr);
      pythiaObjects[iPythia]->settings.flag("Print:quiet", true);
      pythiaObjects[iPythia]->settings.flag("Random:setSeed", true);
      pythiaObjects[iPythia]->settings.mode("Random:seed", seeds[iPythia]);
      pythiaObjects[iPythia]->settings.mode("Parallelism:index", iPythia);

      if (customInit && !customInit(pythiaObjects[iPythia].get()))
        initSuccess = false;
      if (!pythiaObjects[iPythia]->init())
        initSuccess = false;
    });
  }

  // Wait for all initialization threads to finish.
  for (int i = 0; i < numThreads; ++i)
    initThreads[i].join();

  // Set initialization.
  if (!initSuccess) {
    logger.ABORT_MSG("failed to initialize all Pythia objects");
    return false;
  }
  isInit = true;

  // Print warning message and return.
  logger.WARNING_MSG(
    "experimental feature, please send feedback to authors@pythia.org");
  return true;

}

//--------------------------------------------------------------------------

// Run Pythia objects.

vector<long> PythiaParallel::run(long nEvents,
  function<void(Pythia* pythiaPtr)> callback) {

  if (!isInit) {
    logger.ABORT_MSG("not initialized");
    return vector<long>();
  }

  if (nEvents < numThreads)
    logger.WARNING_MSG("more threads than events have been specified");
  int numThreadsNow = nEvents > numThreads ? numThreads : int(nEvents);
  long nShowCount = settings.mode("Next:numberCount");

  mutex callbackMutex;
  vector<long> eventsPerThread(numThreadsNow);
  atomic<long> nStartedEvents(0);
  atomic<long> nFinishedEvents(0);
  vector<thread> threads;

  // Define the thread main that will run for each Pythia object.
  auto threadMain = [&, this, callback](Pythia* pythiaPtr, int iPythia) {

    // If load is balanced, we need the number of events to run on this thread.
    long nLocalEvents = nEvents / numThreadsNow;
    if (iPythia < nEvents - (nLocalEvents * numThreadsNow))
      nLocalEvents += 1;

    // Run the Pythia object.
    while (true) {

      // Check if we're done, depending on whether load should be balanced.
      if (balanceLoad) {
        if (nLocalEvents == 0) break;
        nLocalEvents -= 1;
      }
      else if (nStartedEvents++ >= nEvents) break;

      // Generate the event.
      bool success = !doNext || pythiaPtr->next();

      // Increment counter for number of generated events.
      // Note the use of printf for thread safety.
      eventsPerThread[iPythia] += 1;
      long generatedEventsNow = ++nFinishedEvents;
      if ( nShowCount > 0 && generatedEventsNow % nShowCount == 0
        && generatedEventsNow < nEvents)
        printf("\n PythiaParallel::run(): %ld events have been generated\n",
          generatedEventsNow);

      // Pass the generated event to the callback.
      if (success) {
        if (processAsync) {
          callback(pythiaPtr);
        } else {
          // Lock access to the callback.
          const std::lock_guard<mutex> lock(callbackMutex);
          callback(pythiaPtr);
        }
      }
    }
  }; // end thread main

  // Start all threads.
  for (int iPythia = 0; iPythia < numThreadsNow; ++iPythia)
    threads.emplace_back(threadMain, pythiaObjects[iPythia].get(), iPythia);

  // Zero the counters.
  weightSumSave = 0.;
  sigmaGenSave = 0.;

  // Wait for each thread to finish.
  for (int iPythia = 0; iPythia < numThreadsNow; ++iPythia) {
    threads[iPythia].join();
    logger.errorCombine(pythiaObjects[iPythia]->logger);

    double weightSumNow = pythiaObjects[iPythia]->info.weightSum();
    weightSumSave += weightSumNow;
    sigmaGenSave  += weightSumNow * pythiaObjects[iPythia]->info.sigmaGen();
  }

  // Set generated cross section and return.
  sigmaGenSave /= weightSumSave;
  return eventsPerThread;

}

//--------------------------------------------------------------------------

// Perform the specified action for each Pythia instance.

void PythiaParallel::foreach(function<void(Pythia*)> action) {

  if (!isInit) {
    logger.ERROR_MSG("not initialized");
    return;
  }

  // Perform action in serial.
  for (auto& pythia : pythiaObjects) action(pythia.get());
}

//--------------------------------------------------------------------------

// Perform the specified action for each instance in parallel.

void PythiaParallel::foreachAsync(function<void(Pythia*)> action) {

  if (!isInit) {
    logger.ERROR_MSG("not initialized");
    return;
  }

  // Perform action in parallel.
  vector<thread> threads;
  for (auto& pythiaPtr : pythiaObjects)
    threads.emplace_back(action, pythiaPtr.get());
  for (thread& threadNow : threads)
    threadNow.join();

}

//==========================================================================

}  // end namespace Pythia8
