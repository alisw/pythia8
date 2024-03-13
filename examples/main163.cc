// main163.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.m.utheim@jyu.fi>

// Keywords: parallelism

// This is a simple test program to illustrate the usage of PythiaParallel.
// This program illustrates how to perform both event generation and
// analysis in parallel, using your own mutex object.

#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"
#include <mutex>
#include <thread>

using namespace Pythia8;
int main() {

  // Basic settings.
  PythiaParallel pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Main:numberOfEvents = 10000");

  // This tells PythiaParallel to process events asynchronously.
  // If this is set to off, the program will slow down significantly.
  pythia.readString("Parallelism:processAsync = on");

  // Initialize.
  pythia.init();
  Hist mult("charged multiplicity", 100, -0.5, 799.5);

  // This mutual exclusion (mutex) object controls access to histogram.
  mutex histMutex;

  // Generate events.
  pythia.run([&](Pythia* pythiaPtr) {

    // Find number of all final charged particles and fill histogram.
    int nCharged = 0;
    for (int i = 0; i < pythiaPtr->event.size(); ++i)
      if (pythiaPtr->event[i].isFinal() && pythiaPtr->event[i].isCharged())
        ++nCharged;

    // Simulate a slow analysis by delaying for 20 milliseconds.
    std::this_thread::sleep_for(std::chrono::milliseconds(20));

    // Lock mutex. The above part of the analysis can be done in parallel,
    // but two threads must not write to the histogram at the same time.
    // If this line is removed, the output will be wrong.
    std::lock_guard<mutex> lock(histMutex);

    // Fill histogram
    mult.fill( nCharged );

    // The mutex will be released when the lock_guard goes out of scope.
  });

  pythia.stat();
  cout << mult;
  return 0;
}
