// main31.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: matching; merging; powheg

// Example how to perform matching with POWHEG-BOX events,
// based on the code found in include/Pythia8Plugins/PowhegHooks.h.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Generator
  Pythia pythia;

  // Load configuration file
  pythia.readFile("main31.cmnd");

  // Read in main settings.
  int nEvent      = pythia.settings.mode("Main:numberOfEvents");
  int nError      = pythia.settings.mode("Main:timesAllowErrors");
  // Read in key POWHEG matching settings.
  int powhegVeto    = pythia.settings.mode("POWHEG:veto");
  int powhegMPIveto = pythia.settings.mode("POWHEG:MPIveto");
  bool loadHooks  = (powhegVeto > 0 || powhegMPIveto > 0);
  // Read in shower settings.
  int showerModel = pythia.settings.mode("PartonShowers:model");

  // Add in user hooks for shower vetoing.
  shared_ptr<PowhegHooks> powhegHooks;
  if (loadHooks) {

    // For POWHEG:veto >= 1, setup to do vetoed power showers.
    //   1) Setting pTmaxMatch = 2 forces PYTHIA's shower to sweep over the
    //      full phase space.
    //   2) Loading the POWHEG hooks will then veto any shower branchings that
    //      are judged (according to the POWHEG settings) to double-count the
    //      POWHEG one.
    if (powhegVeto > 0) {
      if (showerModel == 1 || showerModel == 3) {
        // For PYTHIA's simple shower (and also for Dire), the FSR and ISR
        // shower starting scales are set by the respective pTmaxMatch values.
        pythia.readString("TimeShower:pTmaxMatch = 2");
        pythia.readString("SpaceShower:pTmaxMatch = 2");
        // Use undamped power showers, except for cases for which there could
        // be an interplay with ISR branchings not simulated by the pure
        // shower, like ISR g->tt.
        pythia.readString("SpaceShower:pTdampMatch = 3");
        pythia.readString("TimeShower:pTdampMatch = 0");
      } else if (showerModel == 2) {
        // Vincia has common settings that apply to both ISR and FSR.
        pythia.readString("Vincia:pTmaxMatch = 2");
        // Use undamped power showers, except for cases for which there could
        // be an interplay with ISR branchings not simulated by the pure
        // shower, like g->tt.
        pythia.readString("Vincia:pTdampMatch = 3");
      }
    }

    // For POWHEG:MPIveto >= 1, also set MPI to start at the kinematical limit.
    if (powhegMPIveto > 0) {
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
    }

    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);
  }

  // Initialise and list settings
  pythia.init();

  // Counters for number of ISR/FSR emissions vetoed
  unsigned long int nISRveto = 0, nFSRveto = 0;

  // Begin event loop; generate until nEvent events are processed
  // or end of LHEF file
  int iEvent = 0, iError = 0;
  while (true) {

    // Generate the next event
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop
      if (pythia.info.atEndOfFile()) break;

      // Otherwise count event failure and continue/exit as necessary
      cout << "Warning: event " << iEvent << " failed" << endl;
      if (++iError == nError) {
        cout << "Error: too many event failures... exiting" << endl;
        break;
      }

      continue;
    }

    /*
     * Process dependent checks and analysis may be inserted here
     */

    // Update ISR/FSR veto counters
    if (loadHooks) {
      nISRveto += powhegHooks->getNISRveto();
      nFSRveto += powhegHooks->getNFSRveto();
    }

    // If nEvent is set, check and exit loop if necessary
    ++iEvent;
    if (nEvent != 0 && iEvent == nEvent) break;

  } // End of event loop.

  // Statistics, histograms and veto information
  pythia.stat();
  cout << "Number of ISR emissions vetoed: " << nISRveto << endl;
  cout << "Number of FSR emissions vetoed: " << nFSRveto << endl;
  cout << endl;

  // Done.
  return 0;
}
