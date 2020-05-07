// main89.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Stefan Prestel <stefan.prestel@thep.lu.se>.

// Keywords: matching; merging; leading order; NLO; powheg; madgraph; aMC@NLO;
// CKKW-L; UMEPS; NL3; UNLOPS; FxFx; MLM; userhooks; LHE file; hepmc;

// This program illustrates how to do run PYTHIA with LHEF input, allowing a
// sample-by-sample generation of
// a) Non-matched/non-merged events
// b) MLM jet-matched events (kT-MLM, shower-kT, FxFx)
// c) CKKW-L and UMEPS-merged events
// d) UNLOPS NLO merged events
// see the respective sections in the online manual for details.
//
// An example command is
//     ./main89 main89ckkwl.cmnd hepmcout89.dat
// where main89.cmnd supplies the commands and hepmcout89.dat is the
// output file. This example requires HepMC 3.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include <unistd.h>

// Include UserHooks for Jet Matching.
#include "Pythia8Plugins/CombineMatchingInput.h"
// Include UserHooks for randomly choosing between integrated and
// non-integrated treatment for unitarised merging.
#include "Pythia8Plugins/aMCatNLOHooks.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging.

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 3) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Output file for HepMC events" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // New setting to allow processing of multiple input LHEFs.
  pythia.settings.addMode("LHEFInputs:nSubruns",0,true,false,0,100);

  // Input parameters:
  pythia.readFile(argv[1],0);

  // Interface for conversion from Pythia8::Event to HepMC one.
  HepMC3::Pythia8ToHepMC3 toHepMC;
  // Specify file where HepMC events will be stored.
  HepMC3::WriterAscii ascii_io(argv[2]);
  // Switch off warnings for parton-level events.
  toHepMC.set_print_inconsistency(false);
  toHepMC.set_free_parton_warnings(false);
  // Do not store cross section information, as this will be done manually.
  toHepMC.set_store_pdf(false);
  toHepMC.set_store_proc(false);
  toHepMC.set_store_xsec(false);

  // Check if jet matching should be applied.
  bool doMatch   = pythia.settings.flag("JetMatching:merge");

  // Check if internal merging should be applied.
  bool doMerge   = !(pythia.settings.word("Merging:Process").compare("void")
    == 0);

  // Currently, only one scheme at a time is allowed.
  if (doMatch && doMerge) {
    cerr << " Jet matching and merging cannot be used simultaneously.\n"
         << " Program stopped.";
  }

  // Get number of subruns.
  int nMerge = pythia.mode("LHEFInputs:nSubruns");
  bool doMatchMerge = true;
  if (nMerge == 0) { nMerge = 1; doMatchMerge = false; }

  // Number of events. Negative numbers mean all events in the LHEF will be
  // used.
  long nEvent = pythia.settings.mode("Main:numberOfEvents");
  if (nEvent < 1) nEvent = 1000;

  // For jet matching, initialise the respective user hooks code.
  //shared_ptr<UserHooks> matching;

  // Allow to set the number of addtional partons dynamically.
  shared_ptr<amcnlo_unitarised_interface> setting;
  if ( doMerge ) {
    // Store merging scheme.
    int scheme = ( pythia.settings.flag("Merging:doUMEPSTree")
                || pythia.settings.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( pythia.settings.flag("Merging:doUNLOPSTree")
                || pythia.settings.flag("Merging:doUNLOPSSubt")
                || pythia.settings.flag("Merging:doUNLOPSLoop")
                || pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
    setting = make_shared<amcnlo_unitarised_interface>(scheme);
    pythia.setUserHooksPtr(setting);
  }

  // For jet matching, initialise the respective user hooks code.
  CombineMatchingInput combined;
  if (doMatch) combined.setHook(pythia);

  vector<double> xss;

  // Allow usage also for non-matched configuration.
  if(!doMatchMerge) {
    // Loop over subruns with varying number of jets.
    for (int iMerge = 0; iMerge < nMerge; ++iMerge) {
      // Read in file for current subrun and initialize.
      pythia.readFile(argv[1], iMerge);
      // Initialise.
      pythia.init();
      // Start generation loop
      while( pythia.info.nSelected() < nEvent ){
        // Generate next event
        if( !pythia.next() ) {
          if ( pythia.info.atEndOfFile() ) break;
          else continue;
        }
      } // end loop over events to generate.
      // print cross section, errors
      pythia.stat();
      xss.push_back(pythia.info.sigmaGen());
    }
  }

  // Cross section and error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  // Allow abort of run if many errors.
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  int  iAbort  = 0;
  bool doAbort = false;

  cout << endl << endl << endl;
  cout << "Start generating events" << endl;

  // Loop over subruns with varying number of jets.
  for (int iMerge = 0; iMerge < nMerge; ++iMerge) {

    double sigmaSample = 0., errorSample = 0.;

    // Read in name of LHE file for current subrun and initialize.
    pythia.readFile(argv[1], iMerge);

    // If the process string is "guess", temporarily set it to something safe
    // for initialization.
    bool doGuess = pythia.settings.word("Merging:process") == "guess";
    if (doMerge && doGuess) pythia.settings.word("Merging:process","pp>e+e-");
    // Initialise.
    pythia.init();
    // Reset the process string to "guess" if necessary.
    if (doGuess) pythia.settings.word("Merging:process","guess");

    // Get the inclusive x-section by summing over all process x-sections.
    double xs = 0.;
    for (int i=0; i < pythia.info.nProcessesLHEF(); ++i)
      xs += pythia.info.sigmaLHEF(i);

    if (!doMatchMerge) xs = xss[iMerge];

    // Start generation loop
    while( pythia.info.nSelected() < nEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if ( pythia.info.atEndOfFile() ) break;
        else if (++iAbort > nAbort) {doAbort = true; break;}
        else continue;
      }

      // Get event weight(s).
      double evtweight         = pythia.info.weight();
      // Additional PDF/alphaS weight for internal merging.
      if (doMerge) evtweight  *= pythia.info.mergingWeightNLO()
      // Additional weight due to random choice of reclustered/non-reclustered
      // treatment. Also contains additional sign for subtractive samples.
                                *setting->getNormFactor();

      // Do not print zero-weight events.
      if ( evtweight == 0. ) continue;
      // Construct new empty HepMC event.
      HepMC3::GenEvent hepmcevt;

      // Work with weighted (LHA strategy=-4) events.
      double normhepmc = 1.;
      if (abs(pythia.info.lhaStrategy()) == 4)
        normhepmc = 1. / double(1e9*nEvent);
      // Work with unweighted events.
      else
        normhepmc = xs / double(1e9*nEvent);

      // Set event weight
      hepmcevt.weights().push_back(evtweight*normhepmc);
      // Fill HepMC event
      toHepMC.fill_next_event( pythia, &hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal  += evtweight*normhepmc;
      sigmaSample += evtweight*normhepmc;
      errorTotal  += pow2(evtweight*normhepmc);
      errorSample += pow2(evtweight*normhepmc);
      // Report cross section to hepmc.
      shared_ptr<HepMC3::GenCrossSection> xsec;
      xsec = make_shared<HepMC3::GenCrossSection>();
      xsec->set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt.set_cross_section( xsec );
      // Write the HepMC event to file.
      ascii_io.write_event(hepmcevt);

    } // end loop over events to generate.
    if (doAbort) break;

    // print cross section, errors
    pythia.stat();

    cout << endl << " Contribution of sample " << iMerge
         << " to the inclusive cross section : "
         << scientific << setprecision(8)
         << sigmaSample << "  +-  " << sqrt(errorSample)  << endl;

  }

  cout << endl << endl << endl;
  if (doAbort)
    cout << " Run was not completed owing to too many aborted events" << endl;
  else
    cout << "Inclusive cross section: " << scientific << setprecision(8)
         << sigmaTotal << "  +-  " << sqrt(errorTotal) << " mb " << endl;
  cout << endl << endl << endl;

  // Done
  return 0;

}
