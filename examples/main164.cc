// main164.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Stefan Prestel, Christian T. Preuss

// Contact: Christian T. Preuss <christian.preuss@uni-goettingen.de>

// Keywords: matching; merging; leading order; NLO; powheg; madgraph; aMC@NLO;
//           CKKW-L; UMEPS; NL3; UNLOPS; FxFx; MLM;
//           userhooks; LHE file; HDF5 file; LHEH5; hepmc; rivet

// This program illustrates how to do run PYTHIA with LHEF input, allowing a
// sample-by-sample generation of
// a) Non-matched/non-merged events
// b) NLO matched events (with MadGraph5_aMC@NLO or POWHEG-BOX)
// c) MLM jet-matched events (kT-MLM, shower-kT, FxFx)
// d) CKKW-L and UMEPS LO merged events
// e) UNLOPS NLO merged events
// see the respective sections in the online manual for details.
//
// An example command is
//     ./main164 -c main164ckkwl.cmnd
// where main164ckkwl.cmnd supplies the commands.
// This example requires HepMC2 or HepMC3 and optionally RIVET.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"
#if defined(HEPMC3)
#include "Pythia8Plugins/HepMC3.h"
#elif defined(HEPMC2)
#include "Pythia8Plugins/HepMC2.h"
#endif
#ifdef RIVET
#include "Pythia8Plugins/Pythia8Rivet.h"
#endif
#ifdef HDF5
#include "Pythia8Plugins/LHAHDF5v2.h"
#endif

// Include UserHooks for POWHEG vetos.
#include "Pythia8Plugins/PowhegHooks.h"
// Include UserHooks for Jet Matching.
#include "Pythia8Plugins/CombineMatchingInput.h"
// Include UserHooks for randomly choosing between integrated and
// non-integrated treatment for unitarised merging.
#include "Pythia8Plugins/aMCatNLOHooks.h"

using namespace Pythia8;

//==========================================================================

// General example program for matching and merging in PYTHIA.

int main(int argc, char** argv){

  // Set up command line options.
  InputParser ip("Illustrates how to do matching and merging.",
    {"./main164 -c main164ckkwl.cmnd",
        "./main164 -c main164amcatnlo.cmnd",
        "./main164 -c main164fxfx.cmnd",
        "./main164 -c main164mlm.cmnd",
        "./main164 -c main164umeps.cmnd",
        "./main164 -c main164mess.cmnd",
        "./main164 -c main164unlops.cmnd",
        "./main164 -c main164dis.cmnd",
        "./main164 -c main164powheg.cmnd",});
  ip.require("c", "Use this user-written command file.", {"-cmnd"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Input file.
  string cmndFile = ip.get<string>("c");

  // Generator.
  Pythia pythia;

  // New settings for HepMC interface.
  pythia.settings.addFlag("Main:HepMC", false);
  pythia.settings.addWord("HepMC:output", "events.hepmc");

  // New settings for Rivet interface.
  // Inactive if Rivet not used.
  pythia.settings.addFlag("Main:Rivet", true);
  pythia.settings.addWord("Rivet:output", "analysis.yoda");
  pythia.settings.addWVec("Rivet:analyses", vector<string>());

  // Input parameters:
  pythia.readFile(cmndFile, 0);

  // Optionally HepMC interface.
#if defined(HEPMC2) || defined(HEPMC3)
  const bool doHepMC = pythia.flag("Main:HepMC");
  Pythia8ToHepMC* toHepMCPtr = nullptr;
  if (doHepMC) {
    toHepMCPtr = new Pythia8ToHepMC(pythia.word("HepMC:output"));
    // Switch off warnings for parton-level events.
    toHepMCPtr->set_print_inconsistency(false);
    toHepMCPtr->set_free_parton_warnings(false);
    // Do not store the following information.
    toHepMCPtr->set_store_pdf(false);
    toHepMCPtr->set_store_proc(false);
  }
#endif

  // Optionally Rivet interface if Pythia is linked to Rivet.
#ifdef RIVET
  const bool doRivet = pythia.flag("Main:Rivet");
  Pythia8Rivet rivet(pythia, pythia.word("Rivet:output"));
  vector<string> analyses = pythia.settings.wvec("Rivet:analyses");
  for (const string& analysis : analyses) rivet.addAnalysis(analysis);
#endif

  // Save which shower we are using.
  int showerModel = pythia.mode("PartonShowers:model");

  // Check if PowhegHooks should be used for NLO matching.
  int  pwhgVetoMode     = pythia.mode("POWHEG:veto");
  int  pwhgVetoModeMPI  = pythia.mode("POWHEG:MPIveto");
  bool doPowhegMatching = (pwhgVetoMode>0 || pwhgVetoModeMPI>0);

  // Check if jet matching should be applied.
  bool doJetMatching = pythia.flag("JetMatching:merge");

  // Check if internal merging should be applied.
  bool doMerging = !(pythia.word("Merging:Process").compare("void")==0);

  // Currently, only one scheme at a time is allowed.
  if (doPowhegMatching + doJetMatching + doMerging > 1) {
    cerr << " Error in " << argv[0]
         << ": matching and merging schemes cannot be combined" << endl;
  }

  // Set UserHooks for POWHEG vetos.
  shared_ptr<PowhegHooks> powhegHooks;
  int nVetoISR = 0, nVetoFSR = 0;
  if (doPowhegMatching) {
    // Set showers to start at the kinematical limit.
    if (pwhgVetoMode > 0) {
      if (showerModel == 1 || showerModel == 3) {
        pythia.readString("SpaceShower:pTmaxMatch = 2");
        pythia.readString("TimeShower:pTmaxMatch = 2");
      } else if (showerModel == 2) {
        pythia.readString("Vincia:tune = 0");
        pythia.readString("Vincia:pTmaxMatch = 2");
      }
    }
    // Set MPI to start at the kinematical limit.
    if (pwhgVetoModeMPI > 0)
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
    // Load POWHEG hooks.
    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);
  }

  // Set UserHooks for jet matching.
  CombineMatchingInput jetMatchingHook;
  if (doJetMatching) jetMatchingHook.setHook(pythia);

  // Set UserHooks for unitarised merging schemes.
  shared_ptr<amcnlo_unitarised_interface> mergingHooks;
  if (doMerging) {
    // Store merging scheme.
    int scheme = ( pythia.flag("Merging:doUMEPSTree")
                || pythia.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( pythia.flag("Merging:doUNLOPSTree")
                || pythia.flag("Merging:doUNLOPSSubt")
                || pythia.flag("Merging:doUNLOPSLoop")
                || pythia.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
    // Load merging hooks.
    mergingHooks = make_shared<amcnlo_unitarised_interface>(scheme);
    pythia.setUserHooksPtr(mergingHooks);
  }

  // Get number of subruns and information about external events.
  int nMerge = pythia.mode("Main:numberOfSubruns");
  if (nMerge == 0) nMerge = 1;
  bool useLHA  = (pythia.mode("Beams:frameType") >= 4);
#ifdef HDF5
  bool useHDF5 = (pythia.mode("Beams:frameType") == 5);
#endif

  // Allow abort of run if many errors.
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  int  iAbort  = 0;
  bool doAbort = false;

  // Loop over subruns with varying number of jets.
  for (int iMerge = 0; iMerge<nMerge; ++iMerge) {

    // Read in name of LHE file for current subrun and initialize.
    pythia.readFile(cmndFile, iMerge);

    // Set number of events.
    long nEvent = pythia.mode("Main:numberOfEvents");

    // Set up LHAupPtr for this run when using HDF5 event files.
#ifdef HDF5
    if (useHDF5) {
      HighFive::File file(pythia.word("Beams:LHEF"), HighFive::File::ReadOnly);
      size_t readSize = size_t(nEvent);
      size_t eventOffset = 0;
      shared_ptr<LHAupH5v2> lhaUpPtr =
        make_shared<LHAupH5v2>(&file, eventOffset, readSize, true);
      pythia.setLHAupPtr(lhaUpPtr);
    }
#endif

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Get the inclusive cross section by
    // summing over all process cross sections.
    double xs = 0.;
    if (useLHA) {
      for (int i=0; i<pythia.info.nProcessesLHEF(); ++i)
        xs += pythia.info.sigmaLHEF(i);
    }

    // Start generation loop.
    while (pythia.info.nSelected() < nEvent) {

      // Generate next event.
      if (!pythia.next()) {
        if (pythia.info.atEndOfFile()) break;
        else if (++iAbort > nAbort) {doAbort = true; break;}
        else continue;
      }

      // For internal events, get current cross section.
      if (!useLHA) xs = pythia.info.sigmaGen();

      // For POWHEG matching, count vetos.
      if (doPowhegMatching) {
        nVetoISR += powhegHooks->getNISRveto();
        nVetoFSR += powhegHooks->getNFSRveto();
      }

      // Get event weight.
      // Includes additional weight in unitarised merging due to random
      // choice of reclustered/non-reclustered treatment and additional
      // sign for subtractive samples.
      double weight = pythia.info.weightValueByIndex();
      // Do not print zero-weight events.
      if (weight == 0.) continue;

      // Do not print broken/empty events.
      if (pythia.event.size() < 3) continue;

      // Work with unweighted events.
      double norm = xs/double(nEvent);
      // Work with weighted (LHA strategy=-4) events.
      if (abs(pythia.info.lhaStrategy()) == 4)
        norm = 1./double(nEvent);
      if (useLHA) norm /= MB2PB;

      // Accumulate cross section, including norm.
      pythia.info.weightContainerPtr->accumulateXsec(norm);

#if defined(HEPMC2) || defined(HEPMC3)
      // Optionally write HepMC events.
      if (doHepMC) {
        // Copy the weight names to HepMC.
        toHepMCPtr->setWeightNames(pythia.info.weightNameVector());
        // Fill HepMC event.
        toHepMCPtr->writeNextEvent(pythia);
      }
#endif

#ifdef RIVET
      // Optionally pass event to Rivet.
      if (doRivet) rivet();
#endif

    }

    // Break out of loop over iMerge if aborting.
    if (doAbort) break;

    // Print cross section and errors.
    pythia.stat();

    // Get cross section statistics for sample.
    double sigmaSample = pythia.info.weightContainerPtr->getSampleXsec()[0];
    double errorSample = pythia.info.weightContainerPtr->getSampleXsecErr()[0];

    cout << endl << " Cross section of sample " << iMerge << ": "
         << scientific << setprecision(8)
         << sigmaSample << " +- " << errorSample  << endl << endl;
  }

  // For Powheg matching, print veto information.
  if (doPowhegMatching) {
    cout << " PowhegHooks: ISR vetos: " << nVetoISR << endl
         << " PowhegHooks: FSR vetos: " << nVetoFSR << endl << endl;
  }

  // Get cross section statistics for total run.
  double sigmaTotal = pythia.info.weightContainerPtr->getTotalXsec()[0];
  double errorTotal = pythia.info.weightContainerPtr->getSampleXsecErr()[0];
  if (doAbort)
    cout << " Run was not completed owing to too many aborted events" << endl;
  else
    cout << " Inclusive cross section:   " << scientific << setprecision(8)
         << sigmaTotal << "  +-  " << errorTotal << " mb" << endl << endl;

  // Optionally finalise Rivet analysis.
#ifdef RIVET
  if (doRivet) rivet.done();
#endif

  // Optionally delete HepMC converter pointer.
#if defined(HEPMC2) || defined(HEPMC3)
  if (doHepMC) delete toHepMCPtr;
#endif

  // Done.
  return 0;

}
