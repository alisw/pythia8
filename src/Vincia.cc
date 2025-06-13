// Vincia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the headers) for the Vincia class.

#include "Pythia8/Vincia.h"
#include "Pythia8/Merging.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/Plugins.h"

namespace Pythia8 {

using namespace VinciaConstants;

//==========================================================================

// Vincia parton shower class.

//--------------------------------------------------------------------------

// Initialize.

bool Vincia::init(MergingPtr mrgPtrIn, MergingHooksPtr mrgHooksPtrIn,
                  PartonVertexPtr partonVertexPtrIn,
                  WeightContainer* weightContainerPtrIn) {

  // Verbosity output. Allow local debug level, else just use Print:verbosity.
  verbose = flag("Vincia:debug") ? VinciaConstants::DEBUG :
    mode("Print:verbosity");
  settingsPtr->addMode("Vincia:verbose",verbose,true,true,0,
    VinciaConstants::DEBUG);
  if (verbose >= VinciaConstants::DEBUG) {
    settingsPtr->mode("Print:verbosity",VinciaConstants::DEBUG);
    printOut(__METHOD_NAME__, "begin", DASHLEN);
  }

  // Create diagnostics pointer.
  diagnosticsPtr = make_shared<VinciaDiagnostics>();
  diagnosticsPtr->initInfoPtr(*infoPtr);
  if (verbose >= Logger::REPORT) diagnosticsPtr->start(__METHOD_NAME__);

  // Clear Vincia's register of PhysicsBase objects.
  subObjects.clear();

  bool vinciaOn      = (settingsPtr->mode("PartonShowers:model") == 2);
  bool doCutMerging  = flag("Merging:doCutBasedMerging");
  bool doKTMerging   = flag("Merging:doKTMerging");
  bool doMGMerging   = flag("Merging:doMGMerging");
  bool doDynMerging  = flag("Merging:doDynamicMerging");
  doMerging          = flag("Merging:doMerging");
  if ((doCutMerging || doKTMerging || doMGMerging || doDynMerging)
    && !doMerging) {
    doMerging = true;
    settingsPtr->flag("Merging:doMerging",true);
  }
  doMerging          = ( doMerging && vinciaOn );

  // Setup Vincia's merging if requested.
  if (doMerging) {
    // Ensure consistency in settings with merging.
    if (mode("Vincia:ewMode") > 2) {
      loggerPtr->WARNING_MSG(
        "reverting to default QED mode;"
        " EW shower not yet supported by merging");
      // Use readString so change is reapplied after Vincia tune setting.
      int ewModeDef = settingsPtr->modeDefault("Vincia:ewMode");
      settingsPtr->readString("Vincia:ewMode = "+to_string(ewModeDef));
    }
    if (flag("Vincia:interleaveResDec")) {
      loggerPtr->WARNING_MSG(
        "switching off interleaved resonance decays;"
        " not yet supported by merging.");
      // Must switch both Vincia and TimeShower flags off, since PartonLevel
      // uses the TimeShower one.
      settingsPtr->readString("Vincia:interleaveResDec = off");
      settingsPtr->readString("TimeShower:interleaveResDec = off");
    }

    // Set and register merging pointers
    mergingHooksPtr = make_shared<VinciaMergingHooks>();
    registerSubObject(*mergingHooksPtr);
    mrgHooksPtrIn = mergingHooksPtr;
    mergingPtr = make_shared<VinciaMerging>();
    registerSubObject(*mergingPtr);
    mrgPtrIn = mergingPtr;

    // Initialise Vincia's mergingHookPtr.
    mergingHooksPtr->init();

    if (!mergingHooksPtr->initSuccess()) {
      loggerPtr->ERROR_MSG("initialisation of MergingHooks failed");
      return false;
    }

    // Create Vincia's own userhook.
    shared_ptr<MergeResScaleHook> mergeResHookPtr =
      make_shared<MergeResScaleHook>(mergingHooksPtr);

    // Update userHooksPtr.
    if ( !userHooksPtr )
      userHooksPtr = mergeResHookPtr;
    else {
      shared_ptr<UserHooksVector> uhv =
        dynamic_pointer_cast<UserHooksVector>(userHooksPtr);
      if ( !uhv ) {
        uhv = make_shared<UserHooksVector>();
        uhv->hooks.push_back(userHooksPtr);
        userHooksPtr = uhv;
      }
      uhv->hooks.push_back(mergeResHookPtr);
    }

    // Update infoPtr's pointer to userhooks.
    infoPtr->userHooksPtr = userHooksPtr;
  }

  // Set weightContainerPtr and tell weightContainer where to find our weights.
  weightContainerPtr = weightContainerPtrIn;
  if (vinciaOn) weightContainerPtr->weightsShowerPtr = &vinWeights;

  // Create EW/QED Shower module(s).
  int ewMode = settingsPtr->mode("Vincia:EWmode");
  // Create the QED and EW shower pointers.
  ewShowerPtr      = std::make_shared<VinciaEW>();
  qedShowerHardPtr = std::make_shared<VinciaQED>();
  qedShowerSoftPtr = std::make_shared<VinciaQED>();

  if (vinciaOn && ewMode >= 3 && settingsPtr->flag("Vincia:EWOverlapVeto")) {
    // Initialize the overlap veto
    shared_ptr<VinciaEWVetoHook> EWvetoPtr = make_shared<VinciaEWVetoHook>();
    registerSubObject(*EWvetoPtr);
    EWvetoPtr->init(dynamic_pointer_cast<VinciaEW>(ewShowerPtr));

    // Update userHooksPtr.
    if ( !userHooksPtr ) {
      userHooksPtr = EWvetoPtr;
    }
    else {
      shared_ptr<UserHooksVector> uhv =
        dynamic_pointer_cast<UserHooksVector>(userHooksPtr);
      if ( !uhv ) {
        uhv = make_shared<UserHooksVector>();
        uhv->hooks.push_back(userHooksPtr);
        userHooksPtr = uhv;
      }
      uhv->hooks.push_back(EWvetoPtr);
    }

    // Update infoPtr's pointer to userhooks.
    infoPtr->userHooksPtr = userHooksPtr;
  }

  // Create and register VinciaFSR and VinciaISR instances.
  timesPtr = make_shared<VinciaFSR>() ;
  registerSubObject(*timesPtr);
  spacePtr = make_shared<VinciaISR>() ;
  registerSubObject(*spacePtr);
  timesDecPtr = timesPtr;

  // Set pointers in showers.
  timesPtr->initPtrs( mergingHooksPtr, partonVertexPtrIn,
    weightContainerPtr);
  spacePtr->initPtrs( mergingHooksPtr, partonVertexPtrIn,
    weightContainerPtr);

  // Pass verbose settings to members
  setVerbose(verbose);
  if (verbose >= Logger::REPORT) printOut(__METHOD_NAME__,
    "setting Vincia pointers...");

  // Init FSR shower pointers and default settings, beyond those set
  // by the non-virtual TimeShower::initPtr().
  timesPtr->initVinciaPtrs(&colour, spacePtr, &mecs,
    &resolution, &vinCom, &vinWeights);
  timesPtr->setDiagnosticsPtr(diagnosticsPtr);

  // Init ISR shower pointers and default settings, beyond those set
  // by the non-virtual SpaceShower::initPtr().
  spacePtr->initVinciaPtrs(&colour, timesPtr, &mecs,
    &resolution, &vinCom, &vinWeights);
  spacePtr->setDiagnosticsPtr(diagnosticsPtr);

  // FSR and ISR antenna sets.
  antennaSetFSR.initPtr(infoPtr, &dglap);
  antennaSetISR.initPtr(infoPtr, &dglap);

  // Hand antenna set pointers to shower and matching objects.
  timesPtr->initAntPtr(&antennaSetFSR);
  spacePtr->initAntPtr(&antennaSetISR);
  mecs.initAntPtr(&antennaSetFSR, &antennaSetISR);

  // Set SLHA pointer
  slhaPtr = coupSUSYPtr->slhaPtr;
  if (slhaPtr == nullptr)
    printOut(__METHOD_NAME__, "Warning: SLHA pointer is null pointer.");

  // Load the matrix element correction plugin.
  string melib = settingsPtr->word("Vincia:MEplugin");
  if (melib.size() > 0)
    mg5mes = make_plugin<ExternalMEs>(
      "libpythia8mg5" + melib + ".so", "ExternalMEsMadgraph",
      nullptr, settingsPtr, loggerPtr);
  if (mg5mes == nullptr) mg5mes = make_shared<ExternalMEs>();

  // Pass pointers on to objects that require them.
  rambo.initPtr(rndmPtr);
  vinCom.initPtr(infoPtr);
  resolution.initPtr(settingsPtr, infoPtr, &vinCom);
  if (mg5mes != nullptr) mg5mes->initPtrs(infoPtr);
  mecs.initPtr(infoPtr, mg5mes, &vinCom, &resolution);
  colour.initPtr(infoPtr);
  vinWeights.initPtr(infoPtr, &vinCom);

  // Initialize pointers in EW shower modules.
  // Set EW/QED Shower module in timesPtr and spacePtr.
  // QED shower for hard interaction + resonance decays.
  qedShowerHardPtr->initPtr(infoPtr, &vinCom);
  timesPtr->setQEDShowerHardPtr(qedShowerHardPtr);
  spacePtr->setQEDShowerHardPtr(qedShowerHardPtr);

  // QED shower for MPI and hadronisation.
  qedShowerSoftPtr->initPtr(infoPtr, &vinCom);
  timesPtr->setQEDShowerSoftPtr(qedShowerSoftPtr);
  spacePtr->setQEDShowerSoftPtr(qedShowerSoftPtr);

  // Electroweak shower.
  ewShowerPtr->initPtr(infoPtr, &vinCom);
  // Save some information on resonances locally,
  // and modify particleDataPtr if doing resonance decays.
  if (ewMode >= 3) ewShowerPtr->load();
  timesPtr->setEWShowerPtr(ewShowerPtr);
  spacePtr->setEWShowerPtr(ewShowerPtr);

  // If Vincia is on, allow to override some Pythia settings by
  // Vincia-specific ones.
  if (vinciaOn) {
    // PartonLevel only checks TimeShower:interleaveResDec, so set that to
    // agree with the corresponding Vincia flag.
    bool interleaveResDec = settingsPtr->flag("Vincia:interleaveResDec");
    settingsPtr->flag("TimeShower:interleaveResDec",interleaveResDec);
  }

  // Initialise Vincia auxiliary classes (showers initialised by Pythia).
  vinCom.init();
  resolution.init();
  colour.init();
  vinWeights.init( doMerging );

  // MECs depend on Pythia/SLHA Couplings.
  mecs.init();
  if (!mecs.isInitialised()) {
    loggerPtr->ERROR_MSG("failed to initialise MECs");
    return false;
  }

  // Print VINCIA header and list of parameters
  if (verbose >= Logger::NORMAL && vinciaOn) timesPtr->header();

  // Diagnostics
  if (verbose >= Logger::REPORT) diagnosticsPtr->stop(__METHOD_NAME__);

  // Verbose output.
  if (verbose >= VinciaConstants::DEBUG)
    printOut(__METHOD_NAME__, "end", DASHLEN);
  return true;

}

//--------------------------------------------------------------------------

// Automatically set verbose level in all members.

void Vincia::setVerbose(int verboseIn) {

  verbose = verboseIn;
  vinCom.setVerbose(verboseIn);
  resolution.setVerbose(verboseIn);
  timesPtr->setVerbose(verboseIn);
  spacePtr->setVerbose(verboseIn);
  colour.setVerbose(verboseIn);
  mecs.setVerbose(verboseIn);
  if (doMerging) {
    mergingHooksPtr->setVerbose(verboseIn);
    mergingPtr->setVerbose(verboseIn);
  }
  if (ewShowerPtr != nullptr) ewShowerPtr->setVerbose(verboseIn);
  if (qedShowerHardPtr != nullptr) qedShowerHardPtr->setVerbose(verboseIn);
  if (qedShowerSoftPtr != nullptr) qedShowerSoftPtr->setVerbose(verboseIn);

  // If Vincia:Debug is on, also set Logger verbosity to REPORT.
  if (verbose >= VinciaConstants::DEBUG)
    loggerPtr->setVerbosity(Logger::REPORT);

}

//==========================================================================

} // end namespace Pythia8
