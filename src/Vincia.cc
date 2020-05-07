// VinciaCommon.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the headers) for the Vincia class.

#include "Pythia8/Vincia.h"
#include "Pythia8/Merging.h"
#include "Pythia8/MergingHooks.h"

namespace Pythia8 {

//==========================================================================

// Vincia parton shower class.

//--------------------------------------------------------------------------

// Initialize.

bool Vincia::init(MergingPtr mrgPtrIn, MergingHooksPtr mrgHooksPtrIn,
                  PartonVertexPtr partonVertexPtrIn,
                  WeightContainer* weightContainerPtrIn) {

  // Clear Vincia's register of PhysicsBase objects
  subObjects.clear();

  // Set and register merging pointers
  mergingPtr = mrgPtrIn;
  if ( mergingPtr ) registerSubObject(*mergingPtr);
  mergingHooksPtr = mrgHooksPtrIn;
  if ( mergingHooksPtr ) registerSubObject(*mergingHooksPtr);

  // Create and register VinciaFSR and VinciaISR instances
  timesPtr = make_shared<VinciaFSR>() ;
  registerSubObject(*timesPtr);
  spacePtr = make_shared<VinciaISR>() ;
  registerSubObject(*spacePtr);
  timesDecPtr = timesPtr;

  // Set pointers in showers.
  timesPtr->initPtrs( mergingHooksPtr, partonVertexPtrIn,
    weightContainerPtrIn);
  spacePtr->initPtrs( mergingHooksPtr, partonVertexPtrIn,
    weightContainerPtrIn);

  // Verbosity level.
  setVerbose(settingsPtr->mode("Vincia:verbose"));
  if (verbose >= quiteloud) printOut(__METHOD_NAME__, "setting pointers...");

  // Init FSR shower pointers and default settings, beyond those set
  // by the non-virtual TimeShower::initPtr().
  timesPtr->initVinciaPtrs(&colour,spacePtr,&qedShower,&mecs,
    &resolution, &vinCom,&vinWeights);

  // Init ISR shower pointers and default settings, beyond those set
  // by the non-virtual SpaceShower::initPtr().
  spacePtr->initVinciaPtrs(&colour,timesPtr,&qedShower,&mecs,
    &resolution, &vinCom,&vinWeights);

  // FSR and ISR antenna sets.
  antennaSetFSR.initPtr(infoPtr, &dglap);
  antennaSetISR.initPtr(infoPtr, &dglap);

  // QED Shower module.
  qedShower.initPtr(infoPtr, &vinCom);

  // Hand antenna set pointers to shower and matching objects.
  timesPtr->initAntPtr(&antennaSetFSR);
  spacePtr->initAntPtr(&antennaSetISR);
  mecs.initAntPtr(&antennaSetFSR, &antennaSetISR);

  // Set SLHA pointer
  slhaPtr = coupSUSYPtr->slhaPtr;

  // Pass pointers on to objects that require them.
  resolution.initPtr(settingsPtr);
  rambo.initPtr(rndmPtr);
  vinCom.initPtr(infoPtr);
  mg5mes.initPtr(infoPtr, slhaPtr, &vinCom);
  mecs.initPtr(infoPtr, &mg5mes, &vinCom);
  colour.initPtr(infoPtr);
  vinWeights.initPtr(infoPtr, &vinCom);

  // Now set tune parameters
  bool vinciaOn = settingsPtr->mode("PartonShowers:model") == 2;
  int baseTune = settingsPtr->mode("Vincia:Tune");
  if (vinciaOn && baseTune >= 0) {
    // Store user-specified settings before overwriting with tune parameters
    vector<string> userSettings = settingsPtr->getReadHistory();
    if (initTune(baseTune)) {
      // Reapply user settings
      for (int i=0; i<(int)userSettings.size(); ++i) {
        string lineNow      = userSettings[i];
        string lineNowLower = toLower(lineNow);
        if (lineNowLower.find("tune:ee") == string::npos &&
          lineNowLower.find("tune:pp") == string::npos)
          settingsPtr->readString(lineNow);
      }
    }
  }

  // Initialise Vincia auxiliary classes (showers initialised by Pythia)
  vinCom.init();
  resolution.init();
  colour.init();
  vinWeights.init();

  // MECs depend on Pythia/SLHA Couplings
  mecs.init();

  // Print VINCIA header and list of parameters
  if (verbose >= 1 && vinciaOn) timesPtr->header();

  // Verbose output
  if(verbose >= veryloud) printOut(__METHOD_NAME__, "end --------------");
  return true;

}

//--------------------------------------------------------------------------

// Vincia tune settings.

bool Vincia::initTune(int iTune) {

  // iTune = 0 : default Vincia tune from Pythia 8.302
  if (iTune == 0) {
    // Z fractions in string breaks
    settingsPtr->parm("StringZ:aLund            ", 0.55 );
    settingsPtr->parm("StringZ:bLund            ", 0.78 );
    settingsPtr->parm("StringZ:aExtraDiquark    ", 0.90 );
    // Z fractions for heavy quarks
    settingsPtr->parm("StringZ:rFactC           ", 1.15 );
    settingsPtr->parm("StringZ:rFactB           ", 0.85 );
    // pT in string breaks
    settingsPtr->parm("StringPT:sigma",            0.305);
    settingsPtr->parm("StringPT:enhancedFraction", 0.01);
    settingsPtr->parm("StringPT:enhancedWidth",    2.0);
    // String breakup flavour parameters
    settingsPtr->parm("StringFlav:probStoUD     ", 0.205);
    settingsPtr->parm("StringFlav:mesonUDvector ", 0.42 );
    settingsPtr->parm("StringFlav:mesonSvector  ", 0.53 );
    settingsPtr->parm("StringFlav:mesonCvector  ", 1.3  );
    settingsPtr->parm("StringFlav:mesonBvector  ", 2.2  );
    settingsPtr->parm("StringFlav:probQQtoQ     ", 0.077);
    settingsPtr->parm("StringFlav:probSQtoQQ    ", 1.0  );
    settingsPtr->parm("StringFlav:probQQ1toQQ0  ", 0.025);
    settingsPtr->parm("StringFlav:etaSup        ", 0.5  );
    settingsPtr->parm("StringFlav:etaPrimeSup   ", 0.1  );
    settingsPtr->parm("StringFlav:decupletSup   ", 1.0  );
    settingsPtr->parm("StringFlav:popcornSpair  ", 0.75 );
    settingsPtr->parm("StringFlav:popcornSmeson ", 0.75 );
    // Primordial kT
    settingsPtr->parm("BeamRemnants:primordialKThard ", 0.4 );
    settingsPtr->parm("BeamRemnants:primordialKTsoft ", 0.25);
    // MB/UE tuning parameters (MPI)
    // Use a "low" alphaS and 2-loop running everywhere, also for MPI
    settingsPtr->parm("SigmaProcess:alphaSvalue ", 0.119);
    settingsPtr->mode("SigmaProcess:alphaSorder ", 2);
    settingsPtr->parm("MultiPartonInteractions:alphaSvalue", 0.119);
    settingsPtr->mode("MultiPartonInteractions:alphaSorder", 2);
    settingsPtr->parm("MultiPartonInteractions:pT0ref     ", 2.24);
    settingsPtr->parm("MultiPartonInteractions:expPow     ", 1.75);
    settingsPtr->parm("MultiPartonInteractions:ecmPow     ", 0.21);
    // Use PYTHIA 8's baseline CR model
    settingsPtr->flag("ColourReconnection:reconnect", true);
    settingsPtr->parm("ColourReconnection:range    ", 1.75);
    // Diffraction: switch off Pythia's perturbative MPI
    // (colours in diffractive systems not yet handled by Vincia)
    settingsPtr->parm("Diffraction:mMinPert", 1000000.0);
    return true;
  }
  // Unknown iTune.
  else return false;
}

//--------------------------------------------------------------------------

// Automatically set verbose level in all members.

void Vincia::setVerbose(int verboseIn) {

  verbose = verboseIn;
  vinCom.setVerbose(verboseIn);
  resolution.setVerbose(verboseIn);
  timesPtr->setVerbose(verboseIn);
  qedShower.setVerbose(verboseIn);
  spacePtr->setVerbose(verboseIn);
  colour.setVerbose(verboseIn);
  mg5mes.setVerbose(verboseIn);
  mecs.setVerbose(verboseIn);

}

//==========================================================================

} // end namespace Pythia8
