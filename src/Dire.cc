// Dire.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Stefan Prestel, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the top level
// Dire class.

#include "Pythia8/Dire.h"

namespace Pythia8 {

//==========================================================================

// The Dire wrapper class.

//--------------------------------------------------------------------------

void Dire::createPointers() {

  // Construct showers.
  if (!weightsPtr)  {
    hasOwnWeights  = true;
    weightsPtr     = new DireWeightContainer(settingsPtr);
  }
  if (!timesPtr) {
    hasOwnTimes    = true;
    timesPtr       = make_shared<DireTimes>( mergingHooksPtr, partonVertexPtr);
  }
  if (!spacePtr) {
    hasOwnSpace    = true;
    spacePtr       = make_shared<DireSpace>( mergingHooksPtr, partonVertexPtr);
  }
  if (!timesDecPtr) {
    hasOwnTimesDec = true;
    timesDecPtr    = make_shared<DireTimes>( mergingHooksPtr, partonVertexPtr);
  }
  if (!mergingPtr) {
    mergingPtr    = make_shared<DireMerging>();
  }
  if (!hardProcessPtr) {
    hasOwnHardProcess = true;
    hardProcessPtr    = new DireHardProcess();
  }
  if (!mergingHooksPtr) {
    hasOwnMergingHooks = true;
    mergingHooksPtr    = make_shared<DireMergingHooks>();
  }

}

//--------------------------------------------------------------------------

void Dire::initTune() {

  initNewSettings = true;

  // Get tune id.
  int iTune = settingsPtr->mode("Dire:Tune");

  // Default tune.
  if (iTune == 1) settingsPtr->readString("include = tunes/Dire.cmnd");

  // For new U(1) splittings, teach Pythia new particles.
  if ( settingsPtr->flag("TimeShower:U1newShowerByL")
    || settingsPtr->flag("TimeShower:U1newShowerByQ")
    || settingsPtr->flag("SpaceShower:U1newShowerByL")
    || settingsPtr->flag("SpaceShower:U1newShowerByQ"))
    settingsPtr->readString("include = tunes/DireU1.cmnd");

}

//--------------------------------------------------------------------------

void Dire::initShowersAndWeights() {

  if (isInitShower) return;

  // Construct showers.
  if (!weightsPtr)  {
    hasOwnWeights  = true;
    weightsPtr     = new DireWeightContainer(settingsPtr);
  }
  if (!timesPtr) {
    hasOwnTimes    = true;
    timesPtr       = make_shared<DireTimes>( mergingHooksPtr, partonVertexPtr);
  }
  if (!spacePtr) {
    hasOwnSpace    = true;
    spacePtr       = make_shared<DireSpace>( mergingHooksPtr, partonVertexPtr);
  }
  if (!timesDecPtr) {
    hasOwnTimesDec = true;
    timesDecPtr    = make_shared<DireTimes>( mergingHooksPtr, partonVertexPtr);
  }
  if (!mergingPtr) {
    mergingPtr    = make_shared<DireMerging>();
  }
  if (!hardProcessPtr) {
    hasOwnHardProcess = true;
    hardProcessPtr    = new DireHardProcess();
  }
  if (!mergingHooksPtr) {
    hasOwnMergingHooks = true;
    mergingHooksPtr    = make_shared<DireMergingHooks>();
  }

  mergingHooksPtr->setHardProcessPtr(hardProcessPtr);
  mergingHooksPtr->init();

  timesPtr->setWeightContainerPtr(weightsPtr);
  spacePtr->setWeightContainerPtr(weightsPtr);
  timesDecPtr->setWeightContainerPtr(weightsPtr);

  isInitShower = true;

}

//--------------------------------------------------------------------------

void Dire::setup(BeamParticle* beamA, BeamParticle* beamB) {

  if (isInit) return;

  // Initialise library of splitting functions.
  if (!splittings) {
    hasOwnSplittings = true;
    splittings       = new DireSplittingLibrary();
  }

  // If Pythia has, for ominous reasons, not initialized the spacelike shower,
  // retry to initialize from timelike shower beams.
  if ( !spacePtr->isInit() && timesPtr->isInit()
    && beamA != 0 && beamB != 0)
    spacePtr->init( beamA, beamB );

  // Reinitialise showers to ensure that pointers are correctly set.
  timesPtr->reinitPtr(infoPtr, mergingHooksPtr, splittings, &direInfo);
  spacePtr->reinitPtr(infoPtr, mergingHooksPtr, splittings, &direInfo);
  timesDecPtr->reinitPtr(infoPtr, mergingHooksPtr, splittings, &direInfo);

  // Reset Pythia masses if necessary.
  if ( settingsPtr->flag("ShowerPDF:usePDFmasses")
    && ( beamA != nullptr || beamB != nullptr) ) {
    for (int i=1; i <= 5; ++i) {
      // Try to get masses from the hadron beams.
      double mPDF = (abs(beamA->id()) > 30)
                  ? beamA->mQuarkPDF(i)
                  : (abs(beamB->id()) > 30)
                    ? beamB->mQuarkPDF(i) : -1.0;
      // If there are no hadron beams, get the masses from either beam.
      if (beamA != nullptr && mPDF < 0.)
        mPDF = beamA->mQuarkPDF(i);
      if (beamB != nullptr && mPDF < 0.)
        mPDF = beamB->mQuarkPDF(i);
      if (mPDF > -1.) {
        stringstream resetMass;
        resetMass << i << ":m0 = " << mPDF;
        settingsPtr->readString(resetMass.str());
      }
    }
  }

  // Switch off all showering and MPI when estimating the cross section,
  if (hooksPtr)
    hooksPtr->initPtr( infoPtr, beamA, beamB);

  splittings->setKernelHooks(hooksPtr);

  // Initialise splitting function library here so that beam pointers
  // are already correctly initialised.
  splittings->init( infoPtr, beamA, beamB, &direInfo, hooksPtr);

  // Feed the splitting functions to the showers.
  splittings->setTimesPtr(timesPtr);
  splittings->setTimesDecPtr(timesDecPtr);
  splittings->setSpacePtr(spacePtr);

  // Initialize splittings in showers again (!), now that splittings are
  // properly set up.
  timesDecPtr->initSplits();
  timesPtr->initSplits();
  spacePtr->initSplits();

  weightsPtr->initPtrs(beamA, beamB, settingsPtr, infoPtr, &direInfo);
  timesDecPtr->initVariations();
  timesPtr->initVariations();
  spacePtr->initVariations();
  if (mergingPtr) mergingPtr->initPtrs( weightsPtr, timesPtr,
    spacePtr, &direInfo);


}

//--------------------------------------------------------------------------

//bool Dire::init(BeamParticle* beamA, BeamParticle* beamB) {
bool Dire::initAfterBeams() {

  if (isInit) return true;

  // Construct showers.
  initShowersAndWeights();

  // Initialize Dire tune settings.
  initTune();

  if ( settingsPtr->flag("Dire:doMerging")
    || settingsPtr->flag("Dire:doMECs")
    || settingsPtr->flag("Dire:doMEM")) {
    settingsPtr->flag("Merging:doMerging",true);
    settingsPtr->flag("Merging:useShowerPlugin",true);
  }

  if ( settingsPtr->flag("Dire:doMECs")
    || settingsPtr->flag("Dire:doMEM")) {
    settingsPtr->parm("Merging:TMS",0.0);
  }

  // No QED radiation by default until properly validated
  settingsPtr->flag("TimeShower:QEDshowerByQ",false);
  settingsPtr->flag("TimeShower:QEDshowerByL",false);
  settingsPtr->flag("SpaceShower:QEDshowerByQ",false);
  settingsPtr->flag("SpaceShower:QEDshowerByL",false);

  // Setup weight container (after user-defined enhance factors have been read)
  weightsPtr->initPtrs(beamAPtr, beamBPtr, settingsPtr, infoPtr, &direInfo);
  weightsPtr->setup();
  setup(beamAPtr, beamBPtr);
  isInit = true;

  printBannerSave = printBannerSave && !settingsPtr->flag("Print:quiet");
  if (printBannerSave) printBanner();
  printBannerSave = false;

  return isInit;

}

void Dire::printBanner() {

  cout << "\n"
       << " *---------------  Welcome to the DIRE parton shower "
       << "  -------------*\n"
       << " |                                                "
       << "                  |\n"
       << " | Please consider citing Eur.Phys.J. C75 (2015)"
       << " 9, 461             |\n"
       << " | if you use this program for scientific purposes."
       << "                 |\n";
  cout << " |                                                "
       << "                  |\n";
  cout << " *----------------------------------------"
       << "--------------------------*" << endl;

}

//==========================================================================

} // end namespace Pythia8
