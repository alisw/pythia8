// Dire.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Stefan Prestel, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the top level Dire class.

#ifndef Pythia8_Dire_H
#define Pythia8_Dire_H

#define DIRE_VERSION "2.002"

// DIRE includes.
#include "Pythia8/DireSplittingLibrary.h"
#include "Pythia8/DireMerging.h"
#include "Pythia8/DireTimes.h"
#include "Pythia8/DireSpace.h"
#include "Pythia8/DireWeightContainer.h"
#include "Pythia8/DireHooks.h"

// Pythia includes.
#include "Pythia8/Info.h"
#include "Pythia8/Settings.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/Basics.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/ShowerModel.h"

#include <iostream>
#include <sstream>

namespace Pythia8 {

//==========================================================================

class Dire : public ShowerModel {

  public:

  Dire() : weightsPtr(nullptr), timesPtr(nullptr), timesDecPtr(nullptr),
    spacePtr(nullptr), splittings(nullptr), hooksPtr(nullptr),
    mergingPtr(nullptr), hardProcessPtr(nullptr), mergingHooksPtr(nullptr),
    hasOwnWeights(false), hasOwnTimes(false), hasOwnTimesDec(false),
    hasOwnSpace(false), hasOwnSplittings(false), hasOwnHooks(false),
    hasUserHooks(false), hasOwnHardProcess(false),
    hasOwnMergingHooks(false), initNewSettings(false), isInit(false),
    isInitShower(false), printBannerSave(true) { createPointers(); }

  Dire( MergingHooksPtr mergingHooksPtrIn, PartonVertexPtr partonVertexPtrIn)
    :  pythiaMergingHooksPtr(mergingHooksPtrIn),
    partonVertexPtr(partonVertexPtrIn), weightsPtr(nullptr),
    timesPtr(nullptr), timesDecPtr(nullptr), spacePtr(nullptr),
    splittings(nullptr), hooksPtr(nullptr),
    mergingPtr(nullptr), hardProcessPtr(nullptr),
    hasOwnWeights(false), hasOwnTimes(false), hasOwnTimesDec(false),
    hasOwnSpace(false), hasOwnSplittings(false), hasOwnHooks(false),
    hasUserHooks(false), hasOwnHardProcess(false),
    hasOwnMergingHooks(false), initNewSettings(false), isInit(false),
    isInitShower(false), printBannerSave(true) { createPointers(); }

 ~Dire() {
    if (hasOwnWeights)      delete weightsPtr;
    if (hasOwnSplittings)   delete splittings;
    if (hasOwnHardProcess)  delete hardProcessPtr;
  }

  // Flexible-use call at the beginning of each event in pythia.next().
  // Currently not used, but should be used for clearing some internal
  // bookkeeping that is otherewise reset in shower prepare functions.
  void beginEvent() {
    return;
  }

  // Flexible-use call at the end of each event in pythia.next().
  // Currently only to accumulate shower weights.
  void endEvent(PhysicsBase::Status status) {
    // No finalize in case of failure.
    if (status == INCOMPLETE) return;
    // Cross section conversion factor to undo annoying conversions done in
    // info.weight() calls.
    double CONVERTMB2PB = 1.;
    if (abs(infoPtr->lhaStrategy()) == 4) CONVERTMB2PB = 1e9;
    // Update the event weight by the Dire shower weight when relevant.
    // Retrieve the shower weight.
    weightsPtr->calcWeight(0.);
    weightsPtr->reset();
    double pswt = weightsPtr->getShowerWeight();
    // Multiply the shower weight to the event weight.
    double wt = infoPtr->weight()/CONVERTMB2PB;
    infoPtr->updateWeight(wt * pswt);
  }

  void createPointers();

  // Initialization function called before beams are set up.
  // Currently only to register objects as PhysicsBase (=initialize ptrs).
  bool init(MergingPtr, MergingHooksPtr, PartonVertexPtr, WeightContainer*) {
    subObjects.clear();
    if (mergingHooksPtr) {
      registerSubObject(*mergingHooksPtr);
    }
    if (mergingPtr) {
      registerSubObject(*mergingPtr);
    }
    if (timesPtr)    registerSubObject(*timesPtr);
    if (timesDecPtr) registerSubObject(*timesDecPtr);
    if (spacePtr)    registerSubObject(*spacePtr);
    return true;
  }

  // Initialization function called after beams are set up, used as main
  // initialization.
  bool initAfterBeams();

  void initTune();
  void initShowersAndWeights();
  void setup(BeamParticle* beamA, BeamParticle* beamB);
  void printBanner();

  TimeShowerPtr  getTimeShower() const   { return timesPtr; }
  TimeShowerPtr  getTimeDecShower() const { return timesDecPtr; }
  SpaceShowerPtr getSpaceShower() const   { return spacePtr; }
  MergingHooksPtr getMergingHooks() const { return mergingHooksPtr; }
  MergingPtr getMerging() const           { return mergingPtr; }

  MergingHooksPtr pythiaMergingHooksPtr;
  PartonVertexPtr partonVertexPtr;

  DireWeightContainer* weightsPtr;
  shared_ptr<DireTimes> timesPtr;
  shared_ptr<DireTimes> timesDecPtr;
  shared_ptr<DireSpace> spacePtr;
  DireSplittingLibrary* splittings;
  DireHooks* hooksPtr;

  DireInfo direInfo;

  // Pointer to Dire merging objects.
  shared_ptr<DireMerging>      mergingPtr;
  DireHardProcess*            hardProcessPtr;
  shared_ptr<DireMergingHooks> mergingHooksPtr;

  bool hasOwnWeights, hasOwnTimes, hasOwnTimesDec, hasOwnSpace,
       hasOwnSplittings, hasOwnHooks, hasUserHooks,
       hasOwnHardProcess, hasOwnMergingHooks;
  bool initNewSettings, isInit, isInitShower, printBannerSave;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_Dire_H
