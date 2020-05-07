// Vincia.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains global header information for Vincia.

#ifndef Pythia8_Vincia_H
#define Pythia8_Vincia_H

// Maths headers.
#include <limits>
#include <cmath>

// Include Pythia 8 headers.
#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/ShowerModel.h"

// Include Vincia headers.
#include "Pythia8/VinciaAntennaFunctions.h"
#include "Pythia8/VinciaCommon.h"
#include "Pythia8/VinciaFSR.h"
#include "Pythia8/VinciaISR.h"
#include "Pythia8/VinciaQED.h"
#include "Pythia8/VinciaMG5MEs.h"

// Define namespace inside which Vincia lives.
namespace Pythia8 {

//==========================================================================

// The Vincia class. Top-level handler class for the Vincia antenna
// shower model.

class Vincia : public ShowerModel {

public:

  // Constructor
  Vincia() = default;

  // Empty virtual destructor
  virtual ~Vincia() = default;

  // Initialize.
  bool init(MergingPtr mrgPtrIn, MergingHooksPtr mrgHooksPtrIn,
            PartonVertexPtr partonVertexPtrIn,
            WeightContainer* weightContainerPtrIn) override;

  // Function called from Pythia after the beam particles have been set up,
  // so that showers may be initialized after the beams are initialized.
  // Currently only dummy dunction.
  bool initAfterBeams() override { return true; }

  // Methods to get
  TimeShowerPtr  getTimeShower() const override { return timesPtr; }
  TimeShowerPtr  getTimeDecShower() const override { return timesDecPtr; }
  SpaceShowerPtr getSpaceShower() const override { return spacePtr; }
  MergingHooksPtr getMergingHooks() const override { return mergingHooksPtr; }
  MergingPtr getMerging() const override { return mergingPtr; }

  // Automatically set verbose level in all members.
  void setVerbose(int verboseIn);

  // Utilities for printing info and internal histograms.
  void printInfo() {
    timesPtr->printInfo(true);
    spacePtr->printInfo(true);
  }
  void printHistos() {
    timesPtr->printHistos();
  }
  void writeHistos(string fileName = "vincia", string lastName = "dat") {
    timesPtr->writeHistos(fileName, lastName);
  }
  const Hist& getDiagnosticHistogram(string name) {
    return timesPtr->getDiagnosticHistogram(name);
  }

  // Public Vincia objects.
  VinciaCommon          vinCom;
  Resolution            resolution;
  QEDShower             qedShower;
  Colour                colour;
  ResScaleHook          resScaleHook;
  VinciaWeights         vinWeights;
  MECs                  mecs;

  // Auxiliary objects.
  VinciaMG5MEs          mg5mes;
  Rambo                 rambo;

  // Vectors of antenna functions.
  DGLAP         dglap;
  AntennaSetFSR antennaSetFSR;
  AntennaSetISR antennaSetISR;

  // Pointers to Pythia classes.
  SusyLesHouches* slhaPtr;

 protected:

  // Method to initialise Vincia tune settings
  bool initTune(int iTune);

  // Members for the FSR and ISR showers.
  shared_ptr<VinciaFSR> timesPtr;
  shared_ptr<VinciaFSR> timesDecPtr;
  shared_ptr<VinciaISR> spacePtr;

 private:

  // Verbosity level.
  int verbose;

};

//==========================================================================

} // end Pythia8 namespace

#endif // end Pythia8_Vincia_H
