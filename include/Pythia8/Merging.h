// Merging.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Merging: Wrapper class to interface matrix element merging schemes with
//          Pythia

#ifndef Pythia8_Merging_H
#define Pythia8_Merging_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/History.h"
#include "Pythia8/Info.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonLevel.h"
#include "Pythia8/PhysicsBase.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"

namespace Pythia8 {

//==========================================================================

// Merging is a wrapper class for the interface of matrix element merging and
// Pythia8.

class Merging : public PhysicsBase {

public:

  // Constructor.
  Merging() : PhysicsBase(), lhaPtr(nullptr), trialPartonLevelPtr(),
    mergingHooksPtr(), tmsNowMin() {}

  // Destructor.
  virtual ~Merging(){}

  // Initialisation function for internal use inside Pythia source code
  void initPtrs( MergingHooksPtr mergingHooksPtrIn,
    PartonLevel* trialPartonLevelPtrIn) {
    trialPartonLevelPtr = trialPartonLevelPtrIn;
    mergingHooksPtr = mergingHooksPtrIn;
  }

  // Initialisation function for internal use inside Pythia source code
  virtual void init();

  // Function to print statistics.
  virtual void statistics();

  // Function to steer different merging prescriptions.
  virtual int mergeProcess( Event& process);

  // Runtime interface functions for communication with aMCatNLO
  // Function to retrieve shower scale information (to be used to set
  // scales in aMCatNLO-LHEF-production.
  virtual void getStoppingInfo(double scales[100][100],
    double masses[100][100]);

  // Function to retrieve if any of the shower scales would not have been
  // possible to produce by Pythia.
  virtual void getDeadzones(bool dzone[100][100]);

  // Function to generate Sudakov factors for MCatNLO-Delta.
  virtual double generateSingleSudakov (double pTbegAll, double pTendAll,
    double m2dip, int idA, int type, double s = -1., double x = -1.);

  LHEF3FromPythia8Ptr lhaPtr;
  void setLHAPtr(LHEF3FromPythia8Ptr lhaUpIn) { lhaPtr = lhaUpIn; }

protected:

  // Make Pythia class friend
  friend class Pythia;

  // Pointer to trial PartonLevel object
  PartonLevel* trialPartonLevelPtr;

  // Pointer to trial MergingHooks object
  MergingHooksPtr mergingHooksPtr;

  // Minimal value found for the merging scale in events.
  double tmsNowMin;
  static const double TMSMISMATCH;

  // Minimum allowed weight value to prevent division by zero.
  static const double MINWGT;

  // Function to perform CKKW-L merging on the event.
  int mergeProcessCKKWL( Event& process);

  // Function to perform UMEPS merging on the event.
  int mergeProcessUMEPS( Event& process);

  // Function to perform NL3 NLO merging on the event.
  int mergeProcessNL3( Event& process);

  // Function to perform UNLOPS merging on the event.
  int mergeProcessUNLOPS( Event& process);

  // Function to apply the merging scale cut on an input event.
  bool cutOnProcess( Event& process);

  // Clear all information stored in the runtime interface to aMCatNLO.
  void clearInfos() {
    stoppingScalesSave.clear();
    mDipSave.clear();
    radSave.clear();
    emtSave.clear();
    recSave.clear();
    isInDeadzone.clear();
  }

  // Store all information required for the runtime interface to aMCatNLO.
  int clusterAndStore(Event& process);

  // Helper function to be able to extract all shower scales by checking
  // all dipoles. Relevant only to runtime aMC@NLO interface.
  void getDipoles( int iRad, int colTag, int colSign,
    const Event& event, vector<pair<int,int> >& dipEnds);

  // Saved information about shower stopping scales, dipole masses,
  // dipole ends, and whether or not a clustering is in the shower
  // deadzone. Relevant only to runtime aMC@NLO interface.
  vector<double> stoppingScalesSave, mDipSave;
  vector<int> radSave, emtSave, recSave;
  vector<bool> isInDeadzone;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Merging_H
