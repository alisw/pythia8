// LowEnergyProcess.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for low-energy hadronic collisions, as used for rescattering.

#ifndef Pythia8_LowEnergyProcess_H
#define Pythia8_LowEnergyProcess_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationSystems.h"
#include "Pythia8/Info.h"
#include "Pythia8/MiniStringFragmentation.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PhysicsBase.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// LowEnergyProcess class.
// Is used to describe the low-energy collision between two hadrons.

class LowEnergyProcess : public PhysicsBase {

public:

  // Constructor. Still to be expanded with further default values.
  LowEnergyProcess() : stringFragPtr(), ministringFragPtr() {}

  // Initialize the class.
  bool init(StringFragmentation* stringFragPtrIn,
    MiniStringFragmentation* ministringFragPtrIn);

  // Produce outgoing primary hadrons from collision of incoming pair.
  bool collide( int i1, int i2, int type, Event& event, Vec4 vtx = Vec4() );

  // Event record to handle hadronization.
  Event         leEvent;

private:

  // Constants: could only be changed in the code itself.
  static const int MAXLOOP;
  static const double MASSREDUCERATE, MDIFFMIN, ALPHAPRIME;

  // Parameters of the generation process.
  double probStoUD, fracEtass, fracEtaPss, xPowMes, xPowBar, xDiqEnhance,
         sigmaQ, mStringMin;

  // Properties of the current collision. 1 or 2 is two incoming hadrons.
  // "c" or "ac" is colour or anticolour component of hadron.
  bool   isBaryon1, isBaryon2;
  int    sizeOld, id1, id2, idc1, idac1, idc2, idac2;
  double m1, m2, eCM, sCM, z1, z2, mT1, mT2, mA, mB,
         mc1, mac1, px1, py1, pTs1, mTsc1, mTsac1, mTc1, mTac1,
         mc2, mac2, px2, py2, pTs2, mTsc2, mTsac2, mTc2, mTac2;

  // Pointer to the generator for normal string fragmentation.
  StringFragmentation* stringFragPtr;

  // Pointer to the generator for special low-mass string fragmentation.
  MiniStringFragmentation* ministringFragPtr;

  // Separate configuration for simple collisions.
  ColConfig     simpleColConfig;

  // Handle inelastic nondiffractive collision.
  bool nondiff();

  // Handle elastic and diffractive collisions.
  bool eldiff( int type);

  // Handle annihilation collisions.
  bool annihilation();

  // Simple version of hadronization for low-energy hadronic collisions.
  bool simpleHadronization(Event& event, bool isDiff = false);

  // Split up hadron A or B into a colour pair, with masses and pT values.
  bool splitA( double redMpT);
  bool splitB( double redMpT);

  // Split a hadron inte a colour and an anticolour part.
  pair< int, int> splitFlav( int id);

  // Choose relative momentum of colour and anticolour constituents in hadron.
  double splitZ( int iq1, int iq2, double mRat1, double mRat2);

  // Estimate lowest possible mass state for flavour combination.
  double mThreshold( int iq1, int iq2);

  // Pick slope b of exp(b * t) for elastic and diffractive events.
  double bSlope( int type);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LowEnergyProcess_H
