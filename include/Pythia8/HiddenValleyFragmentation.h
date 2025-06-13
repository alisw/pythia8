// HiddenValleyFragmentation.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the classes for Hidden-Valley fragmentation.

#ifndef Pythia8_HiddenValleyFragmentation_H
#define Pythia8_HiddenValleyFragmentation_H

#include "Pythia8/FragmentationModel.h"
#include "Pythia8/MiniStringFragmentation.h"
#include "Pythia8/StringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The HVStringFlav class is used to select HV-quark and HV-hadron flavours.

class HVStringFlav : public StringFlav {

public:

  // Constructor.
  HVStringFlav() : separateFlav(), nFlav(), probFlav(), probDiquark(),
    probVector(), probKeepEta1(), sumProbFlav(), probKeepLast(),
    probVecEta1() {}

  // Destructor.
  ~HVStringFlav() {}

  // Initialize data members.
  void init() override;

  // Pick a new flavour (including diquarks) given an incoming one.
  FlavContainer pick(FlavContainer& flavOld, double, double, bool) override;

  // Combine two flavours (including diquarks) to produce a hadron.
  int combine(FlavContainer& flav1, FlavContainer& flav2) override;

  // Lightest flavour-neutral meson.
  int idLightestNeutralMeson() override { return 4900111; }

private:

  // Initialization data, to be read from Settings.
  bool   separateFlav;
  int    nFlav;
  vector<double> probFlav;
  double probDiquark, probVector, probKeepEta1, sumProbFlav, probKeepLast,
    probVecEta1;

};

//==========================================================================

// The HVStringPT class is used to select select HV transverse momenta.

class HVStringPT : public StringPT {

public:

  // Constructor.
  HVStringPT() : setabsigma(), rescalebsigma() {}

  // Destructor.
  ~HVStringPT() {}

  // Feed in extra parameters for HV handling, not part of base class.
  void preinit( int setabsigmaIn, double rescalebsigmaIn);

  // Initialize data members.
  void init() override;

private:

  // Initialization data fed in from HiddenValleyFragmentation::init().
  int    setabsigma;
  double rescalebsigma;

};

//==========================================================================

// The HVStringZ class is used to sample the HV fragmentation function f(z).

class HVStringZ : public StringZ {

public:

  // Constructor.
  HVStringZ() : setabsigma(), rescalebsigma(), mVecRatio(),
    rFactBowler() {}

  // Destructor.
  virtual ~HVStringZ() {}

  // Feed in extra parameters for HV handling, not part of base class.
  void preinit( int setabsigmaIn, double rescalebsigmaIn, double mVecRatioIn);

  // Initialize data members.
  bool init() override;

  // Fragmentation function: top-level to determine parameters.
  double zFrag( int idOld, int idNew = 0, double mT2 = 1.) override;

  // Parameters for stopping in the middle; for now hardcoded.
  virtual double stopMass()    override {return stopM;}
  virtual double stopNewFlav() override {return stopNF;}
  virtual double stopSmear()   override {return stopS;}

private:

  // Initialization data, from preinit, Settings and ParticleData.
  int    setabsigma;
  double rescalebsigma, mVecRatio;
  vector<double> rFactBowler;

};

//==========================================================================

// The HiddenValleyFragmentation class contains the routines
// to fragment a Hidden Valley partonic system.

class HiddenValleyFragmentation : public FragmentationModel {

public:

  // Constructor.
  HiddenValleyFragmentation() : doHVfrag(false), separateFlav(), nFlav(),
    hvOldSize(), hvNewSize(), idEnd1(), idEnd2(), mhvMeson(), mhvMin(),
    mHVvecMin(), mSys(), ihvParton() {}

  // Initialize and save pointers.
  bool init(StringFlav* flavSelPtrIn = nullptr, StringPT* pTSelPtrIn = nullptr,
    StringZ* zSelPtrIn = nullptr, FragModPtr fragModPtrIn = nullptr) override;

  // Fragment the event.
  bool fragment(int iSub, ColConfig& colConfig, Event& event,
    bool isDiff = false, bool systemRecoil = true) override;

protected:

  virtual void onInitInfoPtr() override {
    registerSubObject(hvStringFrag);
    registerSubObject(hvMinistringFrag);
    registerSubObject(hvFlavSel);
    registerSubObject(hvPTSel);
    registerSubObject(hvZSel);
  }

private:

  // Data mambers.
  bool          doHVfrag, separateFlav;
  int           nFlav, hvOldSize, hvNewSize, idEnd1, idEnd2;
  double        mhvMeson, mhvMin[9], mHVvecMin, mSys;
  vector<int>   ihvParton;

  // Configuration of colour-singlet systems.
  ColConfig     hvColConfig;

  // Temporary event record for the Hidden Valley system.
  Event         hvEvent;

  // The generator class for Hidden Valley string fragmentation.
  StringFragmentation hvStringFrag;

  // The generator class for special low-mass HV string fragmentation.
  MiniStringFragmentation hvMinistringFrag;

  // Pointers to classes for flavour, pT and z generation in HV sector.
  HVStringFlav hvFlavSel;
  HVStringPT   hvPTSel;
  HVStringZ    hvZSel;

  // Extract HV-particles from event to hvEvent.
  bool extractHVevent(Event& event);

  // Trace HV-colours of HV-partons.
  bool traceHVcols();

  // Collapse of low-mass system to one HV-meson.
  bool collapseToMeson();

  // Insert HV particles from hvEvent to event.
  bool insertHVevent(Event& event);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HiddenValleyFragmentation_H
