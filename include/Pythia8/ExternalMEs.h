// ExternalMEs.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Peter Skands, Stefan Prestel, Philip Ilten, Torbjorn
// Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the functionality to interface external matrix
// elements for matrix element corrections to parton showers.

#ifndef Pythia8_ExternalMEs_H
#define Pythia8_ExternalMEs_H

// Include Pythia headers.
#include "Pythia8/Basics.h"
#include "Pythia8/PythiaComplex.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyLesHouches.h"

namespace Pythia8 {

//==========================================================================

// Base class for external matrix-element interfaces.

class ExternalMEs {

public:

  // Destructor.
  ExternalMEs() = default;
  virtual ~ExternalMEs() = default;

  // Initialisers for pointers.
  virtual void initPtrs(Info* infoPtrIn);

  // Initialisers.
  virtual bool init() {return false;}
  virtual bool initVincia(Info* /*infoPtrIn*/) {return false;}
  virtual bool initDire(Info* /*infoPtrIn*/, string /*card*/) {return false;}

  // Methods to check availability of matrix elements for list of in/out
  // ID codes, an event (ignoring any event entries before iBeg), or a vector
  // of particles.
  virtual bool isAvailable(vector<int> /*idIn*/, vector<int> /*idOut*/) {
    return false;}
  virtual bool isAvailable(const Event& /*event*/) {return false;}
  virtual bool isAvailable(const Event& /*event*/, const int /*iBeg*/) {
    return false;}
  virtual bool isAvailable(const vector<Particle>& /*state*/) {return false;}

  // Calculate the matrix element squared for a particle state or event
  // (ignoring any event entries before iBeg).
  virtual double calcME2(const vector<Particle>& /*state*/) {return 0;}
  virtual double calcME2(const Event& /*event*/, const int /*iBeg*/) {
    return 0;}

  // Setters.
  virtual void setColourMode(int colModeIn) {
    colMode = colModeIn;}
  virtual void setHelicityMode(int helModeIn) {
    helMode = helModeIn;}
  virtual void setIncludeSymmetryFac(bool doInclIn) {
    inclSymFac = doInclIn;}
  virtual void setIncludeHelicityAvgFac(bool doInclIn) {
    inclHelAvgFac = doInclIn;}
  virtual void setIncludeColourAvgFac(bool doInclIn) {
    inclColAvgFac = doInclIn;}

  // Getters.
  virtual int  colourMode() {return colMode;}
  virtual int  helicityMode() {return helMode;}
  virtual bool includeSymmetryFac() {return inclSymFac;}
  virtual bool includeHelicityAvgFac() {return inclHelAvgFac;}
  virtual bool includeColourAvgFac() {return inclColAvgFac;}
  virtual map<vector<int>, double> getHelicityAmplitudes() {return me2hel;}

protected:

  // Fill a vector of IDs, from an event, starting from entry i = iBeg.
  void fillIds(const Event& event, vector<int>& in, vector<int>& out,
    int iBeg = 3) const;
  // Fill a vector of momenta, from an event, starting from entry i = iBeg.
  void fillMoms(const Event& event, vector<Vec4>& p, int iBeg = 3) const;
  // Fill a vector of colors, from an event, starting from entry i = iBeg.
  void fillCols(const Event& event, vector<int>& colors, int iBeg = 3) const;
  // Return the momenta, from an event, starting from entry i = iBeg.
  vector<vector<double> > fillMoms(const Event& event, int iBeg = 3) const;

  // Colour mode (0: strict LC, 1: LC, 2: LC sum, 3: FC).
  int colMode{1};

  // Saved list of helicity components for last ME evaluated.
  map<vector<int>, double> me2hel;

  // Helicity mode (0: explicit helicity sum, 1: implicit helicity sum).
  int helMode{1};

  // Symmetry and averaging factors.
  bool inclSymFac{false}, inclHelAvgFac{false}, inclColAvgFac{false};

  // Pointers to VINCIA and Pythia 8 objects.
  Info*           infoPtr{};
  Logger*         loggerPtr{};
  CoupSM*         coupSMPtr{};
  ParticleData*   particleDataPtr{};
  Rndm*           rndmPtr{};
  Settings*       settingsPtr{};
  SusyLesHouches* slhaPtr{};

  // Is initialized.
  bool isInitPtr{false}, isInit{false};

};

//==========================================================================

// A helicity sampler using external matrix elements.

class HelicitySampler {

 public:
  // Constructor, destructor, and assignment.
  HelicitySampler() {isInitPtr = false;}
  ~HelicitySampler() = default;

  // Initialise pointers to required Pythia objects.
  void initPtrs(ExternalMEsPtr mePluginPtrIn, Rndm* rndmPtrIn) {
    mePluginPtr = mePluginPtrIn;
    rndmPtr     = rndmPtrIn;
    isInitPtr   = true;}

  // Set helicities for a particle state.
  bool selectHelicities(vector<Particle>& state, bool force);

 private:

  // Pointers to Pythia objects.
  ExternalMEsPtr mePluginPtr;
  Rndm* rndmPtr;

  // Flag whether we have all pointers.
  bool isInitPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_ExternalMEs_H
