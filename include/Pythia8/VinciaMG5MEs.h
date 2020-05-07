// VinciaMG5MEs.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the functionality to interface C++ matrix elements
// generated with the PY8Kernels/MG5MES plugin to MadGraph 5

#ifndef Pythia8_VinciaMG5MEs_H
#define Pythia8_VinciaMG5MEs_H

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

// Include MG5 PY8MEs plugin headers.
#ifdef MG5MES
#include "PY8ME.h"
#include "PY8MEs.h"
#endif

// Include Vincia.
#include "Pythia8/VinciaCommon.h"

namespace Pythia8 {

#ifdef MG5MES
using namespace PY8MEs_namespace;
#endif

//==========================================================================

// Interface to MadGraph 5 matrix elements for matrix element corrections.

class VinciaMG5MEs {

public:

  // Constructor.
  VinciaMG5MEs() {isInitPtr = false; isInit = false; mg5libPtr = nullptr;}

  // Destructor.
  virtual ~VinciaMG5MEs() {if (mg5libPtr != nullptr) delete mg5libPtr;}

  // Set pointers to required PYTHIA 8 objects.
  void initPtr(Info* infoPtrIn, SusyLesHouches* slhaPtrIn,
    VinciaCommon* vinComPtrIn);

  // Initialise the MG5 model, parameters, and couplings.
  bool init();

  // Get pointer to matrix element, e.g. to check if process exists in
  // library. (Enabled when using the --with-mg5mes configure option.)
#ifdef MG5MES
  PY8ME* getProcess(vector<int> idIn, vector<int> idOut, set<int> sChan);
#endif

  // Get the matrix element squared for a particle state. The first
  // nIn particles must be the incoming ones. If MG5MES plugin not
  // linked, pure dummy implementation instead.
  double ME2(vector<Particle> state, int nIn);

  // Use ME2 to set helicities for a state. Takes a reference as
  // input and operates on it.
  bool selectHelicities(vector<Particle>& state, int nIn);

  // Set colour depth.
  void setColourDepth(int colourDepthIn) {colourDepth = colourDepthIn;}
  // Get colour depth.
  int getColourDepth() {return colourDepth;}

  // Convert a process label to a string, e.g. for printing to stdout.
  string makeLabel(vector<int>& id, int nIn, bool convertToNames=false);

  // Set verbosity level.
  void setVerbose(int verboseIn) {verbose = verboseIn;}

private:

  // Is initialized.
  bool isInitPtr, isInit;

  // Saved list of helicity components for last ME evaluated.
  map< vector<int> , double > me2hel;

  // Pointers to VINCIA and Pythia 8 objects.
  Info*           infoPtr;
  CoupSM*         coupSMPtr;
  ParticleData*   particleDataPtr;
  Rndm*           rndmPtr;
  Settings*       settingsPtr;
  VinciaCommon*   vinComPtr;
  SusyLesHouches* slhaPtr;

  // MG5 matrix-element library.
#ifdef MG5MES
  PY8MEs*         mg5libPtr;
#else
  int*            mg5libPtr;
#endif

  // Colour mode (0: leading colour, 1: Vincia colour).
  int colourDepth;

  // Verbosity level.
  int verbose;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_VinciaMG5MEs_H
