// VinciaWeights.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains header information for the VinciaWeights class.

#ifndef Vincia_VinciaWeights_H
#define Vincia_VinciaWeights_H

#include "Pythia8/Info.h"
#include "Pythia8/Settings.h"
#include "Pythia8/VinciaCommon.h"

namespace Pythia8 {

class VinciaWeights {

public:

  // Friends for internal private members.
  friend class VinciaFSR;
  friend class VinciaISR;

  // Constructor.
  VinciaWeights() {}

  // Destructor.
  ~VinciaWeights() {}

  // Initilize pointers.
  bool initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn);

  // Initilize.
  void init();

  // Check if particular weight exists.
  bool existsWeight(int iWeightIn = 0) {
    if (iWeightIn >= 0 && iWeightIn < nWeightsSav) return true;
    else return false;}

  // Get weight of current event weight. Normally 1 for set 0 (user
  // settings). Use iWeightIn > 0 to access uncertainty weights.
  double weight(int iWeightIn = 0) {
    if (existsWeight(iWeightIn)) return weightsSav[iWeightIn];
    return 0.0;}

  // Access the weight labels.
  string weightLabel(int iWeightIn = 0) {
    if (iWeightIn == 0) return "no variation";
    if (existsWeight(iWeightIn) && iWeightIn-1 < (int)varLabels.size())
      return varLabels[iWeightIn-1];
    return "";}

  // Functions to access by VINCIA plugin.
  int  nWeights()            {return nWeightsSav;}
  bool reweightingOccurred() {return didReweight;}

  // Reset the weights, to be called at the beginning of each event.
  void resetWeights(int nAccepted);

  // Scale all event weights.
  void scaleWeightAll(double scaleFacIn);

  // Scale a particular event weight.
  void scaleWeight(double scaleFacIn, int iWeightIn = 0);

  // Scale the uncertainty band weights.
  void scaleWeightVar(vector<double> pAccept, bool accept, bool isHard);

  // Scale the uncertainty band weights if branching is accepted.
  void scaleWeightVarAccept(vector<double> pAccept);

  // Scale the uncertainty band weights if branching is rejected.
  void scaleWeightVarReject(vector<double> pAccept);

  // Enhanced kernels: reweight if branching is accepted.
  void scaleWeightEnhanceAccept(double enhanceFac = 1.);

  // Enhanced kernels: reweight if branching is rejected.
  void scaleWeightEnhanceReject(double pAcceptUnenhanced,
    double enhanceFac = 1.);

  // Helper function for keyword evaluation.
  int doVarNow(string keyIn, int iAntPhys, string type) ;

  // Helper function for antenna function.
  double ant(double antIn, double cNSIn) {return (antIn+cNSIn);}

  // Main function to perform the weighting.
  void doWeighting();

private:

  // Verbosity.
  int verbose;

  // Pointers.
  Info*         infoPtr;
  Settings*     settingsPtr;
  VinciaCommon* vinComPtr;

  // Internal flag.
  bool isInitPtr;

  // Constants.
  static const double TINYANT, PACCEPTVARMAX, MINVARWEIGHT;

  // Parameters taken from settings.
  bool   uncertaintyBands;
  vector<string> varLabels;
  vector< vector<string> > varKeys;
  vector< vector<double> > varVals;

  // Helper parameters.
  vector<string> allKeywords;
  map<int,string> iAntToKeyFSR, iAntToKeyISR;
  double nWeightsSav;
  int    nReportWeight, nReportedWeight;
  bool   doAlphaSvar, doFiniteVar;

  // Weights (from input events, enhancement, MC violations, or variations).
  vector<double> weightsSav;
  vector<double> weightsOld;
  vector<double> weightsMax, weightsMin;
  // Weight sums.
  vector<double> weightSum, weightSum2;
  // Contribution of current event.
  vector<double> contribSum, contribSum2;

  // Counter for total number of weights = nr of events.
  int nTotWeights;
  // Counter for initial and main weight (MC violation only).
  int nNonunityWeight, nNegativeWeight,
    nNonunityInitialWeight, nNegativeInitialWeight;
  // Counter for current event.
  int nNonunityWeightNow, nNegativeWeightNow,
    nNonunityInitialWeightNow, nNegativeInitialWeightNow;

  // Flag to tell if reweighting due to MC violations occured.
  bool didReweight;

  // Flag if doWeighting called for the first time.
  bool firstCall;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_VinciaWeights_H
