// VinciaWeights.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the VinciaWeights class.

#include "Pythia8/VinciaWeights.h"

namespace Pythia8 {

//==========================================================================

// The VinciaWeights class.

//--------------------------------------------------------------------------

// Constants.

const double VinciaWeights::TINYANT       = 1.0e-10;
const double VinciaWeights::PACCEPTVARMAX = 0.99;
const double VinciaWeights::MINVARWEIGHT  = 0.01;

//--------------------------------------------------------------------------

// Initilize pointers.

bool VinciaWeights::initPtr(Info* infoPtrIn,
  VinciaCommon* vinComPtrIn) {
  infoPtr  = infoPtrIn;
  settingsPtr = infoPtr->settingsPtr;
  vinComPtr   = vinComPtrIn;
  isInitPtr   = true;
  return true;
}

//--------------------------------------------------------------------------

// Initialize.

void VinciaWeights::init() {

  // Check initPtr.
  if (!isInitPtr) {
    printOut("VinciaWeights::init","Error! pointers not initialized");
    return;
  }
  verbose          = settingsPtr->mode("Vincia:verbose");
  // TODO: uncertainty bands are disabled in this version of VINCIA.
  uncertaintyBands = false;
  varLabels.resize(0);
  nWeightsSav      = (uncertaintyBands ? 1+varLabels.size() : 1);

  // List of all keywords.
  allKeywords.clear(); allKeywords.resize(0);
  string ffKeys[5] = { ":", ":qqemit:", ":qgemit:", ":ggemit:", ":gxsplit:" };
  string ifKeys[8] = { ":", ":qqemit:", ":gqemit:", ":ggemit:", ":ggemit:",
                       ":qxsplit:", ":gxconv:", ":xgsplit:" };
  string iiKeys[6] = { ":", ":qqemit:", ":gqemit:", ":ggemit:", ":qxsplit:",
                       ":gxconv:" };
  for (int i = 0; i < 5; i++) allKeywords.push_back("ff"+ffKeys[i]+"murfac");
  for (int i = 0; i < 5; i++) allKeywords.push_back("ff"+ffKeys[i]+"cns");
  for (int i = 0; i < 8; i++) allKeywords.push_back("if"+ifKeys[i]+"murfac");
  for (int i = 0; i < 8; i++) allKeywords.push_back("if"+ifKeys[i]+"cns");
  for (int i = 0; i < 6; i++) allKeywords.push_back("ii"+iiKeys[i]+"murfac");
  for (int i = 0; i < 6; i++) allKeywords.push_back("ii"+iiKeys[i]+"cns");

  // Convert iAntPhys to keyword.
  iAntToKeyFSR.clear();
  iAntToKeyFSR[iQQemitFF]  = "qqemit";
  iAntToKeyFSR[iQGemitFF]  = "qgemit";
  iAntToKeyFSR[iGQemitFF]  = "qgemit";
  iAntToKeyFSR[iGGemitFF]  = "ggemit";
  iAntToKeyFSR[iGXsplitFF] = "gxsplit";
  iAntToKeyISR.clear();
  iAntToKeyISR[iQQemitII]  = "qqemit";
  iAntToKeyISR[iGQemitII]  = "gqemit";
  iAntToKeyISR[iGGemitII]  = "ggemit";
  iAntToKeyISR[iQXsplitII] = "qxsplit";
  iAntToKeyISR[iGXconvII]  = "gxconv";
  iAntToKeyISR[iQQemitIF]  = "qqemit";
  iAntToKeyISR[iQGemitIF]  = "qgemit";
  iAntToKeyISR[iGQemitIF]  = "gqemit";
  iAntToKeyISR[iGGemitIF]  = "ggemit";
  iAntToKeyISR[iQXsplitIF] = "qxsplit";
  iAntToKeyISR[iGXconvIF]  = "gxconv";
  iAntToKeyISR[iXGsplitIF] = "xgsplit";

  // Clean up the names of requested variations.
  for (int i = 0; i < (int)varLabels.size(); i++) {
    int strBegin = varLabels[i].find_first_not_of(" ");
    if ((i == 0) && (varLabels[i].find("{") != string::npos)) strBegin += 1;
    int strEnd   = varLabels[i].find_last_not_of(" ");
    if ((i == (int)varLabels.size()-1) && (varLabels[i].find("}") !=
        string::npos)) strEnd -= 1;
    int strRange = strEnd - strBegin + 1;
    varLabels[i] = toLower(varLabels[i].substr(strBegin, strRange));
  }

  // Parse names of requested variations and check for keywords.
  varKeys.clear(); varKeys.resize(varLabels.size());
  varVals.clear(); varVals.resize(varLabels.size());
  for (int i = 0; i < (int)varLabels.size(); i++) {
    varKeys[i].clear(); varKeys[i].resize(0);
    varVals[i].clear(); varVals[i].resize(0);
    string var      = varLabels[i];
    int    iEndName = var.find(" ", 0);
    string varName  = var.substr(0, iEndName);
    string left     = var.substr(iEndName, var.length());
    left     = left.substr(left.find_first_not_of(" "), left.length());
    varLabels[i] = varName;
    while (left.length() > 1) {
      int    iEnd = left.find("=", 0);
      int    iAlt = left.find(" ", 0);
      if ( (iAlt < iEnd) && (iAlt > 0) ) iEnd = iAlt;
      string key  = left.substr(0, iEnd);
      if (1+key.find_last_not_of(" ") < key.length())
        key = key.substr(0, 1+key.find_last_not_of(" "));
      varKeys[i].push_back(key);
      string val  = left.substr(iEnd+1, left.length());
      val  = val.substr(val.find_first_not_of(" "), val.length());
      val  = val.substr(0, val.find_first_of(" "));
      varVals[i].push_back(atof(val.c_str()));
      left = left.substr(iEnd+1, left.length());
      if (left.find_first_of(" ") >= left.length()) break;
      left = left.substr(left.find_first_of(" "), left.length());
      if (left.find_first_not_of(" ") >= left.length()) break;
      left = left.substr(left.find_first_not_of(" "), left.length());
    }
  }

  if (uncertaintyBands && (verbose >= 4)) {
    printOut("VinciaWeights", "List of variations, keywords and values:");
    for (int i = 0; i < (int)varLabels.size(); i++) {
      cout << "  " << varLabels[i] << " : ";
      for (int j=0; j<(int)varVals[i].size(); j++) {
        cout << " " << varKeys[i][j] << " -> " << varVals[i][j];
        if (j < (int)varVals[i].size() - 1) cout << ",";
      }
      cout << endl;
    }
  }

  // Safety check for non-sensible input.
  if (uncertaintyBands)
    for (int i = 0; i < (int)varKeys.size(); i++)
      for (int j = 0; j < (int)varKeys[i].size(); j++) {
        // Check input keywords against known ones.
        bool foundValidKey = false;
        for (int k = 0; k < (int)allKeywords.size(); k++)
          if (allKeywords[k] == varKeys[i][j]) {
            foundValidKey = true;
            break;
          }
        if (!foundValidKey)
          printOut("VinciaWeights", "Error! Unknown key " +
            varKeys[i][j] + " found, please check!");
      }

  nReportWeight     = 5;
  nReportedWeight   = 0;

  // Resize weight arrays to the relevant number of elements.
  weightsSav.resize(nWeightsSav);
  weightsOld.resize(nWeightsSav);
  weightsMax.resize(nWeightsSav);
  weightsMin.resize(nWeightsSav);
  weightSum.resize(nWeightsSav);
  weightSum2.resize(nWeightsSav);
  contribSum.resize(nWeightsSav);
  contribSum2.resize(nWeightsSav);

  // Initialize counters and flags.
  nTotWeights               = 0;
  nNonunityWeight           = 0;
  nNegativeWeight           = 0;
  nNonunityInitialWeight    = 0;
  nNegativeInitialWeight    = 0;
  nNonunityWeightNow        = 0;
  nNegativeWeightNow        = 0;
  nNonunityInitialWeightNow = 0;
  nNegativeInitialWeightNow = 0;
  didReweight               = false;
  firstCall                 = false;

  // Initialize weights to 0 (or 1 respectively).
  for (int iWeight=0; iWeight<nWeightsSav; iWeight++) {
    weightsSav[iWeight]  = 1.0;
    weightsOld[iWeight]  = 0.0;
    weightsMax[iWeight]  = 0.0;
    weightsMin[iWeight]  = 1.0e10;
    weightSum[iWeight]   = 0.0;
    weightSum2[iWeight]  = 0.0;
    contribSum[iWeight]  = 0.0;
    contribSum2[iWeight] = 0.0;
  }

}

//--------------------------------------------------------------------------

// Reset the weights, to be called at the beginning of each event.

void VinciaWeights::resetWeights(int nAccepted) {

  for (int iWeight = 0; iWeight < nWeightsSav; iWeight++) {
    weightsSav[iWeight]  = 1.0;
    weightsOld[iWeight]  = 0.0;
  }

  // Next doWeighting call will be the first one.
  firstCall = true;

  // No reweighting so far.
  didReweight = false;

  // Check if Pythia/Vincia vetoed event, restore counters.
  if (nAccepted < nTotWeights) {
    nTotWeights--;
    nNonunityWeight        -= nNonunityWeightNow;
    nNegativeWeight        -= nNegativeWeightNow;
    nNonunityInitialWeight -= nNonunityInitialWeightNow;
    nNegativeInitialWeight -= nNegativeInitialWeightNow;

    for (int iWeight = 0; iWeight < nWeightsSav; iWeight++) {
      weightSum[iWeight]  -= contribSum[iWeight];
      weightSum2[iWeight] -= contribSum2[iWeight];
    }
  }
  nNonunityWeightNow        = 0;
  nNegativeWeightNow        = 0;
  nNonunityInitialWeightNow = 0;
  nNegativeInitialWeightNow = 0;

  for (int iWeight = 0; iWeight<nWeightsSav; iWeight++) {
    contribSum[iWeight]  = 0.0;
    contribSum2[iWeight] = 0.0;
  }

}

//--------------------------------------------------------------------------

// Scaling functions.

void VinciaWeights::scaleWeightAll(double scaleFacIn) {
  for (int iWeight=0; iWeight<nWeightsSav; iWeight++)
    weightsSav[iWeight] *= scaleFacIn;
}

void VinciaWeights::scaleWeight(double scaleFacIn, int iWeightIn) {
  if (existsWeight(iWeightIn)) weightsSav[iWeightIn] *= scaleFacIn;
}

void VinciaWeights::scaleWeightVar(vector<double> pAccept, bool accept,
  bool isHard) {
  if (!uncertaintyBands) return;
  if (nWeightsSav <= 1) return;
  // Variations only pertain to hard process and resonance decays.
  if (!isHard) return;
  if (accept) scaleWeightVarAccept(pAccept);
  else scaleWeightVarReject(pAccept);
}

void VinciaWeights::scaleWeightVarAccept(vector<double> pAccept) {
  for (int iWeight = 1; iWeight<nWeightsSav; iWeight++) {
    double pAcceptVar = pAccept[iWeight];
    if (pAcceptVar > PACCEPTVARMAX) pAcceptVar = PACCEPTVARMAX;
    scaleWeight( pAcceptVar/pAccept[0], iWeight );
  }
}

void VinciaWeights::scaleWeightVarReject(vector<double> pAccept) {
  for (int iWeight = 1; iWeight<nWeightsSav; iWeight++) {
    double pAcceptVar = pAccept[iWeight];
    if (pAcceptVar > PACCEPTVARMAX) pAcceptVar = PACCEPTVARMAX;
    double reWeight = (1.0-pAcceptVar)/(1.0-pAccept[0]);
    if (reWeight < MINVARWEIGHT) reWeight = MINVARWEIGHT;
    scaleWeight( reWeight, iWeight );
  }
}

void VinciaWeights::scaleWeightEnhanceAccept(double enhanceFac) {
  if (enhanceFac == 1.0) return;
  else scaleWeightAll(1./enhanceFac);
}

void VinciaWeights::scaleWeightEnhanceReject(double pAcceptUnenhanced,
  double enhanceFac) {
  if (enhanceFac == 1.0) return;
  if (enhanceFac > 1.0) {
    double rRej =
      (1. - pAcceptUnenhanced/enhanceFac)/(1. - pAcceptUnenhanced);
    scaleWeightAll(rRej);
  } else {
    double rRej =
      (1. - pAcceptUnenhanced)/(1. - enhanceFac*pAcceptUnenhanced);
    scaleWeightAll(rRej);
  }
}

//--------------------------------------------------------------------------

// Helper function for keyword evaluation

int VinciaWeights::doVarNow(string keyIn, int iAntPhys, string type) {

  // Check variation for all branching types.
  string asKey = ":murfac", nsKey = ":cns";
  if (type + asKey == keyIn) return 1;
  if (type + nsKey == keyIn) return 2;
  // Check variation for specific branching type.
  map<int,string> keyCvt = (type == "ff" ? iAntToKeyFSR : iAntToKeyISR);
  if (type + ":" + keyCvt[iAntPhys] + asKey == keyIn) return 1;
  if (type + ":" + keyCvt[iAntPhys] + nsKey == keyIn) return 2;
  return -1;

}

//--------------------------------------------------------------------------

// Main function to perform the weighting.

void VinciaWeights::doWeighting() {

  // Check for non-unity and negative initial weights.
  double inW = 1.0;
  if (firstCall) {
    inW = infoPtr->weight();
    bool report = false;
    bool isLast = false;
    if ( (abs(inW-1.0) > TINY) || (inW < 0.0) ) {
      nReportedWeight++;
      if (nReportedWeight <= nReportWeight) {
        report = true;
        if (nReportedWeight == nReportWeight) isLast = true;
      }
      string printMessage = "Nonunity initial weight occurred, w = ";
      // Updated counters.
      if (inW < 0.0) {
        printMessage = "Negative initial weight occurred, w = ";
        nNegativeInitialWeight++;
        nNegativeInitialWeightNow++;
      } else {
        nNonunityInitialWeight++;
        nNonunityInitialWeightNow++;
      }
      // Always report in debugging mode.
      if ((report && verbose >= 3) || verbose >= 4)
        printOut("VinciaWeights", printMessage + num2str(inW) + (isLast ?
            ": further output suppressed" : ""));
    }
  }

  // Add weight correction to cumulative sum (!= w0 for systems in series).
  // Take into account Pythia's weight (only != 1 for first time in event).
  double w0       = inW*weightsSav[0];
  double wMC0     = weightsSav[0];
  double wOld     = weightsOld[0];
  weightSum[0]   += (w0 - wOld);
  weightSum2[0]  += (pow2(w0) - pow2(wOld));
  // Contribution of current event
  contribSum[0]  += (w0 - wOld);
  contribSum2[0] += (pow2(w0) - pow2(wOld));

  if (firstCall) {
    firstCall = false;
    // Add one to number of events.
    nTotWeights++;
  }

  // Check for non-unity and negative weights due to MC violations.
  bool foundWeightToReport = false;
  if (abs(wMC0-1.0) > TINY) {
    // Only check if wOld = 0 or if weight changed.
    if ( wOld < TINY || (abs(wOld) > TINY && (abs(wOld-w0) > TINY)) ) {
      didReweight         = true;
      foundWeightToReport = true;
      // Update counters.
      nNonunityWeight++;
      nNonunityWeightNow++;
    }
  }
  if (wMC0 < 0.0 && wOld >= 0.0) {
    didReweight         = true;
    foundWeightToReport = true;
    // Update counters.
    nNegativeWeight++;
    nNegativeWeightNow++;
  }
  bool report = false;
  bool isLast = false;
  if ( foundWeightToReport && (nReportedWeight <= nReportWeight) ) {
    report = true;
    if (nReportedWeight == nReportWeight) isLast = true;
  }
  // Always report in debugging mode.
  if ((report && verbose >= 3) || verbose >= 4) {
    string printMessage = (wMC0 > 0.0 ? "Nonunity MC reweight occurred, w = "
      : "Negative MC reweight occurred, w = ");
    printOut("VinciaWeights", printMessage + num2str(w0) + (isLast ?
        ": further output suppressed" : ""));
  }

  // Keep track of minimum and maximum weights.
  weightsMax[0] = max(weightsMax[0], w0);
  weightsMin[0] = min(weightsMin[0], w0);

  // Update the old weight.
  weightsOld[0] = w0;

  // Save uncertainty band weights.
  if (nWeightsSav > 1) {
    for (int iWeight = 1; iWeight < nWeightsSav; iWeight++) {
      double wVar          = inW*weightsSav[iWeight];
      double wVarOld       = weightsOld[iWeight];
      // Add weight correction to cumulative sum.
      weightSum[iWeight]   += (wVar - wVarOld);
      weightSum2[iWeight]  += (pow2(wVar) - pow2(wVarOld));
      // Contribution of current event.
      contribSum[iWeight]  += (wVar - wVarOld);
      contribSum2[iWeight] += (pow2(wVar) - pow2(wVarOld));
      // Keep track of minimum and maximum weights.
      weightsMax[iWeight]  = max(weightsMax[iWeight], wVar);
      weightsMin[iWeight]  = min(weightsMin[iWeight], wVar);
      // Update old weight.
      weightsOld[iWeight]  = wVar;
    }
  }

}

//==========================================================================


} // end namespace Pythia8
