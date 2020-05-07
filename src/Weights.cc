// Weights.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Weight classes.

#include "Pythia8/Weights.h"
#include <limits>

namespace Pythia8 {

//==========================================================================

// WeightsBase class.

//--------------------------------------------------------------------------

// Function to return processed weights to weight container, e.g. if
// weights should be combined before proceeding.
void WeightsBase::collectWeightValues(vector<double>& outputWeights,
  double norm) {
  for (int iwt=0; iwt < getWeightsSize(); ++iwt) {
    double value = getWeightsValue(iwt)*norm;
    outputWeights.push_back(value);
  }
  return;
}

//--------------------------------------------------------------------------

// Similar function to return processed weight names.
void WeightsBase::collectWeightNames(vector<string>& outputNames) {
  for (int iwt=0; iwt < getWeightsSize(); ++iwt) {
    string name  = getWeightsName(iwt);
    outputNames.push_back(name);
  }
  return;
}

//==========================================================================

// WeightsShower class.

//--------------------------------------------------------------------------

// Reset all internal values;

void WeightsShower::clear() {
  for (size_t i=0; i < weightValues.size(); ++i) weightValues[i] = 1.;
}

//--------------------------------------------------------------------------

// Store the current event information.

void WeightsShower::init(vector<double> weights, vector<string> names) {

  // Remember the nominal weight, since this might be needed for splitting
  // enhancement hanlding.
  bookWeight("Baseline");

  replaceWhitespace(names);

  for (size_t i=0; i < weights.size(); ++i) bookWeight(names[i], weights[i]);

}

//--------------------------------------------------------------------------

// Replace whitespace with underscore in wieght names, so that the names
// transferred to HepMC do not contain whitespace.
void WeightsShower::replaceWhitespace( vector<string>& namesIn) {
  vector<string> ret;
  for (size_t i=0; i < namesIn.size(); ++i) {
    string name=namesIn[i];
    replace(name.begin(), name.end(), ' ', '_');
    ret.push_back(name);
    namesIn[i] = name;
  }
}

//--------------------------------------------------------------------------

// Functions to set values of weights.

void WeightsShower::reweightValueByIndex(int iPos, double val) {
  weightValues[iPos] *= val;
}

//--------------------------------------------------------------------------

void WeightsShower::reweightValueByName(string name, double val) {
  // Use existing functions: Find index of name, then set by index.
  int iPos = findIndexOfName(name);
  reweightValueByIndex(iPos, val);
}


//==========================================================================

// WeightsHeavyIon class.

//--------------------------------------------------------------------------

// Reset all internal values;

void WeightsHeavyIon::clear() {
  for (size_t i=0; i < weightValues.size(); ++i) weightValues[i] = 1.;
}

//--------------------------------------------------------------------------

// Store the current event information.

void WeightsHeavyIon::init(vector<double> weights, vector<string> names) {
  for (size_t i=0; i < weights.size(); ++i) bookWeight(names[i], weights[i]);
}

//--------------------------------------------------------------------------

// Functions to set values of weights.

void WeightsHeavyIon::reweightValueByIndex(int iPos, double val) {
  weightValues[iPos] *= val;
}

//--------------------------------------------------------------------------

void WeightsHeavyIon::reweightValueByName(string name, double val) {
  // Use existing functions: Find index of name, then set by index.
  int iPos = findIndexOfName(name);
  reweightValueByIndex(iPos, val);
}

//==========================================================================

// WeightsLHEF class.

// Reset all internal values;
void WeightsLHEF::clear() {
  weightValues.resize(0);
  weightNames.resize(0);
}

//--------------------------------------------------------------------------

// Store the current event information.
void WeightsLHEF::init(vector<double> weights_detailed_vecIn,
  vector<string> weights_detailed_name_vecIn){
  weightValues = weights_detailed_vecIn;
  weightNames  = weightnames_lhef2hepmc(weights_detailed_name_vecIn);
}

//--------------------------------------------------------------------------

// Function to return processed weights to weight container.
void WeightsLHEF::collectWeightValues(vector<double>& ret, double) {

  // Attach the LHEF weights, starting with well-defined MUF and MUR
  // variations, and then followed by any other LHEF weight.
  for (int iwt=0; iwt < getWeightsSize(); ++iwt) {
    double value = getWeightsValue(iwt);
    string name  = getWeightsName(iwt);
    if (name.find("MUR") == string::npos || name.find("MUF") == string::npos)
      continue;
    ret.push_back(value);
  }
  for (int iwt=0; iwt < getWeightsSize(); ++iwt) {
    double value = getWeightsValue(iwt);
    string name  = getWeightsName(iwt);
    if (name.find("MUR") != string::npos || name.find("MUF") != string::npos)
      continue;
    ret.push_back(value);
  }

  // Done.
  return;
}

//--------------------------------------------------------------------------

// Function to return processed weight names to weight container.
void WeightsLHEF::collectWeightNames(vector<string>& ret) {

  // Attach the LHEF weights, starting with well-defined MUF and MUR
  // variations, and then followed by any other LHEF weight.
  for (int iwt=0; iwt < getWeightsSize(); ++iwt) {
    string name  = getWeightsName(iwt);
    if (name.find("MUR") == string::npos || name.find("MUF") == string::npos)
      continue;
    ret.push_back(name);
  }
  for (int iwt=0; iwt < getWeightsSize(); ++iwt) {
    string name  = getWeightsName(iwt);
    if (name.find("MUR") != string::npos || name.find("MUF") != string::npos)
      continue;
    ret.push_back(name);
  }

  // Done.
  return;
}

//--------------------------------------------------------------------------

// Convert weight names in MadGraph5 convention to the convention outlined
// in https://arxiv.org/pdf/1405.1067.pdf, page  162ff.
vector<string> WeightsLHEF::weightnames_lhef2hepmc(
  vector<string> weights_detailed_name_vecIn) {
  vector<string> ret;
  for (size_t i=0; i < weights_detailed_name_vecIn.size(); ++i) {
    string name=weights_detailed_name_vecIn[i];
    if (name=="1001") name="MUR1.0_MUF1.0";
    if (name=="1002") name="MUR1.0_MUF2.0";
    if (name=="1003") name="MUR1.0_MUF0.5";
    if (name=="1004") name="MUR2.0_MUF1.0";
    if (name=="1005") name="MUR2.0_MUF2.0";
    if (name=="1006") name="MUR2.0_MUF0.5";
    if (name=="1007") name="MUR0.5_MUF1.0";
    if (name=="1008") name="MUR0.5_MUF2.0";
    if (name=="1009") name="MUR0.5_MUF0.5";
    ret.push_back(name);
  }
  return ret;
}

//==========================================================================

// The WeightContainer class.

//--------------------------------------------------------------------------

void WeightContainer::setWeightNominal(double weightNow) {
  weightNominal = weightNow;
}

//--------------------------------------------------------------------------

// Functions to retrieve the stored information.

int WeightContainer::numberOfWeights() {
  // Currently, only LHEF weights considered, plus one entry for the
  // nominal Pythia weight.
  return (1 + weightsLHEF.getWeightsSize()
            + weightsPS.getWeightsSize()
            + weightsHI.getWeightsSize());
}

double WeightContainer::weightValueByIndex(int key) {
  vector<double> values = weightValueVector();
  return values[key];
}

string WeightContainer::weightNameByIndex(int key) {
  vector<string> names = weightNameVector();
  return names[key];
}

//--------------------------------------------------------------------------

// Function to return the vector of weight values, combining all weights from
// all subcategories.
// Currently, only the nominal weight and LHEF weights are
// considered. Internal Pythia weights should also be included eventually.
vector<double> WeightContainer::weightValueVector() {
  vector<double> ret;

  // The very first entry in the vector should always be the nominal weight.
  ret.push_back(weightNominal);

  // Let all weights attach the weight values to the return vector.
  weightsLHEF.collectWeightValues(ret);
  weightsPS.collectWeightValues(ret,weightNominal);
  weightsHI.collectWeightValues(ret);

  // Note: Here, we could still manipulate the weight vector, e.g. to
  // combine weights from different sources.

  // Done
  return ret;

}

// Function to return the vector of weight names, combining all names from
// all subcategories, cf. weightValueVector function.
vector<string> WeightContainer::weightNameVector() {
  vector<string> ret;

   // The very first entry in the vector should always be the nominal weight.
  ret.push_back("Weight");

  // Let all weights attach the weight names to the return vector.
  weightsLHEF.collectWeightNames(ret);
  weightsPS.collectWeightNames(ret);
  weightsHI.collectWeightNames(ret);

  // Note: Here, we could still manipulate the weight name vector, e.g. to
  // combine weights from different sources.

  // Done
  return ret;

}

//--------------------------------------------------------------------------

// Reset all members to default status.

void WeightContainer::clear() {
  weightNominal = 1.;
  weightsLHEF.clear();
  weightsPS.clear();
}

//==========================================================================

} // end namespace Pythia8
