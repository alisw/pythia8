// Weights.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains classes that keep track of event weights.

#ifndef Pythia8_Weights_H
#define Pythia8_Weights_H

#include "Pythia8/Basics.h"
#include "Pythia8/LHEF3.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/SharedPointers.h"

namespace Pythia8 {

//==========================================================================

// This is a base class to store weight information in a way that allows
// unified access in the structure that contains all event generation weights
// information (WeightContainer below). The main purpuse of this class is to
// supply convenience features to make defining new categories of weights easy.
// All weights should inherit from this base class. The specialized classes
// may then contain specialized functions, or only act as a glorified
// book-keeping struncture.

class WeightsBase {

  public:

  // Reset all internal values;
  virtual void clear() { return; };

  // Store the current event information.
  virtual void init(vector<double> /*weightValues*/,
    vector<string> /*weightNames*/) { return; }

  // Function to return processed weights to weight container, e.g. if
  // weights should be combined before proceeding.
  virtual void collectWeightValues(vector<double>& outputWeights,
    double norm = 1.);

  // Similar function to return processed weight names.
  virtual void collectWeightNames(vector<string>& outputNames);

  // Get the stored information.
  // Direcly use storage members here in the base class,
  // and access those through non-virtual getters.
  // Note: NOT opting for a map data structure, since information will anyway
  // have to be serialized in output.
  vector<double> weightValues;
  vector<string> weightNames;
  string getWeightsName(int iPos)  { return weightNames[iPos];   }
  double getWeightsValue(int iPos) { return weightValues[iPos];  }
  int getWeightsSize()             { return weightValues.size(); }

  // Function to create a new, synchronized, pair of weight name and value.
  void bookWeight(string name, double defaultValue = 1.)  {
    weightNames.push_back(name);
    weightValues.push_back(defaultValue);
  }

  // Functions to set values of weights.
  void setValueByIndex(int iPos, double val) {
    weightValues[iPos] = val;
  }
  void setValueByName(string name, double val) {
    // Use existing functions: Find index of name, then set by index.
    int iPos = findIndexOfName(name);
    setValueByIndex(iPos, val);
  }

  // Function to find the index of a given entry in the weightNames vector,
  // e.g., to be able to synchronize with the weightValues vector.
  int findIndexOfName(string name) {
    vector<string>::iterator it
      = find(weightNames.begin(), weightNames.end(), name);
    return distance(weightNames.begin(), it);
  }

};

//==========================================================================

// This is a short example class to collect information on parton shower
// weights into a container class that can be part of Weight, which
// in turn is part of InfoHub.

class WeightsShower : public WeightsBase {

  public:

  // Reset all internal values;
  void clear();

  // Store the current event information.
  void init(vector<double> weights, vector<string> names);

  //// Function to group weights and return processed weights to container.
  //void collectWeightValues(vector<double>& outputWeights);
  //void collectWeightNames(vector<string>& outputNames);

  // Functions to set values of weights.
  void reweightValueByIndex(int iPos, double val);
  void reweightValueByName(string name, double val);

  void replaceWhitespace( vector<string>& namesIn);

};

//==========================================================================

// This class collects information on weights generated in the heavy ion
// framework, now limited to a single weight (per event) from
// impact parameter sampling.

class WeightsHeavyIon : public WeightsBase {

  public:

  // Reset all internal values;
  void clear();

  // Store the current event information.
  void init(vector<double> weights, vector<string> names);

  // Functions to set values of weights.
  void reweightValueByIndex(int iPos, double val);
  void reweightValueByName(string name, double val);

};

//==========================================================================

// This is a short example class to collect information on Les Houches Event
// weights into a container class that can be part of Weight, which
// in turn is part of InfoHub.

class WeightsLHEF : public WeightsBase {

  public:

  // Reset all internal values;
  void clear();

  // Store the current event information.
  void init(vector<double> weights_detailed_vecIn,
    vector<string> weights_detailed_name_vecIn);

  // Function to return processed weights to weight container, e.g. if
  // weights should be combined before proceeding.
  void collectWeightValues(vector<double>& outputWeights,
     double norm = 1.);
  void collectWeightNames(vector<string>& outputNames);

  // Convert weight names in MadGraph5 convention to the convention outlined
  // in https://arxiv.org/pdf/1405.1067.pdf, page  162ff.
  vector<string> weightnames_lhef2hepmc(
    vector<string> weights_detailed_name_vecIn);

};

//==========================================================================

// This is a container class to collect all event generation weight
// information into a wrapper which is in turn is part of InfoHub. In this
// way, we could avoid cluttering InfoHub.

class WeightContainer {

  public:

  // Default constructor only ensures that members are initialized with
  // sensible default values.
  WeightContainer() : weightNominal(1.0) {}

  // The nominal Pythia weight
  double weightNominal;
  void setWeightNominal( double weightNow );

  // First example of a weight subcategory.
  WeightsLHEF          weightsLHEF;

  // Other possible sub-categories:
  WeightsShower        weightsPS;

  WeightsHeavyIon     weightsHI;

  // Other possible sub-categories:
  //WeightsMerging       weightInfoMerging;
  //WeightsHadronization weightInfoHadronization;

  // Functions to retrieve information stored in the subcategory members.
  int numberOfWeights();
  double weightValueByIndex(int key=0);
  string weightNameByIndex(int key=0);

  // Function to return the vector of weight values, combining all weights from
  // all subcategories.
  // Currently, only the nominal weight and LHEF weights are
  // considered. Internal Pythia weights should also be included eventually.
  vector<double> weightValueVector();

  // Function to return the vector of weight names, combining all names from
  // all subcategories, cf. weightValueVector function.
  vector<string> weightNameVector();

  // Reset all members to default stage.
  void clear();

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Weights_H
