// FlavorVariations.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Stephen Mrenna, Christian Bierlich, Philip Ilten,
// Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef Pythia8_FlavorVariations_H
#define Pythia8_FlavorVariations_H

// Include Pythia headers.
#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// Calculate the weight for an event given a different set of flavor
// parameters used in the hadronization.

class FlavorVariations {

public:

  // Constructor, given an intialized Pythia object.
  FlavorVariations(Settings &settings) : FlavorVariations(
    settings.parm("StringFlav:ProbQQtoQ"),
    settings.parm("StringFlav:ProbStoUD"),
    settings.parm("StringFlav:ProbSQtoQQ"),
    settings.parm("StringFlav:ProbQQ1toQQ0")) {}

  // Constructor, given the default base parameters.
  FlavorVariations(double xi, double rho, double x, double y) :
    pythia("", false) {
    pythia.settings.flag("ProcessLevel:all", false);
    pythia.settings.flag("Print:quiet", true);
    pythia.settings.flag("VariationFrag:flav", true);
    pythia.settings.parm("StringFlav:ProbQQtoQ", xi);
    pythia.settings.parm("StringFlav:ProbStoUD", rho);
    pythia.settings.parm("StringFlav:ProbSQtoQQ", x);
    pythia.settings.parm("StringFlav:ProbQQ1toQQ0", y);
    pythia.settings.addMVec(key, vector<int>(14, 0), false, false, 0, 0);
    pythia.init();
  }

  // Read and write string break counts.
  vector<int> read(string breaks) {
    pythia.settings.readString(key + " = " + breaks);
    return pythia.settings.mvec(key);}
  string write(const vector<int>& breaks) {
    string out = "{";
    for (const int& val : breaks) out += toString(val) + ",";
    return out.substr(0, out.length() - 1) + "}";}

  // Calculate the derived parameters.
  vector<double> parms(double xi, double rho, double x, double y) {
    return pythia.info.weightContainerPtr
      ->weightsFragmentation.flavParms(xi, rho, x, y);}

  // Calculate the weight.
  double weight(const vector<double>& parms, const vector<int>& breaks) {
    return pythia.info.weightContainerPtr
      ->weightsFragmentation.flavWeight(parms, breaks);}

private:

  // Pythia object.
  Pythia pythia;

  // Key to serialize the string breaks.
  string key{"VariationFrag:breaks"};

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_FlavorVariations_H
