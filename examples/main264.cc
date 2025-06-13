// main264.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>
//          Stephen Mrenna <mrenna@fnal.gov>
//          Philip Ilten <philten@cern.ch>

// Keywords: hadronization; reweighting; tuning

// Demonstrates how to reweight an event for different flavor
// parameters, after the event has been produced, i.e. post-hoc rather
// than in-situ reweighting. The result here should be equivalent to
// the in-situ flavor reweighitng of main263. The power of this method
// is that by saving a minimal set of information per event (10
// doubles per string), the entire event can be reweighted to whatever
// flavor parameters are requested by the user.

// Pythia includes.
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FlavorVariations.h"
#include "Pythia8/PythiaParallel.h"

using namespace Pythia8;

int main() {

  // Define the new set of flavor parameters that we wish to reweight
  // to. This can be set at the event level, but we define it here so
  // that we can compare with the in-situ reweighting.
  double rho = 0.2;  // StringFlav:ProbStoUD
  double xi  = 0.1;  // StringFlav:ProbQQtoQ
  double x   = 0.9;  // StringFlav:ProbSQtoQQ
  double y   = 0.04; // StringFlav:ProbQQ1toQQ0

  // Create, configure, and initialize Pythia.
  PythiaParallel pythia;
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("Beams:eCM = 91.189");
  pythia.readString("PDF:lepton = off");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");
  pythia.readString("HadronLevel:Decay = off");
  pythia.readString("StringFlav:popcornRate = 0");
  //pythia.readString("Parallelism:numThreads = 1");
  pythia.readString("VariationFrag:List = {flav frag:rho=" + to_string(rho)
    + " frag:xi=" + to_string(xi) + " frag:x=" + to_string(x)
    + " frag:y=" + to_string(y) + "}");
  pythia.init();

  // Define the plot title.
  string title = "default: (" +
    toString(pythia.settings.parm("StringFlav:ProbStoUD")) + ", " +
    toString(pythia.settings.parm("StringFlav:ProbQQtoQ")) + ", " +
    toString(pythia.settings.parm("StringFlav:ProbSQtoQQ")) + ", " +
    toString(pythia.settings.parm("StringFlav:ProbQQ1toQQ0")) + "), " +
    "variation: (" + toString(rho) + ", " + toString(xi) + ", " +
    toString(x) + ", " + toString(y) + ") ";

  // Create the reweighting tool given a settings object. Here, since
  // the parallel framework is being used, the settings from the
  // helper Pythia object which is used to intialize the other Pythia
  // versions is passed.
  FlavorVariations vars(pythia.pythiaHelper.settings);
  // Alternatively, the tool can be created by passing the default parameters.
  // FlavorVariations vars(xiDefault, rhoDefault, xDefault, yDefault);

  // Calculate the derived parameters for a given xi, rho, x, and y.
  // This operation is expensive, so should be done up front and only
  // once for each set of primary parameters.
  vector<double> parms = vars.parms(xi, rho, x, y);

  // Identified final state hadrons to include in the histograms.
  vector<int> hadrons = {
    211, 221, 331, 213, 223, 321, 311, 333, 2212, 2112, 2214, 2224, 3222,
    3212, 3122, 3322, 3334};

  // Define multiplicity histograms.
  // default:     the default parameters in Pythia
  // posthoc-wgt: post-hoc reweighted using the Pythia event
  // posthoc-str: post-hoc reweighted using the saved string break
  // insitu:      in-situ reweighted
  // rerun:       a run with the varied parameters
  vector<string> names = {
    "default", "posthoc-wgt", "posthoc-str", "insitu", "rerun"};
  map<string, Hist> hists;
  for (string &name : names)
    hists[name] = Hist(name, hadrons.size(), 0, hadrons.size());

  // Track the weights.
  map<string, double> wgts, sumWgts, sumWgt2s;
  for (string &name : names)
    wgts[name] = sumWgts[name] = sumWgt2s[name] = 0;
  names.pop_back();

  // Run events.
  int nEvent = 1e6;
  // This defines a lambda function that acts as a callback.
  // This function is called for each event generated.
  // The argument is a pointer to the instance that generated the event.
  pythia.run( nEvent,[&](Pythia* pythiaPtr) {

    Event &event = pythiaPtr->event;

    // The necessary information for reweighting later can be saved to
    // a string. Note, other serialization implementations can be
    // used, and could then be implemented with symmetric
    // FlavorVariations::read and FlavorVariations::write methods.
    string saved = vars.write(pythiaPtr->info.weightContainerPtr
      ->weightsFragmentation.flavBreaks);

    // For the default parameters, the weight is just 1.
    wgts["default"] = 1;

    // Calculate the weight for the event, assuming we already have
    // the weight container and associated string breaks.
    wgts["posthoc-wgt"] = vars.weight(parms,
      pythiaPtr->info.weightContainerPtr->weightsFragmentation.flavBreaks);

    // If instead we have saved the breaks to a string, as we did
    // above, we can calculate the weight from the saved string.
    wgts["posthoc-str"] = vars.weight(parms, vars.read(saved));

    // We can also use the in-situ reweighting.
    wgts["insitu"] = pythiaPtr->info.weightValueByIndex(1);

    // Keep track of the weights.
    for (string &name : names) {
      sumWgts[name]  += wgts[name];
      sumWgt2s[name] += pow2(wgts[name]);
    }

    // Fill the histograms.
    for (const Particle &prt : event) {
      if (!prt.isFinal()) continue;
      int pid = prt.idAbs();
      int idx = -1;
      for (int iHad = 0; iHad < (int)hadrons.size(); ++iHad)
        if (pid == hadrons[iHad]) {idx = iHad; break;}
      if (idx < 0) continue;
      for (string &name : names) hists[name].fill(idx, wgts[name]);
    }
  });
  pythia.stat();

  // Rerun Pythia with the varied parameters.
  pythia.settings.parm("StringFlav:ProbStoUD",    rho);
  pythia.settings.parm("StringFlav:ProbQQtoQ",    xi);
  pythia.settings.parm("StringFlav:ProbSQtoQQ",   x);
  pythia.settings.parm("StringFlav:ProbQQ1toQQ0", y);
  pythia.settings.wvec("VariationFrag:List",      {});
  pythia.init();

  pythia.run( nEvent, [&](Pythia* pythiaPtr) {

    Event &event = pythiaPtr->event;
    sumWgts["rerun"]  += 1;
    sumWgt2s["rerun"] += 1;
    for (const Particle &prt : event) {
      if (!prt.isFinal()) continue;
      int pid = prt.idAbs();
      int idx = -1;
      for (int iHad = 0; iHad < (int)hadrons.size(); ++iHad)
        if (pid == hadrons[iHad]) {idx = iHad; break;}
      if (idx >= 0) hists["rerun"].fill(idx, 1.);
    }
  });
  pythia.stat();

  // Normalize the histograms.
  for (auto &hist : hists) hist.second /= sumWgts[hist.first];

  // Print the histogram ratios.
  string xlabel;
  for (int iHad = 0; iHad < (int)hadrons.size(); ++iHad) {
    string name = pythia.particleData.name(hadrons[iHad]);
    cout << left << setw(3) << iHad << ": " << name << "\n";
    xlabel += " " + name + "(" + to_string(iHad) + ")";
  }
  for (auto &hist : hists)
    cout << "\n" << hist.first << hist.second/hists["default"];

  // Print the reweighting stats.
  // The 1 - mu should be statistically consistent with zero if the
  // reweighting has proper coverage.
  // The n_eff gives the statistical power of the reweighted sample.
  for (string &name : names) {
    double w(sumWgts[name]), w2(sumWgt2s[name]), n(sumWgts["default"]);
    cout << name << "\n"
         << "\t1 - mu = " << scientific << setprecision(3) << abs(1. - w/n)
         << " +- "<< abs(1. - sqrt((w2/n - pow2(w/n))*n/(n - 1)))/sqrt(n)
         << "\n\tn_eff  = " << scientific << setprecision(3) << w/sqrt(n*w2)
         << "\n";
  }

  // Create Python plot.
  HistPlot hpl("main264plot");
  hpl.frame("main264plot", title, xlabel, "n(variation)/n(default)");
  for (string &name : names)
    hpl.add(hists[name]/hists["default"], "e", name);
  hpl.add(hists["rerun"]/hists["default"], "e", "rerun");
  hpl.plot();

  return 0;
}
