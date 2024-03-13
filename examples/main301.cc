// main301.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>
//          Stephen Mrenna <mrenna@fnal.gov>
//          Philip Ilten <philten@cern.ch

// Keywords: hadronization; reweighting; tuning; parallelism; matplotlib

// Demonstrates how to reweight an event for different kinematic or flavor
// hadronization parameters using in-situ reweighting.

// Pythia includes.
#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"

using namespace Pythia8;

int main() {

  // Choose to reweight kinematic (0) or flavor (1) hadronization
  // parameters. Both can be reweighted simultaneously, but the
  // separation is kept here for illustrative purposes.
  int type = 0;

  // Number of events to generate per run.
  int nEvent = 1e6;

  // Define the new set of kinematic parameters that we wish to reweight to.
  double aLund   = 0.6;  // StringZ:aLund, default 0.68
  double bLund   = 0.9;  // StringZ:bLund, default 0.98
  double rFactC  = 1.3;  // StringZ:rFactC, default 1.32
  double rFactB  = 0.9;  // StringZ:rFactB, default 0.855
  double ptSigma = 0.3;  // StringPT:sigma, default 0.335

  // Define the new set of flavor parameters that we wish to reweight to.
  double rho = 0.2;  // StringFlav:ProbStoUD, default 0.217
  double xi  = 0.1;  // StringFlav:ProbQQtoQ, default 0.081
  double x   = 0.9;  // StringFlav:ProbSQtoQQ, default 0.915
  double y   = 0.04; // StringFlav:ProbQQ1toQQ0, default 0.0275

  // Create and configure Pythia.
  PythiaParallel pythia;
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("Beams:eCM = 91.189");
  pythia.readString("PDF:lepton = off");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");

  // Configure the in-situ kinematic or flavor reweighting.
  if (type == 0) {
    pythia.readString("VariationFrag:List = {kineVar frag:aLund="
      + to_string(aLund) + " frag:bLund=" + to_string(bLund) + " frag:rFactC="
      + to_string(rFactC) + " frag:rFactB=" + to_string(rFactB)
      + " frag:ptSigma=" + to_string(ptSigma) + "}");
  } else if (type == 1) {
    pythia.readString("VariationFrag:List = {flavVar frag:rho="
      + to_string(rho) + " frag:xi=" + to_string(xi) + " frag:x="
      + to_string(x) + " frag:y=" + to_string(y) + "}");
    pythia.readString("StringFlav:popcornRate = 0");
    pythia.readString("HadronLevel:Decay = off");
  }

  // Initialize Pythia.
  pythia.init();

  // Define the plot title.
  string title = "default: (";
  if (type == 0)
    title +=
      toString(pythia.settings.parm("StringZ:aLund")) + ", " +
      toString(pythia.settings.parm("StringZ:bLund")) + ", " +
      toString(pythia.settings.parm("StringZ:rFactC")) + ", " +
      toString(pythia.settings.parm("StringZ:rFactB")) + ", " +
      toString(pythia.settings.parm("StringPT:sigma")) + "), " +
      "variation: (" + toString(aLund) + ", " + toString(bLund) + ", " +
      toString(rFactC) + ", " + toString(rFactB) + ", " +
      toString(ptSigma) + ") ";
  else if (type == 1)
    title +=
      toString(pythia.settings.parm("StringFlav:ProbStoUD")) + ", " +
      toString(pythia.settings.parm("StringFlav:ProbQQtoQ")) + ", " +
      toString(pythia.settings.parm("StringFlav:ProbSQtoQQ")) + ", " +
      toString(pythia.settings.parm("StringFlav:ProbQQ1toQQ0")) + "), " +
      "variation: (" + toString(rho) + ", " + toString(xi) + ", " +
      toString(x) + ", " + toString(y) + ") ";

  // Identified final state hadrons to include in the flavor histograms.
  vector<int> hadrons = {
    211, 221, 331, 213, 223, 321, 311, 333, 2212, 2112, 2214, 2224, 3222,
    3212, 3122, 3322, 3334};

  // Define multiplicity histograms. For kinematics, we look at
  // charged multiplicity while for flavor we look at multiplicity per
  // species.
  // default: the default parameters in Pythia
  // insitu:  in-situ reweighted
  // rerun:   a run with the varied parameters
  vector<string> names = {"default", "insitu", "rerun"};
  map<string, Hist> hists;
  for (string &name : names) {
    if (type == 0)
      hists[name] = Hist(name, 25, 2, 51);
    else if (type == 1)
      hists[name] = Hist(name, hadrons.size(), 0, hadrons.size());
  }

  // Track the weights.
  map<string, double> wgts, sumWgts, sumWgt2s;
  for (string &name : names)
    wgts[name] = sumWgts[name] = sumWgt2s[name] = 0;
  names.pop_back();

  // Run events.
  // This defines a lambda function that acts as a callback.
  // This function is called for each event generated.
  // The argument is a pointer to the instance that generated the event.
  // This is neccesary to use PythiaParallel (multi-core).
  pythia.run( nEvent, [&](Pythia* pythiaPtr) {

    // For the default parameters, the weight is just 1.
    wgts["default"] = 1;

    // The weight given by the in-situ reweighting.
    wgts["insitu"] = pythiaPtr->info.getGroupWeight(0);

    // Keep track of the weights.
    for (string &name : names) {
      sumWgts[name]  += wgts[name];
      sumWgt2s[name] += pow2(wgts[name]);
    }

    // Fill the histograms.
    int mult = 0;
    for (const Particle &prt : pythiaPtr->event) {
      if (!prt.isFinal()) continue;
      if (type == 0) {
        if (prt.isCharged()) ++mult;
      } else if (type == 1) {
        int pid = prt.idAbs();
        int idx = -1;
        for (int iHad = 0; iHad < (int)hadrons.size(); ++iHad)
          if (pid == hadrons[iHad]) {idx = iHad; break;}
        if (idx < 0) continue;
        for (string &name : names) hists[name].fill(idx, wgts[name]);
      }
    }
    if (type == 0)
      for (string &name : names) hists[name].fill(mult, wgts[name]);
  });
  pythia.stat();

  // Rerun Pythia with the varied parameters.
  if (type == 0) {
    pythia.settings.parm("StringZ:aLund",  aLund);
    pythia.settings.parm("StringZ:bLund",  bLund);
    pythia.settings.parm("StringZ:rFactC", rFactC);
    pythia.settings.parm("StringZ:rFactB", rFactB);
    pythia.settings.parm("StringPT:sigma", ptSigma);
  } else if (type == 1) {
    pythia.settings.parm("StringFlav:ProbStoUD",    rho);
    pythia.settings.parm("StringFlav:ProbQQtoQ",    xi);
    pythia.settings.parm("StringFlav:ProbSQtoQQ",   x);
    pythia.settings.parm("StringFlav:ProbQQ1toQQ0", y);
  }
  pythia.settings.wvec("VariationFrag:List", {});
  pythia.init();

  // Run events.
  // This defines a lambda function that acts as a callback.
  // This function is called for each event generated.
  // The argument is a pointer to the instance that generated the event.
  // This is neccesary to use PythiaParallel (multi-core).
  pythia.run( nEvent, [&](Pythia* pythiaPtr) {
    sumWgts["rerun"]  += 1;
    sumWgt2s["rerun"] += 1;
    int mult = 0;
    for (const Particle &prt : pythiaPtr->event) {
      if (!prt.isFinal()) continue;
      if (type == 0) {
        if (prt.isCharged()) ++mult;
      } else if (type == 1) {
        int pid = prt.idAbs();
        int idx = -1;
        for (int iHad = 0; iHad < (int)hadrons.size(); ++iHad)
          if (pid == hadrons[iHad]) {idx = iHad; break;}
        if (idx >= 0) hists["rerun"].fill(idx, 1.);
      }
    }
    if (type == 0) hists["rerun"].fill(mult, 1);
  });
  pythia.stat();

  // Normalize the histograms.
  for (auto &hist : hists) hist.second /= sumWgts[hist.first];

  // Print the histogram ratios.
  string xlabel;
  if (type == 0) {
    xlabel = "multiplicity";
  } else if (type == 1) {
    for (int iHad = 0; iHad < (int)hadrons.size(); ++iHad) {
      string name = pythia.particleData.name(hadrons[iHad]);
      cout << left << setw(3) << iHad << ": " << name << "\n";
      xlabel += " " + name + "(" + to_string(iHad) + ")";
    }
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
         << "\n\tn_eff  = " << scientific << setprecision(3) << w*w/(n*w2)
         << "\n";
  }

  // Create the Python plot and return.
  HistPlot hpl("main301plot");
  hpl.frame("main301plot", title, xlabel, "n(variation)/n(default)");
  for (string &name : names)
    hpl.add(hists[name]/hists["default"], "e", name);
  hpl.add(hists["rerun"]/hists["default"], "e", "rerun");
  hpl.plot();
  return 0;

}
