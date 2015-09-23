// main63.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Example how you can use UserHooks to enhance rare emission rates,
// in this case q -> q gamma.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Write own derived UserHooks class.

class EnhanceHooks : public UserHooks {

public:  

  // Constructor and destructor do nothing.
  EnhanceHooks() {}
  ~EnhanceHooks() {}

  // Enhance real-emission rate. Thus no trial-emission enhancement. 
  bool canEnhanceEmission() { return true;}
  bool canEnhanceTrial()    { return false;}

  // Function to return the weight enhance factor.
  double enhanceFactor(string name) {
    if (name == "isr:Q2QA") return 50.;
    return 1.0;
  }

  // Function to return the vetoing probability.
  double vetoProbability(string name) {
    if (name == "isr:Q2QA") return 0.5;
    return 0.0;
  }

};

//==========================================================================

int main() {

  // Number of events to generate.
  int nEvent = 1000;

    // Histogram pT spectrum of photons and event weights.
    Hist gamNoEnh(   "gamma pT spectrum, no enhancement",   100, 0., 100.); 
    Hist gamWithEnh( "gamma pT spectrum, with enhancement", 100, 0., 100.); 
    Hist gamBefWt(   "gamma pT spectrum, without weight",   100, 0., 100.); 
    Hist eventWt(   "log10(event weight)",               100, -7., 3.);

  // Compare generation without and with enhanced q -> q gamma emission.
  for (int iCase = 0; iCase < 2; ++iCase) {

    // Generator. Default empty user hook.
    Pythia pythia;
    UserHooks* enhanceHooks = 0;

    //  Process selection. No need to study Z0 decay products or hadron level.
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:mMin = 80.");
    pythia.readString("23:mayDecay = off");
    pythia.readString("HadronLevel:all = off");

    // No event printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    // Set up a user hook and send it in.
    if (iCase == 1) {
      enhanceHooks = new EnhanceHooks();
      pythia.setUserHooksPtr( enhanceHooks);
    }

    // 8 TeV LHC initialization.
    pythia.readString("Beams:eCM = 8000.");
    pythia.init();

    // Begin event loop.
    double sumWt = 0.;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Find and histogram event weight.
      pythia.next();
      double weight = (iCase == 1) 
                    ? enhanceHooks->getEnhancedEventWeight() : 1.;
      if (iCase == 1) eventWt.fill( log10(weight) );
      sumWt += weight;

      // Find all final-state photons and histogram them.
      for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].isFinal() && pythia.event[i].id() == 22) {
        double pT = pythia.event[i].pT();
        if (iCase == 0) gamNoEnh.fill(   pT, 1.);
        if (iCase == 1) gamBefWt.fill(   pT, 1.);
        if (iCase == 1) gamWithEnh.fill( pT, weight);
      }


    // End of event loop.
    }

    // Statistics.
    pythia.stat();
    cout << "\n Average event weight = " << scientific 
         << sumWt / nEvent << endl; 

    // End of case loop. 
    if (iCase == 1) delete enhanceHooks;
  }

  // Histograms and done.
  cout << gamNoEnh << gamWithEnh << gamBefWt << eventWt;
  return 0;
}
