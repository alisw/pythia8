// main427.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// Keywords: heavy ions; impact parameter; angantyr; event weights; reweighting

// This test program shows how one can replace the impact parameter
// generator in Angantyr, through a UserHook. Specifically, the
// example replaces the existing generator with one that samples
// impact parameters with unit weights. Unit weights can also easier
// be obtained with HeavyIon:forceUnitWeight = on.

#include "Pythia8/Pythia.h"
#include "Pythia8/HIInfo.h"

using namespace Pythia8;

//==========================================================================

// The UnitImpactParameterGenerator generates impact parameters with unit
// weights. It will sample from b = 0 to HeavyIon::bWidth. For realistic
// events, this number should be larger than twice the radius of the
// considered nuclei.

class UnitImpactParameterGenerator : public ImpactParameterGenerator {

public:
  // Constructor.
  UnitImpactParameterGenerator() {}

  // Return the impact parameter.
  Vec4 generate(double& weight) const override {
    // Set unit weight.
    weight = 1.;
    // Generate b.
    double b = sqrt(rndmPtr->flat() * width() * width());
    // Generate a random angle.
    double phi = 2.0*M_PI*rndmPtr->flat();
    // Convert to cartesian coordinates.
    return Vec4(b*sin(phi), b*cos(phi), 0., 0.);
  }

  // Return the cross section scale.
  double xSecScale() const override {
    return M_PI*pow2(width());
  }

};

//==========================================================================

// MyHIUserHooks works as a wrapper around the impact parameter class (above),
// functioning to expose it to Angantyr.

class MyHIUserHooks : public HIUserHooks {

public:
  // Constructor. Create the impact parameter generator.
  MyHIUserHooks() {
    impactGeneratorPtr = make_shared<UnitImpactParameterGenerator>();
  }

  // Tell Angantyr that you want to replace the impact parameter generator.
  bool hasImpactParameterGenerator() const override {
    return true;
  }

  // Return a pointer to the impact parameter generator.
  shared_ptr<ImpactParameterGenerator> impactParameterGenerator()
    const override {
      return impactGeneratorPtr;
    }

private:
  // A shared pointer to the new impact parameter generator lives here as a
  // member. The object is kept alive as long as the user hook is kept alive.
  shared_ptr<UnitImpactParameterGenerator> impactGeneratorPtr;
};

//==========================================================================

// Read common parameters into the two Pythia objects.
void readCommonConfiguration(Pythia& pythia) {
  // Set beams. We choose Pb.
  pythia.readString("Beams:idA = 1000822080");
  pythia.readString("Beams:idB = 1000822080");

  // Collision energy.
  pythia.readString("Beams:eCM = 200");
  pythia.readString("Beams:frameType = 1");
  pythia.readString("SoftQCD:all = on");

  // Initialize the Angantyr model cross sections.
  // Warning: if the model changes, inactivate lines below to retune.
  pythia.readString("HeavyIon:SigFitDefAvNDb = 0.78");
  pythia.readString("HeavyIon:SigFitDefPar = 1.24,3.19,0.72");
  pythia.readString("HeavyIon:SigFitNGen = 0");

}

//==========================================================================

int main() {
  // We make two Pythia objects, one for each impact parameter
  // selector, so we can compare the results.
  Pythia pythiaDefault;
  readCommonConfiguration(pythiaDefault);
  Pythia pythiaCustom;
  readCommonConfiguration(pythiaCustom);
  // Make the user hook.
  shared_ptr<MyHIUserHooks> hiHookPtr = make_shared<MyHIUserHooks>();
  // Make the hook available to Pythia.
  pythiaCustom.setHIHooks(hiHookPtr);
  // Set the width (i.e. maximal b) of the new impact parameter sampler.
  pythiaCustom.readString("HeavyIon:bWidth = 20");

  // Histogram for the impact parameter. One for the default selector.
  Hist impactDefault("bDef", 15, 0., 20.);
  // One for the default selector, not filled with event weights.
  Hist impactDefaultUnw("bDefUnweighted", 15, 0., 20.);
  // One for the custom selector.
  Hist impactCustom("bCustom", 15, 0., 20.);

  // If Pythia fails to initialize, exit with error.
  if (!pythiaDefault.init()) return 1;
  if (!pythiaCustom.init()) return 1;

  // Loop over events.
  int nEvents = 5000;
  double sumWDefault = 0;
  double sumNDefault = 0;
  double sumNCustom = 0;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    // Event with default impact parameter.
    if (pythiaDefault.next()) {
      // Fill the histograms.
      impactDefault.fill(pythiaDefault.info.hiInfo->b(),
        pythiaDefault.info.weight());
      sumWDefault += pythiaDefault.info.weight();
      impactDefaultUnw.fill(pythiaDefault.info.hiInfo->b());
      sumNDefault += 1.0;
    }

    // Event with custom impact parameter.
    if (pythiaCustom.next()) {
      // No need to use event weights to fill this histogram now,
      // as they are set to 1 in the custom class.
      impactCustom.fill(pythiaCustom.info.hiInfo->b());
      sumNCustom += 1.0;
    }
  }

  // Statistics.
  pythiaDefault.stat();
  pythiaCustom.stat();

  // Write histograms.
  cout << "Run is done. Writing histograms to file." << endl;
  HistPlot hpl("plot427");
  hpl.frame("fig427", "Impact parameter distribution" );
  hpl.add(impactDefault/sumWDefault, "-","Default");
  hpl.add(impactDefaultUnw/sumNDefault, "-","Default (unweighted)");
  hpl.add(impactCustom/sumNCustom, "-","Custom (unweighted)");
  hpl.plot(0, 20, 0, 0.2, false);

  return 0;
}
