// main426.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// Keywords: heavy ions; glauber; Npart; Ncoll; fixed target;
// p-Ne; angantyr; total cross sections

// This test program is intended to showcase how initial state
// quantities from the heavy ion model Angantyr can be extracted.
// In the example, Neon is added as a target, Angantyr is run
// in fixed target mode, and cross section and Glauber quantities
// are extracted from the program, as well as calculated directly,
// to illustrate how the numbers are obtained.

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;

//==========================================================================

int main() {

  Pythia pythia;

  // We want to do fixed target, proton on neon. Since 20Ne is not
  // a standard particle in Pythia, we add it to the database before
  // setting the beam.
  pythia.particleData.addParticle(1000100200, "20Ne", 6, 30, 0, 19.992440);

  // Set up beams, use the newly added Ne.
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 1000100200");

  // Run fixed target. When beam energy is lower than the mass per nucleon,
  // it is assumed at rest.
  pythia.readString("Beams:eA = 2500");
  pythia.readString("Beams:eB = 0");
  pythia.readString("Beams:frameType = 2");

  // Parameters for cross section model.
  pythia.readString("Angantyr:CollisionModel = 2");
  pythia.readString("HeavyIon:SigFitDefPar = "
    "21.06,27.80,0.15");
  pythia.readString("HeavyIon:SigFitNGen = 10");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Count the number of inelastic non-diffractive subcollisions.
  // In proton-Ion, this is the same as wounded nucleons in the target.
  // Compare to the ones output by built-in methods.
  Hist nColl1("nColl1", 15, 0, 15);
  Hist nColl2("nColl2", 15, 0, 15);
  Hist nCollGlauber("nCollGlauber", 15, 0, 15);
  // Number of participants in the target.
  Hist nPart("nPart", 15, 0, 15);
  // Sum of weights.
  double sumW = 0;

  // Loop over events.
  for ( int iEvent = 0; iEvent < 100000; ++iEvent ) {
    if ( !pythia.next() ) continue;
    // Short-hand for the pointer to Heavy Ion Info.
    auto hiPtr = pythia.info.hiInfo;
    // Short-hand for the pointer to SubCollisions.
    auto scPtr = hiPtr->subCollisionsPtr();
    // nColl are all the subcollisions which ended up generating
    // particles, i.e. did not fail due to non-conservation of
    // momentum.
    int nColl = 0;
    // nCollGlauber are all subcollisions from the Glauber
    // calculation. Not all of them will be realized.
    int nCollG = 0;
    // Count subcollisions by type. Inelastic non-diffractive are
    // called "absorptive", code ABS.
    for (auto sc : *scPtr) {
      if (sc.type == SubCollision::ABS) ++nCollG;
      if (sc.type == SubCollision::ABS && !sc.failed) ++nColl;
    }
    // The event weight is needed to fill histograms.
    const double weight = pythia.info.weight();
    sumW += weight;
    nColl1.fill(nColl, weight);
    nColl2.fill(hiPtr->nCollND(), weight);
    nCollGlauber.fill(nCollG, weight);
    nPart.fill(hiPtr->nAbsTarg(), weight);
  }

  pythia.stat();
  // Extract the incoherent, inelastic proton-ion cross sections,
  // and print it, i.e. the cross section of everything that produces
  // particles in the final state. It is the sum of all inelastic
  // nucleon-nucleon channels.
  double sigmaInel = 0.;
  double sigmaInel2Err = 0.;
  for (int i : {101, 103, 104, 105}) {
    sigmaInel += pythia.info.sigmaGen(i);
    sigmaInel2Err += pow2(pythia.info.sigmaErr(i));
  }
  // Convenient short-hands.
  const string projName = pythia.particleData.name(pythia.event[1].id());
  const string targName = pythia.particleData.name(pythia.event[2].id());
  const double eCM = pythia.info.eCM();
  cout << "\n=======Output from main426=====================" << endl;
  cout << "Incoherent, inelastic cross section for:" << endl;
  cout << projName << " on " << targName << " at sqrt(s_NN) = ";
  cout << eCM << " GeV." << endl;
  cout << "sigma = " << sigmaInel << " +/- " << sqrt(sigmaInel2Err) <<
  " mb." <<  endl;
  cout << "===============================================" << endl;

  // Plot the histograms with matplotlib.
  HistPlot hpl("plot426");
  hpl.frame("fig426", "Glauber statistics "+projName+" on "+targName+
    " $\\sqrt{s_{\\mathrm{NN}}}$ = "+to_string(eCM)+" GeV", "$N$", "$P(N)$" );
  hpl.add(nColl1/sumW, "-,blue","Subcollisions from main426");
  hpl.add(nColl2/sumW, "x,green","Subcollisions from hiInfo");
  hpl.add(nCollGlauber/sumW, "-,red", "Subcollisions from main426, "
    "including failed");
  hpl.add(nPart/sumW, "*,orange", "Participants from hiInfo");
  hpl.plot(0, 15, 1e-6, 1, true);

  return 0;
}
