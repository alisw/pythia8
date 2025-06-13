// main162.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Christian T. Preuss <christian.preuss@uni-goettingen.de>

// Keywords: merging; CKKW-L; MESS; UMEPS; NL3; UNLOPS; NLO;

// It illustrates how to do merging, see the Matrix Element
// Merging page in the online manual. An example command is
//     ./main162 -c main162ckkwl.cmnd
// where main162ckkwl.cmnd supplies the commands.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/InputParser.h"

using namespace Pythia8;

//==========================================================================

// Example main program to illustrate multi-jet merging with PYTHIA.

int main(int argc, char** argv) {

  // Set up command line options.
  InputParser ip("Illustrates how to do merging.",
    {"./main162 -c main162ckkwl.cmnd",
        "./main162 -c main162mess.cmnd",
        "./main162 -c main162nl3.cmnd",
        "./main162 -c main162umeps.cmnd",
        "./main162 -c main162unlops.cmnd"});
  ip.require("c", "Use this user-written command file.", {"-cmnd"});

  // Initialize the parser and exit if necessary.
  InputParser::Status status = ip.init(argc, argv);
  if (status != InputParser::Valid) return status;

  // Generator.
  string cmndFile = ip.get<string>("c");
  Pythia pythia;
  pythia.readFile(cmndFile);

  // Number of events.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  // Number of subruns.
  int nSubruns = pythia.mode("Main:numberOfSubruns");

  // Histograms combined over all jet multiplicities.
  int nBins = 100;
  double xMin = 0., xMax = 200.;
  Hist pTWsum("pT of W, summed over all subruns", nBins, xMin, xMax);
  double binWidth = (xMax-xMin)/nBins;

  // Initialise merged cross section and error.
  double sigmaTotal = 0., errorTotal = 0.;
  vector<double> sigmaExcTree, sigmaExcLoop;
  vector<double> sigmaExcTreeSubt, sigmaExcLoopSubt;

  // Start loop over subruns.
  for (int iMerge = 1; iMerge <= nSubruns; ++iMerge) {

    // Read section of command file for this subrun and initialise.
    pythia.readFile(cmndFile, iMerge);

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Determine which merging type this subrun corresponds to.
    bool isTree     = pythia.mergingHooksPtr->doCKKWLMerging()
      || pythia.mergingHooksPtr->doUMEPSTree()
      || pythia.mergingHooksPtr->doNL3Tree()
      || pythia.mergingHooksPtr->doUNLOPSTree();
    bool isLoop     = pythia.mergingHooksPtr->doNL3Loop()
      || pythia.mergingHooksPtr->doUNLOPSLoop();
    bool isTreeSubt = pythia.mergingHooksPtr->doUMEPSSubt()
      || pythia.mergingHooksPtr->doNL3Subt()
      || pythia.mergingHooksPtr->doUNLOPSSubt();
    bool isLoopSubt = pythia.mergingHooksPtr->doUNLOPSSubtNLO();

    // In UNLOPS skip zero-jet tree-level sample if requested.
    if (pythia.mergingHooksPtr->doUNLOPSTree()
      && pythia.mergingHooksPtr->nRequested() == 0) {
      sigmaExcTree.push_back(0.);
      continue;
    }

    // Get the inclusive cross section of this sample.
    double sigmaInc = 0.;
    for (int i = 0; i < pythia.info.nProcessesLHEF(); ++i)
      sigmaInc += pythia.info.sigmaLHEF(i)/MB2PB;

    // Histograms for current jet multiplicity.
    Hist pTWnow("pT of W, current subrun", 100, 0., 200.);

    // Exclusive cross section of this multiplicity.
    double sigmaNow = 0., errorNow = 0.;

    // Start event-generation loop.
    for (int iEvent=0; iEvent<nEvent; ++iEvent) {

      // Generate next event
      if (!pythia.next()) {
        if (pythia.info.atEndOfFile()) break;
        else continue;
      }

      // Get event weight(s).
      double weight        = pythia.info.weight();
      double weightMerging = pythia.info.mergingWeight();
      weight *= weightMerging;
      // Swap sign for counter events (only in UMEPS and UNLOPS).
      if (isTreeSubt || isLoopSubt) weight *= -1.;

      // Do nothing for vanishing event weight.
      if (weight == 0.) continue;

      // Add to current exclusive cross section.
      sigmaNow += weight;
      errorNow += pow2(weight);

      // Find the final copy of the W+, which is after the full shower.
      int iW = 0;
      for (int i = 1; i < pythia.event.size(); ++i)
        if (pythia.event[i].id() == 24) iW = i;
      // Fill the pTW histogram, including merging weight.
      double pTW = pythia.event[iW].pT();
      pTWnow.fill(pTW, weight);

    }

    // Print cross section and errors.
    pythia.stat();

    // Calculate event normalisation, depending on whether events have
    // unit weight (LHA strategy < 4) or come weighted (LHA strategy 4).
    double norm = 1. / pythia.info.nSelected();
    if (abs(pythia.info.lhaStrategy()) != 4) norm *= sigmaInc;

    // Normalise current cross section and add to total cross section.
    sigmaNow   *= norm;
    errorNow   *= pow2(norm);
    sigmaTotal += sigmaNow;
    errorTotal += errorNow;

    // Save sample cross section for output.
    if (isTree) sigmaExcTree.push_back(sigmaNow);
    if (isLoop) sigmaExcLoop.push_back(sigmaNow);
    if (isTreeSubt) sigmaExcTreeSubt.push_back(sigmaNow);
    if (isLoopSubt) sigmaExcLoopSubt.push_back(sigmaNow);

    // Normalise histograms and add to the combined ones.
    pTWnow *= MB2PB * norm / binWidth;
    pTWsum += pTWnow;

  }

  // Print combined pTW histogram to file.
  pTWsum.table("pTWsum.dat");

  // Print cross section information.
  cout << endl << endl;
  cout << " *-------  MEPS Cross Sections  ---------------------*" << endl;
  cout << " |                                                   |" << endl;
  if (sigmaExcTree.size() > 0) {
    cout << " | Exclusive LO cross sections (mb):                 |" << endl;
    for (int i(0); i<int(sigmaExcTree.size()); ++i)
      cout << " |     " << sigmaExcTree.size()-i-1 << "-jet:  "
           << setw(17) << scientific << setprecision(6)
           << sigmaExcTree[i] << "                     |" << endl;
    cout << " |                                                   |" << endl;
  }
  if (sigmaExcLoop.size() > 0) {
    cout << " | Exclusive NLO cross sections (mb):                |" << endl;
    for (int i(0); i<int(sigmaExcLoop.size()); ++i)
      cout << " |     " << sigmaExcLoop.size()-i-1 << "-jet:  "
           << setw(17) << scientific << setprecision(6)
           << sigmaExcLoop[i] << "                     |" << endl;
    cout << " |                                                   |" << endl;
  }
  if (sigmaExcTreeSubt.size() > 0) {
    cout << " | Exclusive LO subtractive cross sections (mb):     |" << endl;
    for (int i(0); i<int(sigmaExcTreeSubt.size()); ++i)
      cout << " |     " << sigmaExcTreeSubt.size()-i << "-jet:  "
           << setw(17) << scientific << setprecision(6)
           << sigmaExcTreeSubt[i] << "                     |" << endl;
    cout << " |                                                   |" << endl;
  }
  if (sigmaExcLoopSubt.size() > 0) {
    cout << " | Exclusive NLO subtractive cross sections (mb):    |" << endl;
    for (int i(0); i<int(sigmaExcLoopSubt.size()); ++i)
      cout << " |     " << sigmaExcLoopSubt.size()-i << "-jet:  "
           << setw(17) << scientific << setprecision(6)
           << sigmaExcLoopSubt[i] << "                     |" << endl;
    cout << " |                                                   |" << endl;
  }
  cout << " |---------------------------------------------------|" << endl;
  cout << " |                                                   |" << endl;
  cout << " | Inclusive merged cross section:                   |" << endl;
  cout << " |                                                   |" << endl;
  cout << " |     " << setw(10) << scientific << setprecision(6)
       << sigmaTotal << " +- " << setw(10) << sqrt(errorTotal) << " mb "
       << "              |" << endl;
  cout << " |                                                   |" << endl;
  cout << " *-------  End MEPS Cross Sections  -----------------*" << endl;
  cout << endl << endl;

  // Done
  return 0;

}
