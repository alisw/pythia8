// VinciaQED.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the QED antenna-shower class and auxiliary
// classes. Main author is Rob Verheyen.

#ifndef Pythia8_VinciaQED_H
#define Pythia8_VinciaQED_H

// Standard library.
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cfloat>
#include <cmath>

// PYTHIA 8 headers.
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/PartonSystems.h"

// VINCIA headers.
#include "Pythia8/VinciaCommon.h"
#include "Pythia8/VinciaWeights.h"

namespace Pythia8 {

//==========================================================================

// Class for the "Hungarian" pairing algorithm. Adapted for Vincia
// from an implementation by M. Buehren and C. Ma, see notices below.

// This is a C++ wrapper with slight modification of a hungarian algorithm
// implementation by Markus Buehren. The original implementation is a few
// mex-functions for use in MATLAB, found here:
// http://www.mathworks.com/matlabcentral/fileexchange/
//    6543-functions-for-the-rectangular-assignment-problem
//
// Both this code and the orignal code are published under the BSD license.
// by Cong Ma, 2016.
//
// Copyright (c) 2014, Markus Buehren
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in
// the documentation and/or other materials provided with the distribution
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

class HungarianAlgorithm {

public:

  // Constructor.
  HungarianAlgorithm() {;}
  // Destructor.
  ~HungarianAlgorithm() {;}

  // Function wrapper for solving assignment problem.
  double solve(std::vector <std::vector<double> >& distMatrix,
    std::vector<int>& assignment);

 private:

  // Solve optimal solution for assignment problem using Munkres algorithm,
  // also known as the Hungarian algorithm.
  void optimal(int *assignment, double *cost, double *distMatrix,
    int nOfRows, int nOfColumns);
  // Build the assignment vector.
  void vect(int *assignment, bool *starMatrix, int nOfRows,
    int nOfColumns);
  // Calculate the assignment cost.
  void calcCost(int *assignment, double *cost, double *distMatrix,
    int nOfRows);
  // Factorized step 2a of the algorithm.
  void step2a(int *assignment, double *distMatrix, bool *starMatrix,
    bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns,
    bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
  // Factorized step 2b of the algorithm.
  void step2b(int *assignment, double *distMatrix, bool *starMatrix,
    bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns,
    bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
  // Factorized step 3 of the algorithm.
  void step3(int *assignment, double *distMatrix, bool *starMatrix,
    bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns,
    bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
  // Factorized step 4 of the algorithm.
  void step4(int *assignment, double *distMatrix, bool *starMatrix,
    bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns,
    bool *coveredRows, int nOfRows, int nOfColumns, int minDim,
    int row, int col);
  // Factorized step 5 of the algorithm.
  void step5(int *assignment, double *distMatrix, bool *starMatrix,
    bool *newStarMatrix, bool *primeMatrix, bool *coveredColumns,
    bool *coveredRows, int nOfRows, int nOfColumns, int minDim);
};

//==========================================================================

// Class for QED emissions.

class QEDemitElemental {

public:

  // Friends for internal private members.
  friend class QEDemitSystem;

  // Constuctor.
  QEDemitElemental() {;}

  // Initialize the pointers.
  void initPtr(Rndm* rndmPtrIn, PartonSystems* partonSystemsPtrIn);
  // Initialize.
  void init(Event &event, int xIn, int yIn, double shhIn,
    double verboseIn);
  // Initialize.
  void init(Event &event, int xIn, vector<int> iRecoilIn, double shhIn,
    double verboseIn);
  // Generate a trial point.
  double generateTrial(Event &event, double q2Start, double q2Low,
    double alphaIn, double cIn);

private:

  // Random pointer.
  Rndm* rndmPtr;

  // Parton system pointer.
  PartonSystems* partonSystemsPtr;

  // Trial variables.
  double q2Sav, zetaSav, phiSav;
  double sxjSav, syjSav;
  double alpha, c;
  bool hasTrial;

  // Particle indices.
  int x, y;
  // Recoiler indices.
  vector<int> iRecoil;
  // IDs.
  int idx, idy;
  // Particle masses.
  double mx2, my2;
  // Particle energies.
  double ex, ey;
  // Antenna invariant mass.
  double m2Ant;
  // Antenna dot product.
  double sAnt;
  // The negative of the product of charges.
  double QQ;

  // Type switches.
  bool isII, isIF, isFF, isRF, isIA, isDip;

  // Hadronic invariant mass.
  double shh;

  // Initialization.
  bool isInitPtr, isInit;
  int verbose;

};

//==========================================================================

// Class for a QED emission system.

class QEDemitSystem {

public:

  QEDemitSystem() :
    iSys(-1), shh(-1.), cMat(0.), eleTrial(nullptr), trialIsVec(false),
    beamAPtr(nullptr), beamBPtr(nullptr), infoPtr(nullptr),
    partonSystemsPtr(nullptr), particleDataPtr(nullptr), rndmPtr(nullptr),
    settingsPtr(nullptr), vinComPtr(nullptr), mode(-1), verbose(-1),
    useFullWkernel(false), q2Cut(-1.), isBelowHad(false),
    emitBelowHad(false), isInitPtr(false), isInit(false), TINYPDF(-1.) {;}

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn);
  // Initialize settings for current run.
  void init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, int verboseIn);
  // Prepare a QCD system.
  void prepare(int iSysIn, Event &event, double q2CutIn, bool isBelowHadIn,
    vector<double> evolutionWindowsIn, AlphaEM alIn);

  // Trial antenna function.
  double aTrial(QEDemitElemental* ele, double sxj, double syj, double sxy);
  // Physical antenna function.
  double aPhys (QEDemitElemental* ele, double sxj, double syj, double sxy);
  // Ratio between PDFs.
  double PDFratio(bool isA, double eOld, double eNew, int id, double Qt2);
  // Set up antenna pairing for incoherent mode.
  void buildSystem(Event &event);
  // Generate a trial scale.
  double generateTrialScale(Event &event, double q2Start);
  // Check the veto.
  bool checkVeto(Event &event);
  // Print the QED emit internal system.
  void print();

private:

  // Event system.
  int iSys;
  double shh;

  // Internal storage.
  vector<vector<QEDemitElemental> > eleMat;
  vector<int> iCoh;
  double cMat;
  vector<QEDemitElemental> eleVec;

  // AlphaEM.
  AlphaEM al;

  // Evolution window.
  vector<double> evolutionWindows;

  // Trial pointer.
  QEDemitElemental* eleTrial;
  bool trialIsVec;

  // Pointers.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  Info*          infoPtr;
  PartonSystems* partonSystemsPtr;
  ParticleData*  particleDataPtr;
  Rndm*          rndmPtr;
  Settings*      settingsPtr;
  VinciaCommon*  vinComPtr;

  // Settings.
  int mode, verbose;
  bool useFullWkernel;
  double q2Cut;
  bool isBelowHad;
  bool emitBelowHad;

  // Initialization.
  bool isInitPtr, isInit;

  // PDF check.
  double TINYPDF;

};

//==========================================================================

// Class for trial QED splittings.

class QEDsplitElemental {

public:

  // Friends for internal private members.
  friend class QEDsplitSystem;

  // Default constructor.
  QEDsplitElemental() = default;

  // Constructor.
  QEDsplitElemental(Event &event, int iPhotIn, int iSpecIn):
    iPhot(iPhotIn), iSpec(iSpecIn), ariWeight(0) {
    m2Ant = m2(event[iPhotIn], event[iSpecIn]);
    sAnt = 2.*event[iPhotIn].p()*event[iSpecIn].p();
    m2Spec = event[iSpecIn].m2();}

  // Kallen function.
  double getKallen() {return m2Ant/(m2Ant - m2Spec);}

private:

  // Internal members.
  int iPhot{}, iSpec{};
  double m2Spec{}, m2Ant{}, sAnt{};
  double ariWeight{};
};

//==========================================================================

// Class for a QED splitting system.

class QEDsplitSystem {

public:

  // Constructor.
  QEDsplitSystem() :
    iSys(-1), totIdWeight(-1.), maxIdWeight(-1.), hasTrial(false),
    q2Trial(-1.), zTrial(-1.), phiTrial(-1.), idTrial(0), eleTrial(nullptr),
    nQuark(-1), nLepton(-1), verbose(-1), q2Max(-1.), q2Cut(-1.),
    isBelowHad(false), beamAPtr(nullptr), beamBPtr(nullptr), infoPtr(nullptr),
    partonSystemsPtr(nullptr), particleDataPtr(nullptr), rndmPtr(nullptr),
    settingsPtr(nullptr), vinComPtr(nullptr), isInitPtr(false), isInit(false)
  {;}

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, VinciaCommon* vinComPtrIn);
  // Initialize.
  void init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, int verboseIn);
  // Prepare list of final-state photons - with recoilers - for splittings.
  void prepare(int iSysIn, Event &event, double q2CutIn, bool isBelowHadIn,
    vector<double> evolutionWindowsIn, AlphaEM alIn);
  // Build the splitting system.
  void buildSystem(Event &event);
  // Generate a scale for the system.
  double generateTrialScale(Event &event, double q2Start);
  // Check the veto.
  bool checkVeto(Event &event);
  // Print the system.
  void print();

private:

  // Event system.
  int iSys;

  // AlphaEM.
  AlphaEM al;

  // Evolution window.
  vector<double> evolutionWindows;

  // Weights for splitting IDs.
  vector<int> ids;
  vector<double> idWeights;
  double totIdWeight, maxIdWeight;

  // Antennae.
  vector<QEDsplitElemental> eleVec;

  // Trial variables.
  bool hasTrial;
  double q2Trial, zTrial, phiTrial, idTrial;
  QEDsplitElemental* eleTrial;

  // Settings.
  int nQuark, nLepton, verbose;
  double q2Max, q2Cut;
  bool isBelowHad;

  // Pointers.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  Info*          infoPtr;
  PartonSystems* partonSystemsPtr;
  ParticleData*  particleDataPtr;
  Rndm*          rndmPtr;
  Settings*      settingsPtr;
  VinciaCommon*  vinComPtr;

  // Initialization.
  bool isInitPtr, isInit;

};

//==========================================================================

// Class for a QED conversion system.

class QEDconvSystem {

public:

  // Constructor.
  QEDconvSystem() : totIdWeight(-1.), maxIdWeight(-1.), iSys(-1), shh(-1.),
      s(-1.), iA(-1), iB(-1), isAPhot(false), isBPhot(false), hasTrial(false),
      iPhotTrial(-1), iSpecTrial(-1), q2Trial(-1.), zTrial(-1.), phiTrial(-1.),
      idTrial(-1), nQuark(-1), verbose(-1), q2Cut(-1.),
      isBelowHad(false), beamAPtr(nullptr), beamBPtr(nullptr),
      infoPtr(nullptr), partonSystemsPtr(nullptr), particleDataPtr(nullptr),
      rndmPtr(nullptr), settingsPtr(nullptr), vinComPtr(nullptr),
      isInitPtr(false), isInit(false), TINYPDF(-1.) {;}

  // Initialize the pointers.
  void initPtr(Info* infoPtrIn, VinciaCommon* vinCluPtrIn);
  // Initialize the system.
  void init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, int verboseIn);
  // Prepare for backwards-evolution of photons.
  void prepare(int iSysIn, Event &event, double q2CutIn, bool isBelowHadIn,
    vector<double> evolutionWindowsIn, AlphaEM alIn);
  // Build the system.
  void buildSystem(Event &event);
  // Generate a trial scale.
  double generateTrialScale(Event &event, double q2Start);
  // Check the veto.
  bool checkVeto(Event &event);
  // Print.
  void print();

private:

  // Trial pdf ratios for conversion.
  map<int,double> Rhat;

  // AlphaEM.
  AlphaEM al;

  // Evolution window.
  vector<double> evolutionWindows;

  // Weights for conversion IDs.
  vector<int> ids;
  vector<double> idWeights;
  double totIdWeight, maxIdWeight;
  int iSys;
  double shh;

  // Antenna parameters.
  double s;
  int iA, iB;
  bool isAPhot, isBPhot;

  // Trial variables.
  bool hasTrial;
  int iPhotTrial, iSpecTrial;
  double q2Trial, zTrial, phiTrial, idTrial;

  // Settings.
  int nQuark, verbose;
  double q2Cut;
  bool isBelowHad;

  // Pointers.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  Info*          infoPtr;
  PartonSystems* partonSystemsPtr;
  ParticleData*  particleDataPtr;
  Rndm*          rndmPtr;
  Settings*      settingsPtr;
  VinciaCommon*  vinComPtr;

  // Initialization.
  bool isInitPtr, isInit;
  double TINYPDF;

};

//==========================================================================

// Class for performing QED showers.

class QEDShower {

public:

  // Friends for internal private members.
  friend class VinciaFSR;

  // Constructor.
  QEDShower() : isInitSav(false) {;}

  // Initialise pointers (called at construction time).
  void initPtr(Info* infoPtrIn, VinciaCommon* vinCluPtrIn);
  // Initialise settings for current run (called as part of Pythia::init())
  void init(BeamParticle* beamAPtrIn = 0, BeamParticle* beamBPtrIn = 0);
  // Prepare to shower a system.
  void prepare(int iSysIn, Event &event, bool isBelowHadIn);
  // Update QED shower system(s) each time something has changed in event.
  void update(Event &event, int iSys);
  // Set verbosity level.
  void setVerbose(int verboseIn) {verbose = verboseIn;}
  // Generate a trial scale.
  double generateTrialScale(Event &event, double q2Start);
  // Check the veto.
  bool checkVeto(Event &event);
  // Check if initialized.
  bool isInit() {return isInitSav;}
  // Return the system window.
  int  sysWin() {return iSysTrial;}
  // Return scales.
  double q2minColoured() {return q2minColouredSav;}
  double q2min() {return q2minSav;}

private:

  // Systems.
  vector<int> iSystems;
  vector<QEDemitSystem> emitSystems;
  vector<QEDsplitSystem> splitSystems;
  vector<QEDconvSystem> convSystems;

  // Settings.
  int  verbose;
  bool doQED;
  bool doEmission;
  int  nGammaToLepton, nGammaToQuark;
  bool doConvertGamma, doConvertQuark;

  // Scales.
  double q2minSav, q2minColouredSav;

  // Trial information.
  int    iSysTrial;
  int    iSysIndexTrial;
  double q2Trial;
  bool   isTrialEmit;
  bool   isTrialSplit;
  bool   isTrialConv;

  // Pointers.
  Info*          infoPtr;
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  ParticleData*  particleDataPtr;
  PartonSystems* partonSystemsPtr;
  Rndm*          rndmPtr;
  Settings*      settingsPtr;
  VinciaCommon*  vinComPtr;

  // AlphaEM.
  AlphaEM al;

  // Evolution windows
  vector<double> evolutionWindows;

  // Initialize.
  bool isInitPtr, isInitSav;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_VinciaQED_H
