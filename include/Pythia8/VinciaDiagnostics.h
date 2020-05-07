// VinciaDiagnostics.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

//File Created 13/11/18 by H Brooks

#ifndef VINCIA_DIAG_H
#define VINCIA_DIAG_H

#include "Pythia8/UserHooks.h"

//completely abstract class for user-defined diagnostics
//to be used like UserHooks - but for diagnostic purposes only
//all functions void and arguments are passed by value or
//by const reference (like const Event &)

namespace Pythia8 {

class VinciaDiagnostics : public UserHooks {

 public:

  //Default constructor
  VinciaDiagnostics(){};

  //Default destructor
  ~VinciaDiagnostics(){};

  virtual void init() = 0;

  virtual void setBranchType(int branchType) = 0;

  virtual void setnBranchSys(int iSys, int nBranch) = 0;

  virtual void checkInvariants(int iSys,int iant, vector<double> invariants,
    bool inPHSP) = 0;

  virtual void checkAnt(int iSys, double ant) = 0;

  virtual void checkAntHel(int iSys, double ant, vector<int> helsIn,
    vector<int> helsOut) = 0;

  virtual void checkpAccept(int iSys, double pAccept) = 0;

  virtual void checkEvent(int iSys,const Event &event,int sizeOld) = 0;

};

} // End namespace Pythia8

#endif
