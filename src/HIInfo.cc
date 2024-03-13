// HIInfo.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HIInfo.h header) for
// the heavy ion classes classes.

#include "Pythia8/Pythia.h"
#include "Pythia8/HIInfo.h"
#include "Pythia8/HINucleusModel.h"
#include "Pythia8/HISubCollisionModel.h"

namespace Pythia8 {

//==========================================================================

// HIInfo functions to collect statistics in an event and in a run.

//--------------------------------------------------------------------------

// Collect statistics for each SubCollision in an event.

int HIInfo::addSubCollision(const SubCollision & c) {
  ++nCollSave[0];
  switch ( c.type ) {
  case SubCollision::ABS:
    return ++nCollSave[1];
  case SubCollision::SDEP:
    return ++nCollSave[2];
  case SubCollision::SDET:
    return ++nCollSave[3];
  case SubCollision::DDE:
    return ++nCollSave[4];
  case SubCollision::CDE:
    return ++nCollSave[5];
  case SubCollision::ELASTIC:
    return ++nCollSave[6];
  default:
    return 0;
  }
}

//--------------------------------------------------------------------------

// Collect statistics for each participating nucleon in an event.

int HIInfo::addProjectileNucleon(const Nucleon & n) {
  ++nProjSave[0];
  switch ( n.status() ) {
  case Nucleon::ABS:
    return ++nProjSave[1];
  case Nucleon::DIFF:
    return ++nProjSave[2];
  case Nucleon::ELASTIC:
    return ++nProjSave[3];
  default:
    return 0;
  }
}

int HIInfo::addTargetNucleon(const Nucleon & n) {
  ++nTargSave[0];
  switch ( n.status() ) {
  case Nucleon::ABS:
    return ++nTargSave[1];
  case Nucleon::DIFF:
    return ++nTargSave[2];
  case Nucleon::ELASTIC:
    return ++nTargSave[3];
  default:
    return 0;
  }
}

//--------------------------------------------------------------------------

// Collect statistics for attemted and accepted impact paramet point
// in an event.

void HIInfo::addAttempt(double T, double bin, double phiin, double bweight) {
  bSave = bin;
  phiSave = phiin;
  nCollSave = nProjSave = nTargSave = vector<int>(10, 0);
  nFailSave = 0;
  weightSave = bweight;
  weightSumSave += weightSave;
  ++NSave;
  double w = 2.0*T*bweight;
  double delta = w - sigmaTotSave;
  sigmaTotSave += delta/double(NSave);
  sigErr2TotSave += (delta*(w - sigmaTotSave) - sigErr2TotSave)/double(NSave);
  w = (2*T - T*T)*bweight;
  delta = w - sigmaNDSave;
  sigmaNDSave += delta/double(NSave);
  sigErr2NDSave += (delta*(w - sigmaNDSave) - sigErr2NDSave)/double(NSave);
}


void HIInfo::accept() {
  int pc = primInfo.code();
  weightSumSave += weightSave;
  ++NAccSave;
  sumPrimW[pc] += weightSave;
  sumPrimW2[pc] += weightSave*weightSave;
  ++NPrim[pc];
  NamePrim[pc] = primInfo.nameProc(pc);
}

//==========================================================================

} // end namespace Pythia8
