// HIInfo.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
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

// Collect Glauber statistics for the event, after it
// has been built.
void HIInfo::glauberStatistics() {
  // Remember the nucleons I have already counted for Npart.
  vector<Nucleon*> proj;
  vector<Nucleon*> targ;
  for (auto sc : *subCollisionsPtrSave) {
    if (sc.failed) continue;
    // Count different types of subcollisions
    ++nCollSave[0];
    switch ( sc.type ) {
      case SubCollision::ABS:
        ++nCollSave[1]; break;
      case SubCollision::SDEP:
        ++nCollSave[2]; break;
      case SubCollision::SDET:
        ++nCollSave[3]; break;
      case SubCollision::DDE:
        ++nCollSave[4]; break;
      case SubCollision::CDE:
        ++nCollSave[5]; break;
      case SubCollision::ELASTIC:
        ++nCollSave[6]; break;
      default:
        break;
    }
    // Check if the projectile/target nucleon is already counted.
    bool counted = false;
    for (int i = 0, N = proj.size(); i < N; ++i)
      if (sc.proj == proj[i]) counted = true;

    // If the (projectile) nucleon has not been counted, check how it is
    // participating and add statistics.
    if (!counted) {
      ++nProjSave[0];
      proj.push_back(sc.proj);
      switch ( sc.proj->status() ) {
        case Nucleon::ABS:
          ++nProjSave[1]; break;
        case Nucleon::DIFF:
          ++nProjSave[2]; break;
        case Nucleon::ELASTIC:
          ++nProjSave[3]; break;
        default:
          break;
      }
    }

    // Repeat for target side.
    counted = false;
    for (int i = 0, N = targ.size(); i < N; ++i)
      if (sc.targ == targ[i]) counted = true;
    if (!counted) {
      ++nTargSave[0];
      targ.push_back(sc.targ);
      switch ( sc.targ->status() ) {
        case Nucleon::ABS:
          ++nTargSave[1]; break;
        case Nucleon::DIFF:
          ++nTargSave[2]; break;
        case Nucleon::ELASTIC:
          ++nTargSave[3]; break;
        default:
          break;
      }
    }
  }
}

//--------------------------------------------------------------------------

// Collect statistics for attempted and accepted impact parameter point
// in an event.

void HIInfo::addAttempt(double T, double bin, double phiin, double bweight,
                        double xSecScale) {
  bSave = bin;
  phiSave = phiin;
  nCollSave = nProjSave = nTargSave = vector<int>(10, 0);
  nFailSave = 0;
  weightSave = bweight;
  weightSumSave += weightSave;
  xSecScaleSave = xSecScale;
  double xweight = bweight*xSecScale;
  ++NSave;
  TSave = T;
  double T11 = T;
  double T12 = subCollisionsPtrSave->T(1);
  double T21 = subCollisionsPtrSave->T(2);
  double T22 = subCollisionsPtrSave->T(3);
  double wtot = 0.5*(T11 + T12 + T21 + T22)*xweight;
  runAvg(sigmaTotSave, sigErr2TotSave, NSave, wtot);
  double wdif = 0.25*(T11*T11 + T12*T12 + T21*T21 + T22*T22)*xweight;
  runAvg(sigmaNDSave, sigErr2NDSave, NSave, wtot - wdif);
  double wel = 0.5*(T11*T22 + T12*T21)*xweight;
  runAvg(sigmaELSave, sigErr2ELSave, NSave, wel);
  double wdift = 0.5*(T11*T21 + T22*T12)*xweight;
  runAvg(sigmaDiffTSave, sigErr2DiffTSave, NSave, wdift - wel);
  double wdifp = 0.5*(T11*T12 + T22*T21)*xweight;
  runAvg(sigmaDiffPSave, sigErr2DiffPSave, NSave, wdifp - wel);
  runAvg(sigmaDDiffSave, sigErr2DDiffSave, NSave, wdif - wdift - wdifp + wel);
  runAvg(sigmaINELSave, sigErr2INELSave, NSave, wtot - wel);
  runAvg(slopeSave, slopeErrSave, NSave, pow2(bSave)*wtot/2.0);

}

//--------------------------------------------------------------------------

// Reset the Glauber statistics.

void HIInfo::glauberReset() {
  sigmaTotSave = sigmaNDSave = sigmaELSave = sigmaINELSave = sigmaDiffPSave =
    sigmaDiffTSave = sigmaDDiffSave = slopeSave = 0.0;
  sigErr2TotSave = sigErr2NDSave = sigErr2ELSave = sigErr2INELSave =
    sigErr2DiffPSave = sigErr2DiffTSave = sigErr2DDiffSave =
    slopeErrSave = 0.0;
  NSave = NAccSave = 0;
}

//--------------------------------------------------------------------------

// Indicate that the last generated collision system was accepted.

void HIInfo::accept() {
  int pc = primInfo.code();
  ++NAccSave;
  sumPrimW[pc] += weightSave*xSecScaleSave;
  sumPrimW2[pc] += pow2(weightSave*xSecScaleSave);
  ++NPrim[pc];
  NamePrim[pc] = primInfo.nameProc(pc);
}

//==========================================================================

} // end namespace Pythia8
