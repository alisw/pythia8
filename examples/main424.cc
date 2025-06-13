// main424.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim
//          Leif Lonnblad <leif.lonnblad@fysik.lu.se>

// Keywords: hadron-ion collisions, optimization

// The purpose of this example is to generate an initialization file that
// can be used to speed up initialization in hadron-hadron or hadron-ion runs.
// By default, it produces data for energies from 10 to 10^6 GeV.
// All hadron-nucleon and hadron-ion interactions are possible.
// In other main programs, after initializing with the main424.cmnd file,
// it is possible to change energy and beam types on a per-event basis.
// Warning: the generation takes a long time, of the order of 120 minutes,
// whereof most of it is spent on the nuclear geometry initialization.
// But, of course, once done you will save time in subsequent runs.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  Pythia pythia;
  // Use Angantyr even when initializing with pp.
  pythia.readString("HeavyIon:mode = 2");

  // Variable energy parameters.
  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.readString("Beams:eCM = 1000000");
  pythia.readString("HeavyIon:varECMMin = 10");
  pythia.readString("HeavyIon:varECMMax = 1000000");
  pythia.readString("HeavyIon:varECMSigFitNPts = 11");

  // Variable beam type parameters.
  pythia.readString("Beams:allowIDASwitch = on");

  // Increased statistics for MPI initialization gives smaller errors.
  //pythia.readString("MultipartonInteractions:nSample = 1000000");

  // Specify where to save. If you set reuseInit = 2, the old files
  // will be replaced/overwritten if they already exist. -1 means that
  // the initialization is stored in the settings, and can be written
  // out.
  pythia.readString("MultipartonInteractions:reuseInit = -1");
  pythia.readString("MultipartonInteractions:initFile = ");
  pythia.readString("HeavyIon:SasdMpiReuseInit = -1");
  pythia.readString("HeavyIon:SasdMpiInitFile = ");
  pythia.readString("HeavyIon:SigFitReuseInit = -1");
  pythia.readString("HeavyIon:SigFitInitFile = ");

  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl
         << " Please validate your settings in main424.cmnd." << endl;
    return -1;
  }

  // Note that writeFile will write out all changed settings, not only
  // the ones containing the initialisations, so some manual clean-up
  // is needed to avoid this.
  pythia.readString("HeavyIon:mode = DEFAULT");
  pythia.readString("Next:numberCount = DEFAULT");
  pythia.readString("Next:numberShowEvent = DEFAULT");
  pythia.readString("Next:numberShowInfo = DEFAULT");
  pythia.readString("Next:numberShowLHA = DEFAULT");
  pythia.readString("Next:numberShowProcess = DEFAULT");
  pythia.readString("PartonLevel:all = DEFAULT");
  pythia.readString("ProcessLevel:all = DEFAULT");
  pythia.readString("Tune:ee = DEFAULT");
  pythia.readString("Tune:pp = DEFAULT");
  pythia.settings.writeFile("main424.cmnd");

  // After initializing, we're done.
  return 0;
}
