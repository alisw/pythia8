// main33.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Philip Ilten <philten@cern.ch>

// Keywords: merging; matching; powheg

// An example where the HVQ POWHEGBOX matrix element binary is
// interfaced directly with PYTHIA. For this example to run correctly
// PYTHIA must be configured with the
// --with-powheg-bin=<path to directory containing only POWHEG binaries>
// option. This will build plugin libraries of the name
// libpythia8powheg<binary name>.so in the library directory.
// For these plugin libraries to build correctly, special compiler flags
// must have been used when building the POWHEGBOX binaries. These are
// "-rdynamic -fPIE -fPIC -pie". The following SED command will correctly
// insert them into the relevant POWHEGBOX Makefile:
//     sed -i "s/F77= gfortran/F77= gfortran -rdynamic -fPIE -fPIC -pie/g"
//     Makefile
// For this specific example the library libpythia8powheghvq.so must
// have been built.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // PYTHIA and the POWHEG user hooks must still be configured, here
  // this is done via main33.cmnd. These settings are sensible
  // defaults, but Powheg:nFinal is dependent upon the POWHEG matrix
  // element being used and so must be changed as appropriate.
  Pythia pythia;
  pythia.readString(
    "Init:plugins = {libpythia8powheghvq.so::LHAupPowheg::main33.cmnd,"
    "libpythia8powhegHooks.so::PowhegHooks::main33.cmnd}");

  // Initialize PYTHIA, based on the specified settings.
  pythia.init();

  // Run PYTHIA. The random numbers are taken from the associated
  // PYTHIA random number generator.
  for (int iEvent = 0; iEvent < 100; ++iEvent) pythia.next();

  // End of run.
  pythia.stat();
  return 0;
}
