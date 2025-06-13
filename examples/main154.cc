// main154.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Philip Ilten <philten@cern.ch>

// Keywords: merging; matching; powheg

// An example where the HVQ POWHEGBOX matrix element binary is
// interfaced directly with PYTHIA. For this example to run correctly
// PYTHIA must be configured with the
// --with-powheg=<path to directory containing only POWHEG libraries>
// option. This will build plugin libraries of the name
// libpythia8powheg<process name>.so in the library directory.
// For these plugin libraries to build correctly, special compiler flags
// must have been used when building the POWHEGBOX binaries. These are
// "-rdynamic -fPIC". The following SED command will correctly
// insert them into the relevant POWHEGBOX Makefile:
//     sed -i "s/F77= gfortran/F77= gfortran -rdynamic -fPIC/g"
//     Makefile
// For this specific example the library libpythia8powheghvq.so must
// have been built. A complete set of libraries can be built using the LCG
// script available here:
//     https://gitlab.cern.ch/sft/lcgcmake/-/blob/master/generators/
//     powheg-box-v2/prepare.sh

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // PYTHIA and the POWHEG user hooks must still be configured, here
  // this is done via main154.cmnd. These settings are sensible
  // defaults, but Powheg:nFinal is dependent upon the POWHEG matrix
  // element being used and so must be changed as appropriate.
  Pythia pythia;
  pythia.readString(
    "Init:plugins = {libpythia8powheghvq.so::LHAupPowheg::set::main154.cmnd,"
    "libpythia8powhegHooks.so::PowhegHooks::set::main154.cmnd}");

  // Initialize PYTHIA, based on the specified settings.
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Run PYTHIA. The random numbers are taken from the associated
  // PYTHIA random number generator.
  for (int iEvent = 0; iEvent < 100; ++iEvent) pythia.next();

  // End of run.
  pythia.stat();
  return 0;
}
