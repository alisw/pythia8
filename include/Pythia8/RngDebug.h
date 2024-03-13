// RngDebug.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file unsets specialized debug versions of random number methods
// from the Rndm class. This header should be included whenever using
// the random number generation debugging, and should follow after the
// inclusion of all PYTHIA headers but before any user headers.

#ifdef RNGDEBUG
#undef flat
#undef xexp
#undef gauss
#undef gamma
#undef phaseSpace2
#endif
