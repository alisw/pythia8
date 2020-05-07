// DireMG5MEs.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Stefan Prestel, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for MG5 matrix element wrapper.

#ifndef Pythia8_DireMG5MEs_H
#define Pythia8_DireMG5MEs_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>

// Include Pythia headers.
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Event.h"

// Include MG5 PY8MEs plugin headers.
#ifdef MG5MES
#include "PY8ME.h"
#include "PY8MEs.h"
#endif

namespace Pythia8 {

//==========================================================================

typedef std::vector<double> vec_double;

#ifdef MG5MES
bool isAvailableME(PY8MEs_namespace::PY8MEs& accessor, vector <int> in,
   vector<int> out);
bool isAvailableME(PY8MEs_namespace::PY8MEs& accessor,
   const Pythia8::Event& event);
double calcME(PY8MEs_namespace::PY8MEs& accessor,
   const Pythia8::Event& event);
#else
bool isAvailableME();
double calcME();
#endif

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_DireMG5MEs_H
