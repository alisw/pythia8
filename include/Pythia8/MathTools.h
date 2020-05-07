// MathTools.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for some mathematics tools, like special functions.
#ifndef Pythia8_MathTools_H
#define Pythia8_MathTools_H

// Header file for the MathTools methods.
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// The Gamma function for real argument.
double GammaReal(double x);

// Modified Bessel functions of the first and second kinds.
double besselI0(double x);
double besselI1(double x);
double besselK0(double x);
double besselK1(double x);

// Integrate f(x) dx over the specified range
bool integrateGauss(double& resultOut, function<double(double)> f,
  double xLo, double xHi, double tol=1e-6);

// Solve f(x) = target for x in the specified range
bool brent(double& solutionOut, function<double(double)> f,
  double target, double xLo, double xHi, double tol=1e-6, int maxIter = 10000);

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_MathTools_H
