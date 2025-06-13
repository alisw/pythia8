// main296Lib.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: Python; total cross section; partial cross sections

// Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

// This code is an example of how Pythia can be used as a shared
// library. The physics case is a study of total cross sections,
// which cannot be exposed through the normal Python
// interface. Instead a class encapsulating the desired functionality
// (PythiaCrossSection) is implemented below, and exposed to calls in
// the corresponding main296.py using ctypes, which comes standard
// with Python. To use this example, the code below must first be
// compiled as a shared library with "make libmain296.so", and then
// the Python code in main296.py run.

#include "Pythia8/Pythia.h"
#include "Pythia8/SigmaTotal.h"

using namespace Pythia8;

// The class exposing cross section calculations.
class PythiaCrossSection{
public:
  // Construct with a specifc model.
  PythiaCrossSection(int modeSig) {
    // Initialize a Pythia object with correct settings.
    pythia.readString("SoftQCD:elastic = on");
    pythia.readString("PartonLevel:all = off");
    pythia.readString("HadronLevel:all = off");
    pythia.settings.mode("SigmaTotal:mode", modeSig);
    pythia.settings.mode("SigmaDiffractive:mode", modeSig);
    pythia.init();

    // Make a copy of the Pythia info object to initialize
    // pointers in SigmaTotal.
    pythiaInfo = pythia.info;

    // Initialize the pointers.
    st.initInfoPtr(pythiaInfo);
    // Initialize the SigmaTotal instance.
    st.init();
  }

  // Calculate cross sections for beam ID selection and sqrt(s).
  void calc(int idA, int idB, double sqrtS){
    st.calc(idA,idB,sqrtS);
  }

  // Output total cross section in mb.
  double sigmaTot() {
    return st.sigmaTot();
  }

  // Output elastic cross section in mb.
  double sigmaEl() {
    return st.sigmaEl();
  }

private:
  Pythia pythia;
  Info pythiaInfo;
  SigmaTotal st;
};

// Since ctypes does not support direct interfacing with C++,
// we go through a C-wrapper.
extern "C" {
    // Return a pointer to an instance of the C++ class.
    PythiaCrossSection* PythiaCrossSection_new(int modeSig){
      return new PythiaCrossSection(modeSig);
    }
    // Garbage collection method.
    void PythiaCrossSection_del(PythiaCrossSection* pxsPtr) {
      delete pxsPtr;
    }
    // Calculate cross sections.
    void PythiaCrossSection_calc(PythiaCrossSection* pxsPtr,
        int idA, int idB, double sqrtS){
      pxsPtr->calc(idA, idB, sqrtS);
    }
    // Return the calculated total cross section.
    double PythiaCrossSection_sigmaTot(PythiaCrossSection* pxsPtr){
      return pxsPtr->sigmaTot();
    }
    // Return the calculated elastic cross section.
    double PythiaCrossSection_sigmaEl(PythiaCrossSection* pxsPtr){
      return pxsPtr->sigmaEl();
    }
}
