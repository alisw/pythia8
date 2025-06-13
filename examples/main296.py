# main296.py is a part of the PYTHIA event generator.
# Copyright (C) 2025 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.

# Keywords: Python; total cross section; partial cross sections

# Authors: Christian Bierlich <christian.bierlich@fysik.lu.se>

# This code is an example of how Pythia can be used as a shared library. 
# The physics case is a study of total cross sections, which cannot be
# exposed through the normal Python interface. Instead a class encapsulating
# the desired functionality (PythiaCrossSection) is implemented in 
# main296.cc, and exposed to calls below using ctypes, which comes 
# standard with Python.
# To use this example, main296Lib.cc must first be compiled as a shared 
# library with make libmain296Lib.so, and then the Python code below run. 

# Import ctypes and load the library.
import ctypes
lib = ctypes.cdll.LoadLibrary('./libmain296Lib.so')

# The Python class, PythiaCrossSection, is the Python interface to the
# PythiaCrossSection C++ class defined in main296.cc.
# The class is written as a resource class, with the actual interface only
# accesible inside the class as PythiaInterfacer, to ensure proper garbage 
# collection.
# The PythiaCrossSection class should therefore only be instantiated using 
# a with statement, which will give direct access to the interface, as seen
# in the example use case below.
class PythiaCrossSection:
    def __init__(self, modeSigma):
        self.modeSigma = modeSigma
    def __enter__(self):
        # The inner class, which provides the interface.
        class PythiaInterfacer:
            # Argument types and return types across the interface must be
            # declared as follows.
            def __init__(self, modeSigma):
                lib.PythiaCrossSection_new.argtypes = [ctypes.c_int]
                lib.PythiaCrossSection_new.restype = ctypes.c_void_p

                lib.PythiaCrossSection_del.argtypes = [ctypes.c_void_p]

                lib.PythiaCrossSection_calc.argtypes = [ctypes.c_void_p, 
                ctypes.c_int, ctypes.c_int, ctypes.c_double]

                lib.PythiaCrossSection_sigmaTot.argtypes = [ctypes.c_void_p]
                lib.PythiaCrossSection_sigmaTot.restype = ctypes.c_double

                lib.PythiaCrossSection_sigmaEl.argtypes = [ctypes.c_void_p]
                lib.PythiaCrossSection_sigmaEl.restype = ctypes.c_double

                self.obj = lib.PythiaCrossSection_new(modeSigma)

            # Calculate cross sections (done in the C++ class).
            def calc(self, idA, idB, sqrtS):
                lib.PythiaCrossSection_calc(self.obj, idA, idB, sqrtS)

            # Return the total cross section in mb.
            def sigmaTot(self):
                return lib.PythiaCrossSection_sigmaTot(self.obj)

            # Return the elastic cross section in mb.
            def sigmaEl(self):
                return lib.PythiaCrossSection_sigmaEl(self.obj)

            # Clean up.
            def cleanup(self):
                lib.PythiaCrossSection_del(self.obj)

        self.interfacer = PythiaInterfacer(self.modeSigma)
        return self.interfacer
    # Clean up the interface by deleting the Pythia object.
    def __exit__(self, exc_type, exc_value, traceback):
        self.interfacer.cleanup()

# The following is the example code illustrating a use case of the above.
# Calculate cross section in the range 10-10000 GeV.
sqrts = [10.*i for i in range(1,1000)]
ppTotal = []
ppElastic = []
# mode = 1 is The Schuler and Sjostrand model.
with PythiaCrossSection(1) as pythiaCross:
    for sqs in sqrts:
        pythiaCross.calc(2212,2212,sqs)
        ppTotal.append(pythiaCross.sigmaTot())
        ppElastic.append(pythiaCross.sigmaEl())

# Use matplotlib to directly plot the results.
import matplotlib.pyplot as plt

plt.plot(sqrts, ppTotal,color='red',label=r'$\sigma_{tot}$ (SaS)')
plt.plot(sqrts, ppElastic,color='blue',label=r'$\sigma_{el}$ (SaS)')

plt.xlabel(r'$\sqrt{s}$ [GeV]')
plt.ylabel(r'$\sigma$ [mb]')
plt.xscale('log')
plt.legend()
plt.show()
