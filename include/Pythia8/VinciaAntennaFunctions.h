// VinciaAntennaFunctions.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains header information for the AntennaFunction base
// class, its derived classes for FF, IF, and II antenna functions,
// the AntennaSetFSR and AntennaSetISR classes, and the MEC class for
// matrix- element corrections.

#ifndef Pythia8_VinciaAntennaFunctions_H
#define Pythia8_VinciaAntennaFunctions_H

// Pythia headers.
#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/PythiaStdlib.h"

// Vincia headers.
#include "Pythia8/VinciaCommon.h"
#include "Pythia8/VinciaMG5MEs.h"

namespace Pythia8 {

//==========================================================================

// A class containing DGLAP splitting functions for limit checking.

class DGLAP {

public:

  // Constructor.
  DGLAP() {;}

  // Helicity-dependent Altarelli-Parisi kernels (mu = m/Q).  Note,
  // Pg2gg is written with standard factor 2 normalization convention.
  double Pg2gg(double z, int hA=9, int hB=9, int hC=9);
  double Pg2qq(double z, int hA=9, int hB=9, int hC=9, double mu=0.);
  double Pq2qg(double z, int hA=9, int hB=9, int hC=9, double mu=0.);
  double Pq2gq(double z, int hA=9, int hB=9, int hC=9, double mu=0.);

  // Wrappers to get unpolarized massive Altarelli-Pariis kernels (mu = m/Q).
  double Pg2qq(double z, double mu) {return Pg2qq(z, 9, 9, 9, mu);}
  double Pq2qg(double z, double mu) {return Pq2qg(z, 9, 9, 9, mu);}
  double Pq2gq(double z, double mu) {return Pq2gq(z, 9, 9, 9, mu);}

  // Altarelli-Parisi kernels with linear in/out polarizations for
  // gluons: pol = +1 for in-plane, -1 for out-of-plane.
  double Pg2ggLin(double z, int polA = 9, int polB = 9, int polC = 9);
  double Pg2qqLin(double z, int polA = 9, int hB = 9, int hC = 9,
                  double mu = 0.);
  double Pq2qgLin(double z, int hA = 9, int hB=9, int polC = 9,
                  double mu = 0.);
  double Pq2gqLin(double z, int hA = 9, int polB=9, int hC = 9,
                  double mu = 0.);

};

//==========================================================================

// The AntennaFunction base class. Base class implementation for all
// AntennaFunction objects.

class AntennaFunction {

public:

  // Constructor.
  AntennaFunction() = default;

  // Destructor
  virtual ~AntennaFunction() {};

  // Names of this antenna, for VINCIA, and for humans.
  virtual string vinciaName() const = 0;

  // Parton types (AB -> 0i 1j 2k): needed by soft- and collinear-limit checks.
  virtual int idA() const = 0;
  virtual int idB() const = 0;
  virtual int id1() const = 0;

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> mNew,
    vector<int> helBef, vector<int> helNew) = 0;

  // Optional implementation of the DGLAP kernels for collinear-limit checks
  // Defined as PI/sij + PK/sjk, i.e. equivalent to antennae.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew) = 0;

  // Default initialization.
  virtual bool init();

  // Construct baseName from idA, idB, and id1.
  virtual string baseName() const {
    return id2str(id1()) + "/" + id2str(idA()) + id2str(idB());}

  // Wrapper that can modify baseName to more human readable form if required.
  virtual string humanName() const {return baseName();}

  // Function to check singularities, positivity, etc.
  virtual bool check();

  // Method to intialise mass values.
  virtual void initMasses(vector<double>* masses) {
    if (masses->size() >= 3) {
      mi = masses->at(0); mj = masses->at(1); mk = masses->at(2);
    } else {mi = 0.0; mj = 0.0; mk = 0.0;}}

  // Method to initialise internal helicity variables.
  virtual int initHel(vector<int>* helBef, vector<int>* helNew);

  // Wrapper for helicity-summed/averaged antenna function.
  double antFun(vector<double> invariants, vector<double> masses) {
    return antFun(invariants, masses, hDum, hDum);}

  // Wrapper for massless, helicity-summed/averaged antenna function.
  double antFun(vector<double> invariants) {
    return antFun(invariants, mDum, hDum, hDum);}

  // Wrapper without helicity assignments.
  double AltarelliParisi(vector<double> invariants, vector<double> masses) {
    return AltarelliParisi(invariants, masses, hDum, hDum);}

  // Wrapper for massless helicity-summed/averaged DGLAP kernels.
  double AltarelliParisi(vector<double> invariants) {
    return AltarelliParisi(invariants, mDum, hDum, hDum);}

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, DGLAP* dglapPtrIn);

  // Get parameters.
  double chargeFac()  {return chargeFacSav;}
  int    kineMap()    {return kineMapSav;}
  double alpha()      {return alphaSav;}
  double sectorDamp() {return sectorDampSav;}

  // Functions to get Altarelli-Parisi energy fractions from invariants.
  double zA(vector<double> invariants) {
    double yij = invariants[1]/invariants[0];
    double yjk = invariants[2]/invariants[0];
    return (1.-yjk)/(1.+yij);}
  double zB(vector<double> invariants) {
    double yij = invariants[1]/invariants[0];
    double yjk = invariants[2]/invariants[0];
    return (1.-yij)/(1.+yjk);}

  // Auxiliary function to translate an ID code to a string.
  string id2str(int id) const;

protected:

  // Is initialized.
  bool isInitPtr{false}, isInit{false};

  // Charge factor, kinematics map, and subleading-colour treatment.
  double chargeFacSav{0.0};
  int kineMapSav{0}, modeSLC{-1};
  bool sectorShower{false};

  // The alpha collinear-partitioning parameter.
  double alphaSav{0.0};

  // The sector-shower collinear dampening parameter.
  double sectorDampSav{0.0};

  // Shorthand for commonly used variable(s).
  double term{}, preFacFiniteTermSav{0.0}, antMinSav{0.0};
  bool   isMinVar{};

  // Variables for internal storage of masses and helicities.
  double mi{0.0}, mj{0.0}, mk{0.0};
  int hA{9}, hB{9}, hi{9}, hj{9}, hk{9};

  // Map to tell whether a given helicity value maps to L- and/or
  // R-handed. Defined by constructor and not to be changed
  // dynamically.
  map<int, bool> LH{{9, true}, {1, false}, {-1, true}};
  map<int, bool> RH{{9, true}, {1, true},  {-1, false}};

  // Verbosity level.
  int verbose{1};

  // Pointers to Pythia8 classes.
  Info*         infoPtr{};
  ParticleData* particleDataPtr{};
  Settings*     settingsPtr{};
  Rndm*         rndmPtr{};

  // Pointer to VINCIA DGLAP class.
  DGLAP* dglapPtr{};

  // Dummy vectors.
  vector<double> mDum{0, 0, 0, 0};
  vector<int> hDum{9, 9, 9, 9};

};

//==========================================================================

// Class QQEmitFF, final-final antenna function.

class QQEmitFF : public AntennaFunction {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:QQEmitFF";};

  // Functions needed by soft- and collinear-limit checks (AB -> 0i 1j 2k).
  virtual int idA() const {return 1;}
  virtual int idB() const {return -1;}
  virtual int id1() const {return 21;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> mNew,
    vector<int> helBef, vector<int> helNew);

  // Function to give Altarelli-Parisi limits of this antenna.
  // Defined as PI/sij + PK/sjk, i.e. equivalent to antennae.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class QGEmitFF, final-final antenna function.

class QGEmitFF : public AntennaFunction {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:QGEmitFF";};

  // Parton types (AB -> 0i 1j 2k): needed by soft- and collinear-limit checks.
  virtual int idA() const {return 1;}
  virtual int idB() const {return 21;}
  virtual int id1() const {return 21;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

  // Function to give Altarelli-Parisi limits of this antenna.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double> /* mNew */, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GQEmitFF, final-final antenna function.

class GQEmitFF : public QGEmitFF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GQEmitFF";};

  // Parton types (AB -> 0i 1j 2k): needed by soft- and collinear-limit checks.
  virtual int idA() const {return 21;}
  virtual int idB() const {return -1;}
  virtual int id1() const {return 21;}

  // The antenna function [GeV^-2] (derived from QGEmit by swapping).
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

  // Function to give Altarelli-Parisi limits of this antenna.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GQEmitFF, final-final antenna function.

class GGEmitFF : public AntennaFunction {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName()    const {return "Vincia:GGEmitFF";};

  // Parton types (AB -> 0i 1j 2k): needed by soft- and collinear-limit checks.
  virtual int idA()    const {return 21;}
  virtual int idB()    const {return 21;}
  virtual int id1()    const {return 21;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

  // Function to give Altarelli-Parisi limits of this antenna.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GXSplitFF, final-final antenna function.

class GXSplitFF : public AntennaFunction {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GXSplitFF";};

  // Parton types (AB -> 0i 1j 2k): needed by soft- and collinear-limit checks.
  virtual int idA() const {return 21;}
  virtual int idB() const {return  0;}
  virtual int id1() const {return -1;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

  // Function to give Altarelli-Parisi limits of this antenna.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class QQEmitFFsec, sector final-final antenna function, identical
// to global one.

class QQEmitFFsec : public QQEmitFF {};

//==========================================================================

// Class QGEmitFFsec, sector final-final antenna function, explicit
// symmetrisation of QGEmitFF.

class QGEmitFFsec : public QGEmitFF {

public:

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GQEmitFFsec, sector final-final antenna function, explicit
// symmetrisation of GQEmitFF.

class GQEmitFFsec : public QGEmitFFsec {

public:

  // Parton types (AB -> 0i 1j 2k): needed by soft- and collinear-limit checks.
  virtual int idA() const {return 21;}
  virtual int idB() const {return -1;}
  virtual int id1() const {return 21;}

  // The antenna function [GeV^-2] (derived from QGEmitFFsec by swapping).
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

  // Function to give Altarelli-Parisi limits of this antenna.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GGEmitFFsec, sector final-final antenna function, explicit
// symmetrisation of GGEmitFF.

class GGEmitFFsec : public GGEmitFF {

public:

  // The dimensionless antenna function.
  virtual double antFun(vector<double> invariants, vector<double> mNew,
    vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GXSplitFFsec, sector final-final antenna function, explicit
// symmetrisation of GXSplitFF.

class GXSplitFFsec : public GXSplitFF {

 public:

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> mNew,
    vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class AntennaFunctionIX, base class for initial-initial and
// initial-final antenna functions which implements II. All derived
// classes for initial-initial antenna functions are the same for
// global and sector cases since even the global ones include sector
// terms representing "emission into the initial state".

class AntennaFunctionIX : public AntennaFunction {

public:

  // Method to initialise (can be different than that of the base class).
  virtual bool init();

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:AntennaFunctionIX";}

  // Parton types AB -> 0a 1j 2b with A,B,a,b initial and j final.
  virtual int idA() const {return 0;}
  virtual int idB() const {return 0;}
  virtual int id0() const {return 0;}
  virtual int id1() const {return 0;}
  virtual int id2() const {return 0;}

  // Functions to get Altarelli-Parisi energy fractions.
  virtual double zA(vector<double> invariants) {double sAB = invariants[0];
    double sjb = invariants[2]; return sAB/(sAB+sjb);}
  virtual double zB(vector<double> invariants) {double sAB = invariants[0];
    double saj = invariants[1]; return sAB/(sAB+saj);}

  // Function to tell if this is an II antenna.
  virtual bool isIIant() {return true;}

  // Function to check singularities, positivity, etc.
  virtual bool check();

};

//==========================================================================

// Class QQEmitII, initial-initial antenna function.

class QQEmitII : public AntennaFunctionIX {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:QQEmitII";}

  // Parton types AB -> 0a 1j 2b with A,B,a,b initial and j final.
  virtual int idA() const {return 1;}
  virtual int idB() const {return -1;}
  virtual int id0() const {return 1;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return -1;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // AP splitting kernel for collinear limit checks.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GQEmitII, initial-initial antenna function.

class GQEmitII : public AntennaFunctionIX {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GQEmitII";}

  // Parton types AB -> 0a 1j 2b with A,B,a,b initial and j final.
  virtual int idA() const {return 21;}
  virtual int idB() const {return 1;}
  virtual int id0() const {return 21;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return 1;}

  // The antenna function.
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // AP splitting kernel for collinear limit checks.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GGEmitII, initial-initial antenna function.

class GGEmitII : public AntennaFunctionIX {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GGEmitII";}

  // Parton types AB -> 0a 1j 2b with A,B,a,b initial and j final.
  virtual int idA() const {return 21;}
  virtual int idB() const {return 21;}
  virtual int id0() const {return 21;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return 21;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // AP splitting kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class QXSplitII, initial-initial antenna function. Splitting is in
// the forwards sense, i.e. quark backwards evolving to a gluon and
// emitting an antiquark in the final state.

class QXSplitII : public AntennaFunctionIX {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const { return "Vincia:QXSplitII";}

  // Parton types AB -> 0a 1j 2b with A,B, a,b initial and j final.
  virtual int idA() const {return 1;}
  virtual int idB() const {return 0;}
  virtual int id0() const {return 21;}
  virtual int id1() const {return -1;}
  virtual int id2() const {return 0;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // AP splitting kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

  // Mark that this function has no zB collinear limit.
  virtual double zB(vector<double>) {return -1.0;}

};

//==========================================================================

// Class GXConvII, initial-initial antenna function. Gluon evolves
// backwards into a quark and emits a quark in the final state.

class GXConvII : public AntennaFunctionIX {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GXConvII";}

  // Parton types AB -> 0a 1j 2b with A,B,a,b initial and j final.
  virtual int idA() const {return 21;}
  virtual int idB() const {return 0;}
  virtual int id0() const {return 2;}
  virtual int id1() const {return 2;}
  virtual int id2() const {return 0;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // AP splitting kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

  // Mark that this function has no zB collinear limit.
  virtual double zB(vector<double>) {return -1.0;}

};

//==========================================================================

// Class AntennaFunctionIF, base class for IF/RF antenna functions
// which implements QQEmitIF. Derived classes are for global
// initial-final and resonance-final antenna functions. The method
// isRFant() distinguishes between the two.

class AntennaFunctionIF : public AntennaFunctionIX {

public:

  // Method to initialise (can be different than that of the base class).
  virtual bool init();

  // Function to check singularities, positivity, etc.
  virtual bool check();

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:AntennaFunctionIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 1;}
  virtual int idB() const {return -1;}
  virtual int id0() const {return 1;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return -1;}

  // Functions to get Altarelli-Parisi energy fractions.
  virtual double zA(vector<double> invariants) {double sAK = invariants[0];
    double sjk = invariants[2]; return sAK/(sAK+sjk);}
  virtual double zB(vector<double> invariants) {double sAK = invariants[0];
    double saj = invariants[1]; return (sAK-saj)/sAK;}

  // Methods to tell II, IF, and RF apart.
  virtual bool isIIant() {return false;}
  virtual bool isRFant() {return false;}

  // Check for resonances.
  virtual bool checkRes();

  // Create the test masses for the checkRes method.
  virtual void getTestMasses(vector<double> &masses) {masses.resize(4, 0.0);}

  // Create the test invariants for the checkRes method.
  virtual bool getTestInvariants(vector<double> &invariants,
    vector<double> masses, double yaj, double yjk);

protected:

  // Massive eikonal factor, n.b. is positive definite - proportional
  // to gram determinant.
  double massiveEikonal(double saj, double sjk, double sak,
    double m_a, double m_k) {return 2.0*sak/(saj*sjk) - 2.0*m_a*m_a/(saj*saj)
      - 2.0*m_k*m_k/(sjk*sjk);}

  // Massive eikonal factor, given invariants and masses.
  double massiveEikonal(vector<double> invariants, vector<double> masses) {
    return massiveEikonal(invariants[1], invariants[2], invariants[3],
                          masses[0], masses[2]);}

  // Return the Gram determinant.
  double gramDet(vector<double> invariants, vector<double> masses) {
    double saj(invariants[1]), sjk(invariants[2]), sak(invariants[3]),
      mares(masses[0]), mjres(masses[1]), mkres(masses[2]);
    return 0.25*(saj*sjk*sak - saj*saj*mkres*mkres -sak*sak*mjres*mjres
      - sjk*sjk*mares*mares + 4.0*mares*mares*mjres*mjres*mkres*mkres);}

  // Wrapper for comparing to AP functions, sums over flipped
  // invariants where appropriate.
  double antFunCollLimit(vector<double> invariants,vector<double> masses);

};

//==========================================================================

// Class QQEmitIF, initial-final antenna function.

class QQEmitIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const { return "Vincia:QQEmitIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 1;}
  virtual int idB() const {return -1;}
  virtual int id0() const {return 1;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return -1;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // The AP kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

  // Functions to get Altarelli-Parisi energy fractions.
  virtual double zA(vector<double> invariants) {double sAK = invariants[0];
    double sjk = invariants[2]; return sAK/(sAK+sjk);}
  virtual double zB(vector<double> invariants) {double sAK = invariants[0];
    double saj = invariants[1]; return (sAK-saj)/sAK;}

};

//==========================================================================

// Class QGEmitIF, initial-final antenna function.

class QGEmitIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:QGEmitIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 1;}
  virtual int idB() const {return 21;}
  virtual int id0() const {return 1;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return 21;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // The AP kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GQEmitIF, initial-final antenna function.

class GQEmitIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GQEmitIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 21;}
  virtual int idB() const {return 1;}
  virtual int id0() const {return 21;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return 1;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // The AP kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GGEmitIF, initial-final antenna function.

class GGEmitIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GGEmitIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 21;}
  virtual int idB() const {return 21;}
  virtual int id0() const {return 21;}
  virtual int id1() const {return 21;}
  virtual int id2() const {return 21;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // The AP kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class QXSplitIF, initial-final antenna function. Splitting is in
// the forwards sense, i.e. quark backwards evolving to a gluon and
// emitting an antiquark in the final state.

class QXSplitIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:QXSplitIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 1;}
  virtual int idB() const {return 0;}
  virtual int id0() const {return 21;}
  virtual int id1() const {return -1;}
  virtual int id2() const {return 0;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  virtual double AltarelliParisi(vector<double> invariants,
    vector<double> /* mNew */, vector<int> helBef, vector<int> helNew);

  // Mark that this function does not have a zB collinear limit.
  virtual double zB(vector<double>) {return -1.0;}

};

//==========================================================================

// Class GXConvIF, initial-final antenna function. Gluon evolves
// backwards into a quark and emits a quark in the final state.

class GXConvIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:GXConvIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 21;}
  virtual int idB() const {return 0;}
  virtual int id0() const {return 2;}
  virtual int id1() const {return 2;}
  virtual int id2() const {return 0;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // The AP kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

  // Mark that this function does not have a zB collinear limit.
  virtual double zB(vector<double>) {return -1.0;}

};

//==========================================================================

// Class XGSplitIF, initial-final antenna function. Gluon splitting in
// the final state.

class XGSplitIF : public AntennaFunctionIF {

public:

  // Names (remember to redefine both for each inherited class).
  virtual string vinciaName() const {return "Vincia:XGsplitIF";}

  // Parton types AB -> 0a 1j 2b with A,a initial and B,b,j final.
  virtual int idA() const {return 0;}
  virtual int idB() const {return 21;}
  virtual int id0() const {return 0;}
  virtual int id1() const {return -1;}
  virtual int id2() const {return 1;}

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> masses,
    vector<int> helBef, vector<int> helNew);

  // The AP kernel, P(z)/Q2.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double>, vector<int> helBef, vector<int> helNew);

  // Mark that this function does not have a zA collinear limit.
  virtual double zA(vector<double>) {return -1.0;}

};

//==========================================================================

// Class QGEmitIFsec, derived class for sector initial-final antenna
// function. Note only the final-state leg needs to be symmetrised,
// as the global IF functions already contain sector terms on their
// initial-state legs to account for the absence of "emission into the
// initial state".

class QGEmitIFsec : public QGEmitIF {

public:

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class GGEmitIFsec, sector initial-final antenna function.

class GGEmitIFsec : public GGEmitIF {

public:

  // The antenna function [GeV^-2].
  virtual double antFun(vector<double> invariants,
    vector<double> mNew, vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class XGSplitIFsec, sector initial-final antenna function. Gluon
// splitting in the final state.

class XGSplitIFsec : public XGSplitIF {

public:

  // The antenna function, just 2*global [GeV^-2].
  virtual double antFun(vector<double> invariants, vector<double> mNew,
    vector<int> helBef, vector<int> helNew);

};

//==========================================================================

// Class QQEmitRF, resonance-final antenna function.

class QQEmitRF : public QQEmitIF {

public:

  // Names (remember to redefine both for each inherited class).
  string vinciaName() const {return "Vincia:QQEmitRF";}

  // Parton types AB -> ijk with A,i initial and B,k,j final.
  int idA() const {return 6;}
  int idB() const {return 5;}
  int id0() const {return 6;}
  int id1() const {return 21;}
  int id2() const {return 5;}

  // Mark that this function does not have a zA collinear limit.
  double zA(vector<double>) {return -1;}

  // Return this is a resonance-final antenna.
  bool isRFant() {return true;}

  // Test masses (top, gluon, bottom, and W mass).
  void getTestMasses(vector<double> &masses) {masses = {particleDataPtr->m0(6),
    0.0, particleDataPtr->m0(5), particleDataPtr->m0(24)};}

  // AP with dummy helicities.
  virtual double AltarelliParisi(vector<double> invariants,
    vector<double> masses, vector<int>, vector<int>) {
    double sjk(invariants[2]), mkres(masses[2]), z(zB(invariants)),
      mu2(mkres*mkres/sjk), Pz(dglapPtr->Pq2gq(z,9,9,9,mu2));
    return Pz/sjk;};

};

//==========================================================================

// Class QGEmitRF, resonance-final antenna function.

class QGEmitRF : public QGEmitIF {

public:

  // Names (remember to redefine both for each inherited class).
  string vinciaName() const {return "Vincia:QGEmitRF";}

  // Parton types AB -> ijk with A,i initial and B,k,j final.
  int idA() const {return 6;}
  int idB() const {return 21;}
  int ida() const {return 6;}
  int idb() const {return 21;}
  int id1() const {return 21;}

  // Mark that this function does not have a zA collinear limit.
  double zA(vector<double>) {return -1;}

  // Return this is a resonance-final antenna.
  bool isRFant() {return true;}

  // Test masses (top, gluon, gluon, X).
  void getTestMasses(vector<double> &masses) {masses = {particleDataPtr->m0(6),
    0.0, 0.0, 0.6*particleDataPtr->m0(6)};}

  // AP with dummy helicities and masses.
  virtual double AltarelliParisi(vector<double> invariants, vector<double>,
    vector<int>, vector<int>) {
    double sjk(invariants[2]), z(zB(invariants)),
      Pz(dglapPtr->Pg2gg(z, 9, 9, 9));
    return Pz/sjk;}

};

//==========================================================================

// Class XGSplitRF, resonance-final antenna function.

class XGSplitRF : public XGSplitIF {

public:

  // Names (remember to redefine both for each inherited class)
  string vinciaName() const {return "Vincia:XGSplitRF";}

  // Mark that this function does not have a zA collinear limit.
  double zA(vector<double>){ return -1;}

  // Return this is a resonance-final antenna.
  bool isRFant() {return true;}

  // Parton types AB -> ijk with A,i initial and B,k,j final.
  int idA() const {return 6;}
  int idB() const {return 21;}
  int ida() const {return 6;}
  int idb() const {return 2;}
  int id1() const {return -2;}

  // Test masses (top, gluon, gluon, X).
  void getTestMasses(vector<double> &masses) {masses = {particleDataPtr->m0(6),
    0.0, 0.0, 0.6*particleDataPtr->m0(6)};}

  // AP with dummy helicities.
  double AltarelliParisi(vector<double> invariants, vector<double> masses,
    vector<int>, vector<int>) {
    double sAK(invariants[0]), saj(invariants[1]), sjk(invariants[2]),
      mkres(masses[2]), m2q(mkres*mkres), Q2(sjk + 2*m2q), mu2(m2q/Q2),
      z((sAK+saj-Q2)/sAK), Pz(dglapPtr->Pg2qq(z, 9, 9, 9, mu2));
    return Pz/Q2;}

};

//==========================================================================

// The AntennaSetFSR class. Simple container of FF and RF antenna functions.

class AntennaSetFSR {

public:

  // Default constructor.
  AntennaSetFSR() = default;

  // Destructor, delete the antennae.
  virtual ~AntennaSetFSR() {
    for (map<int, AntennaFunction*>::iterator it = antFunPtrs.begin();
         it != antFunPtrs.end(); ++it) delete it->second;
    antFunPtrs.clear();}

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, DGLAP* dglapPtrIn);

  // Initialize antenna set (optionally with min or max variation).
  void init();

  // Function to chek if an antenna with the given index exists.
  bool exists(int iAntIn) {return antFunPtrs.count(iAntIn);}

  // Gets an antenna iAntIn from the AntennaSet.
  AntennaFunction* getAnt(int iAntIn) {
    return (exists(iAntIn)) ? antFunPtrs[iAntIn] : nullptr;}

  // Method to return all iAntPhys values that are defined in antFunPtr.
  vector<int> getIant();

  // Get Vincia name, e.g. "Vincia:QQEmitFF".
  string vinciaName(int iAntIn) {
    return exists(iAntIn) ? antFunPtrs[iAntIn]->vinciaName() : "noVinciaName";}

  // Get human name, e.g. "g/qq".
  string humanName(int iAntIn) {
    return exists(iAntIn) ? antFunPtrs[iAntIn]->humanName() : "noHumanName";}

private:

  // Use a map of AntennaFunction pointers, create them with new on
  // initialization.
  map<int,AntennaFunction*> antFunPtrs{};

  // Pointers to Pythia8 classes, needed to initialise antennae.
  bool isInitPtr{false}, isInit{false};
  Info*         infoPtr{};
  ParticleData* particleDataPtr{};
  Settings*     settingsPtr{};
  Rndm*         rndmPtr{};

  // Pointer to VINCIA DGLAP class.
  DGLAP* dglapPtr{};

  // Verbosity level
  int verbose{};

};

//==========================================================================

// The AntennaSetISR class. Simple container of II and IF antenna functions.

class AntennaSetISR {

 public:

  // Default constructor.
  AntennaSetISR() = default;

  // Destructor, delete the antennae.
  ~AntennaSetISR() {
    for (map<int, AntennaFunctionIX*>::iterator it = antFunPtrs.begin();
         it != antFunPtrs.end(); ++it) delete it->second;
    antFunPtrs.clear();}

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, DGLAP* dglapPtrIn);

  // Initialize antenna set.
  void init();

  // Function to chek if an antenna with the given index exists.
  bool exists(int iAntIn) {return antFunPtrs.count(iAntIn);}

  // Gets an antenna iAntIn from the AntennaSetISR.
  AntennaFunctionIX* getAnt(vector<AntennaFunctionIX*>::size_type iAntIn) {
    return (exists(iAntIn)) ? antFunPtrs[iAntIn] : nullptr;}

  // Method to return all iAntPhys values that are defined in antFunPtr.
  vector<int> getIant();

  // Get Vincia name, e.g. "Vincia:QQEmitII".
  string vinciaName(vector<AntennaFunctionIX*>::size_type iAntIn) {
    return exists(iAntIn) ? antFunPtrs[iAntIn]->vinciaName() : "noVinciaName";}

  // Get human name, e.g. "g/qq".
  string humanName(int iAntIn) {
    return exists(iAntIn) ? antFunPtrs[iAntIn]->humanName() : "noHumanName";}

private:

  // Use a map of AntennaFunction pointers, create them with new on
  // initialization.
  map<int,AntennaFunctionIX*> antFunPtrs{};

  // Pointers to Pythia 8 classes, needed to initialise antennae.
  bool isInitPtr{false}, isInit{false};
  Info*         infoPtr{};
  ParticleData* particleDataPtr{};
  Settings*     settingsPtr{};
  Rndm*         rndmPtr{};

  // Pointer to VINCIA DGLAP class
  DGLAP* dglapPtr{};

  // Verbosity level
  int verbose{};

};

//==========================================================================

// Class MECs, for computing matrix-element corrections to antenna
// functions.

class MECs {

public:

  // Constructor.
  MECs() {isInitPtr = false; isInit = false;}

  // Destructor.
  virtual ~MECs() {};

  // Initialize pointers.
  void initPtr(Info* infoPtrIn, VinciaMG5MEs* mg5mesPtrIn,
    VinciaCommon* vinComPtrIn);

  // Initialize pointers to antenna sets.
  void initAntPtr(AntennaSetFSR* antFSRusr, AntennaSetISR* antISRusr) {
    antSetFSR = antFSRusr; antSetISR = antISRusr;}

  // Initialize.
  void init();

  // Function to return ME class (Born type) for a parton
  // configuration. Can be called from either of the ISR::prepare() or
  // FSR::prepare() functions, or from the ISR::branch() or
  // FSR::branch() functions. Returns >= 0 if there an ME for this
  // configuration, associated with the (arbitrary) integer code label
  // denoted by the return value. If return < 0 we don't have an ME /
  // no ME should be used for this system.
  bool prepare(const int iSys, Event& event);

  // Function to assign helicities to particles (using MEs).
  bool polarise(const int iSys, Event& event);

  // Make list of particles as vector<Particle>.
  //   First 1 or 2 entries : incoming particle(s).
  //   Subseqent entries    : outgoing particles.
  // The two last arguments are optional and allow to specify a list
  // of indices to be ignored, and a set of particles to be added, e.g.
  // in the context of setting up a trial state after a branching.
  // The newly added particles are then at the end of the respective
  // lists, i.e. a newly added incoming particle is the last incoming
  // one and newly added outgoing ones are the last among the outgoing
  // ones.
  vector<Particle> makeParticleList(const int iSys, const Event& event,
    const vector<Particle> pNew = vector<Particle>(),
    const vector<int> iOld = vector<int>());

  // Check if state already has helicities. Return true if any
  // particle in state already has a specified helicity.
  bool isPolarised(int iSys, Event& event, bool checkIncoming);

  // Wrapper function to return a specific antenna function.
  AntennaFunction* getAntFSR(const int iAnt) {
    return antSetFSR->getAnt(iAnt); }
  AntennaFunctionIX* getAntISR(const int iAnt) {
    return antSetISR->getAnt(iAnt); }

  // Function to determine if MECs are requested at this order for this system.
  bool doMEC(const int iSys, const int nBranch);

  // Get squared matrix element.
  double getME2(const vector<Particle> parts, const int nIn);
  double getME2(const int iSys, const Event& event);

  // Get matrix element correction factor (TODO: dummy implementation).
  double getMEC(const int, const Event&) {return 1.;}

  // Return number of partons added since Born (as defined by prepare).
  int sizeOutBorn(const int iSys) {return sizeOutBornSav[iSys];}

  // Function to set level of verbose output.
  void setVerbose(const int verboseIn) {verbose = verboseIn;}

  // Header.
  void header();

private:

  // Verbosity level.
  int verbose;

  // Is initialized.
  bool isInitPtr, isInit;

  // Pointers to PYTHIA objects.
  Info*          infoPtr;
  Settings*      settingsPtr;
  Rndm*          rndmPtr;
  PartonSystems* partonSystemsPtr;

  // Pointers to VINCIA objects.
  VinciaMG5MEs*  mg5mesPtr;
  VinciaCommon*  vinComPtr;

  // Antenna sets.
  AntennaSetFSR* antSetFSR;
  AntennaSetISR* antSetISR;

  // Matching settings.
  bool matchingFullColour;
  int  maxMECs2to1, maxMECs2to2, maxMECs2toN, maxMECsResDec, maxMECsMPI;
  int  nFlavZeroMass;

  // Map from iSys to Born multiplicity and ME class; set in prepare().
  map<int,int> sizeOutBornSav;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_VinciaAntennaFunctions_H
