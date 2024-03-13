// HINucleusModel.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the HIUserHooks class and a
// set of other classes used inside Pythia8 to model collisions
// involving heavy ions.
// Nucleon: represents a proton or a neutron inside a necleus.
// NucleusModel: models the Nucleon distribution in a nucleus.
// WoodsSaxonModel: NucleusModel implementing a simple Woods-Saxon.
// GLISSANDOModel: NucleusModel implementing the GLISSANDO prescription.

#ifndef Pythia8_HINucleusModel_H
#define Pythia8_HINucleusModel_H

#include "Pythia8/HIBasics.h"

namespace Pythia8 {

//==========================================================================

// The Nucleon class represents a nucleon in a nucleus. It has an id
// number (proton or neutron) an impact parameter position (absolute
// and relative to the nucleus center), a status and a state to be
// defined and used by a SubCollisionModel.

class Nucleon {

public:

  // Enum for specifying the status of a nucleon.
  enum Status : int {
    UNWOUNDED = 0,  // The nucleon is not wounded.
    ELASTIC = 1,    // The nucleon is elastically scattered.
    DIFF = 2,       // The nucleon is diffractively wounded.
    ABS = 3         // The nucleon is absorptively wounded.
  };

  // The state of a nucleon is a general vector of doubles.
  typedef vector<double> State;

  // The constuctor takes a particle id and a position in impact
  // parameter relative to the nucleus center as arguments.
  Nucleon(int idIn = 0, int indexIn = 0, const Vec4 & pos = Vec4())
    : idSave(idIn), indexSave(indexIn), nPosSave(pos), bPosSave(pos),
      statusSave(UNWOUNDED), eventp(0), isDone(0) {}

  // Accessor functions:

  // The particle id of the nucleon.
  int id() const { return idSave; }

  // The index of the nucleon in the nucleus.
  int index() const { return indexSave; }

  // The position of this nucleon relative to the nucleus center.
  const Vec4 & nPos() const { return nPosSave; }

  // The absolute position in impact parameter space.
  const Vec4 & bPos() const { return bPosSave; }

  // Shift the absolute position in impact parameter space.
  void bShift(const Vec4 & bvec) { bPosSave += bvec; }

  // The status of the nucleon.
  Nucleon::Status status() const { return statusSave; }

  // Check if nucleon has been assigned.
  bool done() const { return isDone; }

  // The event this nucleon is assigned to.
  EventInfo * event() const { return eventp; }

  // The physical state of the incoming nucleon.
  const State & state() const { return stateSave; }

  // Return an alternative state.
  const State & altState(int i = 0) {
    static State nullstate;
    return i < int(altStatesSave.size())? altStatesSave[i]: nullstate;
  }

  // Manipulating functions:

  // Set the status.
  void status(Nucleon::Status s) { statusSave = s; }

  // Set the physical state.
  void state(State s) { stateSave = s; }

  // Add an alternative state.
  void addAltState(State s) { altStatesSave.push_back(s); }

  // Select an event for this nucleon.
  void select(EventInfo & evp, Nucleon::Status s) {
    eventp = &evp;
    isDone = true;
    status(s);
  }

  // Select this nucleon to be assigned to an event.
  void select() { isDone = true; }

  // Print out debugging information.
  void debug();

  // Reset the states and status.
  void reset() {
    statusSave = UNWOUNDED;
    altStatesSave.clear();
    bPosSave = nPosSave;
    isDone = false;
    eventp = 0;
  }

private:

  // The type of nucleon.
  int idSave;

  // The index of this nucleon.
  int indexSave;

  // The position in impact parameter relative to the nucleus center.
  Vec4 nPosSave;

  // The absolute position in impact parameter.
  Vec4 bPosSave;

  // The status.
  Nucleon::Status statusSave;

  // The state of this nucleon.
  State stateSave;

  // Alternative states to be used to understand fluctuations in the
  // state of this nucleon.
  vector<State> altStatesSave;

  // Pointer to the event this nucleon ends up in.
  EventInfo * eventp;

  // True if this nucleon has been assigned to an event.
  bool isDone;

};

//==========================================================================


class Nucleus {

public:

  // Default constructor.
  Nucleus() = default;

  // Constructor with nucleons and impact parameter.
  Nucleus(vector<Nucleon> nucleons, Vec4 bPos) : bPosSave(bPos) {
    nucleonsSave = make_shared<vector<Nucleon>>(nucleons);
    for (Nucleon& nucleon : *nucleonsSave) {
      nucleon.reset();
      nucleon.bShift(bPos);
    }
  }

  // Iterate over nucleons.
  vector<Nucleon>::iterator begin() { return nucleonsSave->begin(); }
  vector<Nucleon>::iterator end() { return nucleonsSave->end(); }
  vector<Nucleon>::const_iterator begin() const {return nucleonsSave->begin();}
  vector<Nucleon>::const_iterator end() const {return nucleonsSave->end();}

private:

  // Saved nucleons and impact parameter.
  shared_ptr<vector<Nucleon>> nucleonsSave;
  Vec4 bPosSave;

};

//==========================================================================

// This class generates the impact parameter distribution of nucleons
// in a nucleus.

class NucleusModel {

public:

  // Default constructor giving the nucleus id and an optional
  // radius (in femtometer).
  NucleusModel() : isProj(true), idSave(0), ISave(0), ASave(0),
     ZSave(0), LSave(0), RSave(0.0), settingsPtr(0),
     rndmPtr(0) {}

  // Virtual destructor.
  virtual ~NucleusModel() {}

  static shared_ptr<NucleusModel> create(int model);

  // Init method.
  void initPtr(int idIn, bool isProjIn, Info& infoIn);
  virtual bool init() { return true; }

  // Set (new) nucleon momentum.
  virtual void setPN(const Vec4 & pNIn) { pNSave = pNIn; }

  // Produce an instance of the incoming nucleon.
  virtual Particle produceIon();

  // Generate a vector of nucleons according to the implemented model
  // for a nucleus given by the PDG number.
  virtual vector<Nucleon> generate() const = 0;

  // Accessor functions.
  int id() const { return idSave; }
  int I() const { return ISave; }
  int A() const { return ASave; }
  int Z() const { return ZSave; }
  int L() const { return LSave; }
  double R() const { return RSave; }

protected:

  // Projectile or target.
  bool isProj;

  // The nucleus.
  int idSave;

  // Cache information about the nucleus.
  int ISave, ASave, ZSave, LSave;

  // The estimate of the nucleus radius.
  double RSave;

  // The mass of the nucleus and its nucleons.
  double mSave{};

  // The nucleon beam momentum.
  Vec4 pNSave{};

  // Pointers to useful objects.
  Info* infoPtr;
  Settings* settingsPtr;
  Rndm* rndmPtr;
  Logger* loggerPtr;

};

//==========================================================================

// A nucleus model defined by an external file to be read in, containing
// x,y,z coordinates of the nucleons.

class ExternalNucleusModel : public NucleusModel {

public:

  // Default constructor.
  ExternalNucleusModel() : fName(""), doShuffle(true), nUsed(0) {}

  // Initialize class. Read in file to buffer.
  bool init() override;

  // Generate a vector of nucleons according to the implemented model
  // for a nucleus given by the PDG number.
  vector<Nucleon> generate() const override;

private:

  // The filename to read from.
  string fName;

  // Shuffle configurations.
  bool doShuffle;

  // The read nucleon configurations. Time component is always zero.
  mutable vector<vector<Vec4> > nucleonPositions;

  // The number of configurations used so far.
  mutable size_t nUsed;

};

//==========================================================================

// A NucleusModel which allows for a hard core, optionally a Gaussian
// hard core. This is an abstract class intended as a base class for
// models with this functionality.

class HardCoreModel : public NucleusModel {

public:

  // Default constructor.
  HardCoreModel() : useHardCore(), gaussHardCore(), hardCoreRadius(0.9) {}

  // Virtual destructor.
  virtual ~HardCoreModel() {}

  // Initialize the parameters for hard core generation.
  // To be called in init() in derived classes.
  void initHardCore();

  // Get the radius of the hard core. If using a Gaussian hard core, the
  // radius is distributed according to a 1D Gaussian.
  double rSample() const {
    if (gaussHardCore) return hardCoreRadius * abs(rndmPtr->gauss());
    return hardCoreRadius;}

protected:

  // Use the hard core or not.
  bool useHardCore;

  // Use a Gaussian hard core.
  bool gaussHardCore;

  // The radius or width of the hard core.
  double hardCoreRadius;

};

//==========================================================================

// A general Woods-Saxon distributed nucleus.

class WoodsSaxonModel : public HardCoreModel {

public:

  // Virtual destructor.
  virtual ~WoodsSaxonModel() {}

  // The default constructor needs a nucleus id, a radius, R, and a
  // "skin width", a (both length in femtometers).
  WoodsSaxonModel(): aSave(0.0), intlo(0.0),
                    inthi0(0.0), inthi1(0.0), inthi2(0.0) {}

  // Initialize parameters.
  bool init() override;

  // Generate all the nucleons.
  vector<Nucleon> generate() const override;

  // Accessor functions.
  double a() const { return aSave; }

protected:

  // Generate the position of a single nucleon. (The time component
  // is always zero).
  Vec4 generateNucleon() const;

  // Calculate overestimates for sampling.
  void overestimates() {
    intlo = R()*R()*R()/3.0;
    inthi0 = a()*R()*R();
    inthi1 = 2.0*a()*a()*R();
    inthi2 = 2.0*a()*a()*a();
  }

protected:

  // The nucleus radius, skin depth parameter, and hard core nucleon radius.
  double aSave;

private:

  // Cashed integrals over the different parts of the over estimating
  // functions.
  double intlo, inthi0, inthi1, inthi2;

};


//==========================================================================

// The GLISSANDOModel is a specific parameterization of a Woods-Saxon
// potential for A>16. It is described in arXiv:1310.5475 [nucl-th].

class GLISSANDOModel : public WoodsSaxonModel {

public:

  // Default constructor.
  GLISSANDOModel() {}

  // Virtual destructor.
  virtual ~GLISSANDOModel() {}

  // Initialize.
  bool init() override;

};

//==========================================================================

// A Harmonic-Oscillator Shell model for light nuclei.

class HOShellModel : public HardCoreModel {

public:

  // Default constructor.
  HOShellModel(): nucleusChR(), protonChR(), C2() {}

  // Destructor.
  virtual ~HOShellModel() {}

  // Initialize, set up parameters.
  virtual bool init() override;

  // Generate a vector of nucleons according to the implemented model
  // for a nucleus given by the PDG number.
  virtual vector<Nucleon> generate() const override;

protected:

  // Generate the position of a single nucleon. (The time component
  // is always zero).
  virtual Vec4 generateNucleon() const;

  // The density function.
  double rho(double r) const {
    double pref = 4./(pow(sqrt(M_PI * C2),3)) * (1 + (A() - 4.)/6. * r*r/C2);
    return pref * exp(-r*r / C2);
  };

  // Nucleus charge radius.
  double nucleusChR;

  // Nucleon charge radius.
  double protonChR;

  // C2 parameter.
  double C2;

  // Maximum rho for these parameters.
  double rhoMax;

};

//==========================================================================

// The Hulthen potential for deuterons.

class HulthenModel : public NucleusModel {

public:

  // Default constructor.
  HulthenModel(): hA(), hB() {}

  // Virtual destructor.
  virtual ~HulthenModel() {}

  virtual bool init() override;

  // Generate a vector of nucleons according to the Hulthen potential.
  virtual vector<Nucleon> generate() const override;

protected:

  // The (normalized) density function.
  double rho(double r) const {
    double pref = (2*hA*hB*(hA + hB))/pow2(hA - hB);
    double exps = exp(-2.*hA*r) + exp(-2.*hB*r) - 2.*exp(-(hA+hB)*r);
    return pref * exps;
  };

  // Parameters of the Hulthen model.
  double hA;
  double hB;

};

//==========================================================================

// A Gaussian distribution for light nuclei.

class GaussianModel : public HardCoreModel {

public:

  // Default constructor.
  GaussianModel(): nucleusChR() {}

  // Destructor.
  virtual ~GaussianModel() {}

  virtual bool init() override;

  // Generate a vector of nucleons according to the implemented model
  // for a nucleus given by the PDG number.
  virtual vector<Nucleon> generate() const override;

protected:

  // Generate the position of a single nucleon. (The time component
  // is always zero).
  virtual Vec4 generateNucleon() const;

  // Nucleus charge radius.
  double nucleusChR;

};

//==========================================================================

// A model for nuclei clustered in smaller nuclei.

class ClusterModel : public HardCoreModel {

public:

  // Contructor.
  ClusterModel() {}

  // Virtual destructor.
  virtual ~ClusterModel() {}

  // Initialize parameters.
  virtual bool init() override;

  // Generate a vector of nucleons. Note that this model
  // is only implemented for XX, YY ZZ.
  virtual vector<Nucleon> generate() const override;

private:

  // The model to generate clusters from.
  unique_ptr<NucleusModel> nModelPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HINucleusModel_H
