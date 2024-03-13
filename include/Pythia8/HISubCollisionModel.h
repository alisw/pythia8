// HISubCollisionModel.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the ImpactParmeterGenerator,
// SubCollision, and SubCollisionModel classes, as well as a set of
// subclasses of SubCollisionModel.
//
// ImpactParameterGenerator: distributes nuclei in impact parameter space.
// SubCollision: a collision between a projectile and a target Nucleon.
// SubCollisionModel: Models the collision probabilities of nucleons.
// BlackSubCollisionModel: A very simple SubCollisionModel.
// NaiveSubCollisionModel: A very simple SubCollisionModel.
// DoubleStrikmanSubCollisionModel: A more advanced SubCollisionModel.

#ifndef Pythia8_HISubCollisionModel_H
#define Pythia8_HISubCollisionModel_H

#include "Pythia8/Pythia.h"
#include "Pythia8/HINucleusModel.h"

namespace Pythia8 {

//==========================================================================

// SubCollision represents a possible collision between a projectile
// and a target nucleon.

class SubCollision {

public:

  // This defines the type of a binary nucleon collison.
  enum CollisionType {
    NONE,       // This is not a collision.
    ELASTIC,    // This is an elastic scattering
    SDEP,       // The projectile is diffractively excited.
    SDET,       // The target is diffractively excited.
    DDE,        // Both projectile and target are diffractively excited.
    CDE,        // Both excited but with central diffraction.
    ABS         // This is an absorptive (non-diffractive) collision.
  };

  // Constructor with configuration.
  SubCollision(Nucleon & projIn, Nucleon & targIn,
               double bIn, double bpIn, CollisionType typeIn)
    : proj(&projIn), targ(&targIn), b(bIn), bp(bpIn), type(typeIn) {}

  // Default constructor.
  SubCollision()
    : proj(0), targ(0), b(0.0), bp(0.0), type(NONE) {}

  // Used to order sub-collisions in a set.
  bool operator< (const SubCollision& s) const { return b < s.b; }

  // Return 0 if neither proj or target are neutrons, 1 if target is
  // neutron, 2 if projectile is neutron, and 3 if both are neutrons.
  int nucleons() const {return ( abs(targ->id()) == 2112? 1: 0 ) +
      ( abs(proj->id()) == 2112? 2: 0 );}

  // The projectile nucleon.
  Nucleon* proj;

  // The target nucleon.
  Nucleon* targ;

  // The impact parameter distance between the nucleons in femtometer.
  double b;

  // The impact parameter distance between the nucleons scaled like
  // in Pythia to have unit average for non-diffractive collisions.
  double bp;

  // The type of collision.
  CollisionType type;

};

//==========================================================================

// The SubCollisionSet gives a set of subcollisions between two nuclei.

class SubCollisionSet {

public:

  // Default constructor.
  SubCollisionSet() = default;

  // Constructor with subcollisions.
  SubCollisionSet(multiset<SubCollision> subCollisionsIn, double TIn)
    : subCollisionsSave(subCollisionsIn), TSave(TIn) {}

  // Reset the subcollisions.
  bool empty() const { return subCollisionsSave.empty(); }

  // The summed elastic amplitude.
  double T() const { return TSave; }

  // Iterators over the subcollisions.
  multiset<SubCollision>::const_iterator begin() const {
    return subCollisionsSave.begin(); }
  multiset<SubCollision>::const_iterator end() const {
    return subCollisionsSave.end(); }

private:

  // Saved subcollisions.
  multiset<SubCollision> subCollisionsSave;
  double TSave;

};

//==========================================================================

// The SubCollisionModel is is able to model the collision between two
// nucleons to tell which type of collision has occurred. The model
// may manipulate the corresponding state of the nucleons.

class SubCollisionModel {

public:

  // Internal class to report cross section estimates.
  struct SigEst {
    // The cross sections (tot, nd, dd, sdp, sdt, cd, el, bslope).
    vector<double> sig;

    // The estimated error (squared)
    vector<double> dsig2;

    // Which cross sections were actually fitted
    vector<bool> fsig;

    // The estimate of the average (and squared error) impact
    // parameter for inelastic non-diffractive collisions.
    double avNDb, davNDb2;

    // Constructor for zeros.
    SigEst(): sig(8, 0.0), dsig2(8, 0.0), fsig(8, false),
              avNDb(0.0), davNDb2(0.0) {}

  };

  // The default constructor is empty.
  SubCollisionModel(int nParm): sigTarg(8, 0.0), sigErr(8, 0.05),
    parmSave(nParm),
    NInt(100000), NPop(20), sigFuzz(0.2), impactFudge(1),
    fitPrint(true), avNDb(1.0*femtometer),
    projPtr(), targPtr(), sigTotPtr(), settingsPtr(), infoPtr(), rndmPtr() {}

  // Virtual destructor.
  virtual ~SubCollisionModel() {}

  // Create a new SubCollisionModel of the given model.
  static shared_ptr<SubCollisionModel> create(int model);

  // Virtual init method.
  virtual bool init(int idAIn, int idBIn, double eCMIn);

  // Initialize the pointers.
  void initPtr(NucleusModel & projIn, NucleusModel & targIn,
               SigmaTotal & sigTotIn, Settings & settingsIn,
               Info & infoIn, Rndm & rndmIn) {
    projPtr = &projIn;
    targPtr = &targIn;
    sigTotPtr = &sigTotIn;
    settingsPtr = &settingsIn;
    infoPtr = &infoIn;
    rndmPtr = &rndmIn;
    loggerPtr = infoIn.loggerPtr;
  }

  // Access the nucleon-nucleon cross sections assumed
  // for this model.

  // The target total nucleon-nucleon cross section.
  double sigTot() const { return sigTarg[0]; }

  // The target elastic cross section.
  double sigEl() const { return sigTarg[6]; }

  // The target central diffractive excitation cross section.
  double sigCDE() const { return sigTarg[5]; }

  // The target single diffractive excitation cross section (both sides).
  double sigSDE() const { return sigTarg[3] + sigTarg[4]; }

  // The target single diffractive excitation cross section (projectile).
  double sigSDEP() const { return sigTarg[3]; }

  // The target single diffractive excitation cross section (target).
  double sigSDET() const { return sigTarg[4]; }

  // The target double diffractive excitation cross section.
  double sigDDE() const { return sigTarg[2]; }

  // The target non-diffractive (absorptive) cross section.
  double sigND() const { return sigTarg[1]; }

  // The target elastic b-slope parameter.
  double bSlope() const { return sigTarg[7]; }

  // Return the average non-diffractive impact parameter.
  double avNDB() const { return avNDb; }

  // Update internally stored cross sections.
  void updateSig();

  // Calculate the Chi2 for the given cross section estimates.
  double Chi2(const SigEst & sigs, int npar) const;

  // Set beam kinematics.
  void setKinematics(double eCMIn);

  // Use a genetic algorithm to fit the parameters.
  bool evolve(int nGenerations, double eCM);

  // Get the number of free parameters for the model.
  int nParms() const { return parmSave.size(); }

  // Set the parameters of this model.
  void setParm(const vector<double>& parmIn) {
    for (size_t i = 0; i < parmSave.size(); ++i)
      parmSave[i] = parmIn[i];
  }

  // Get the current parameters of this model.
  vector<double> getParm() const { return parmSave; }

  // Get the minimum allowed parameter values for this model.
  virtual vector<double> minParm() const = 0;

  // Get the default parameter values for this model.
  virtual vector<double> defParm() const = 0;

  // Get the maximum allowed parameter values for this model.
  virtual vector<double> maxParm() const = 0;

  // Take two nuclei and produce the corresponding subcollisions. The states
  // of the nucleons may be changed if fluctuations are allowed by the model.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) = 0;

  // Calculate the cross sections for the given set of parameters.
  virtual SigEst getSig() const = 0;

private:

  // The nucleon-nucleon cross sections targets for this model
  // (tot, nd, dd, sdp, sdt, cd, el, bslope) and the required precision.
  vector<double> sigTarg, sigErr;

  // Generate parameters based on run settings and the evolutionary algorithm.
  bool genParms();

  // Save/load parameter configuration from disk.
  bool saveParms(string fileName) const;
  bool loadParms(string fileName);

protected:

  // Saved parameters.
  vector<double> parmSave;

  // The parameters stearing the fitting of internal parameters to
  // the different nucleon-nucleon cross sections.
  int NInt, NPop;
  double sigFuzz;
  double impactFudge;
  bool fitPrint;

  // The estimated average impact parameter distance (in femtometer)
  // for absorptive collisions.
  double avNDb;

  // Info from the controlling HeavyIons object.
  NucleusModel* projPtr;
  NucleusModel* targPtr;
  SigmaTotal* sigTotPtr;
  Settings* settingsPtr;
  Info* infoPtr;
  Rndm* rndmPtr;
  Logger* loggerPtr;

  // For variable energies.
  int idASave, idBSave;
  bool doVarECM;
  double eMin{}, eMax{};
  int eCMPts;
  vector<LogInterpolator> subCollParms;

};

//==========================================================================

// The most naive sub-collision model, assuming static nucleons and
// an absorptive cross section equal to the total inelastic. No
// fluctuations, meaning no diffraction.

class BlackSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  BlackSubCollisionModel() : SubCollisionModel(0) {}

  // Virtual destructor.
  virtual ~BlackSubCollisionModel() override {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return vector<double>(); }
  vector<double> defParm() const override { return vector<double>(); }
  vector<double> maxParm() const override { return vector<double>(); }

  // Get cross sections used by this model.
  virtual SigEst getSig() const override {
    SigEst s;
    s.sig[0] = sigTot();
    s.sig[1] = sigND();
    s.sig[6] = s.sig[0] - s.sig[1];
    s.sig[7] = bSlope();
    return s;
  }

  // Take two nuclei and return the corresponding sub-collisions.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

};

//==========================================================================

// A very simple sub-collision model, assuming static nucleons and
// just assuring that the individual nucleon-nucleon cross sections
// are preserved.

class NaiveSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  NaiveSubCollisionModel() : SubCollisionModel(0) {}

  // Virtual destructor.
  virtual ~NaiveSubCollisionModel() override {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return vector<double>(); }
  vector<double> defParm() const override { return vector<double>(); }
  vector<double> maxParm() const override { return vector<double>(); }

  // Get cross sections used by this model.
  virtual SigEst getSig() const override {
    SigEst s;
    s.sig[0] = sigTot();
    s.sig[1] = sigND();
    s.sig[3] = sigSDEP();
    s.sig[4] = sigSDET();
    s.sig[2] = sigDDE();
    s.sig[6] = sigEl();
    s.sig[7] = bSlope();
    return s;
  }

  // Take two nuclei and return the corresponding sub-collisions.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

};

//==========================================================================

// A base class for sub-collision models where each nucleon has a
// fluctuating "radius". The base model has two parameters, sigd and alpha,
// which are used for opacity calculations. Subclasses may have additional
// parameters to describe the radius distributions of that specific model.

class FluctuatingSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  FluctuatingSubCollisionModel(int nParmIn, int modein)
    : SubCollisionModel(nParmIn + 2),
      sigd(parmSave[nParmIn]), alpha(parmSave[nParmIn + 1]),
      opacityMode(modein) {}

  // Virtual destructor.
  virtual ~FluctuatingSubCollisionModel() override {}

  // Take two nuclei and pick specific states for each nucleon,
  // then get the corresponding sub-collisions.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

  // Calculate the cross sections for the given set of parameters.
  virtual SigEst getSig() const override;

protected:

  // Pick a radius for the nucleon, depending on the specific model.
  virtual double pickRadiusProj() const = 0;
  virtual double pickRadiusTarg() const = 0;

private:

  // Saturation scale of the nucleus.
  double& sigd;

  // Power of the saturation scale
  double& alpha;

  // Optional mode for opacity.
  int opacityMode;

  // The opacity of the collision at a given sigma.
  double opacity(double sig) const {
    sig /= sigd;
    if ( opacityMode == 1 ) sig = 1.0/sig;
    return sig > numeric_limits<double>::epsilon() ?
      pow(-expm1(-1.0/sig), alpha) : 1.0;
  }

  // Return the elastic amplitude for a projectile and target state
  // and the impact parameter between the corresponding nucleons.
  double Tpt(const Nucleon::State & p,
             const Nucleon::State & t, double b) const {
    double sig = M_PI*pow2(p[0] + t[0]);
    double grey = opacity(sig);
    return sig/grey > b*b*2.0*M_PI? grey: 0.0;
  }

};

//==========================================================================

// A sub-collision model where each nucleon has a fluctuating
// "radius" according to a Strikman-inspired distribution.

class DoubleStrikmanSubCollisionModel : public FluctuatingSubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  DoubleStrikmanSubCollisionModel(int modeIn = 0)
    : FluctuatingSubCollisionModel(1, modeIn), k0(parmSave[0]) {}

  // Virtual destructor.
  virtual ~DoubleStrikmanSubCollisionModel() override {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return {  0.01,  1.0,  0.0  }; }
  vector<double> defParm() const override { return {  2.15, 17.24, 0.33 }; }
  vector<double> maxParm() const override { return { 60.00, 60.0, 20.0  }; }

protected:

  double pickRadiusProj() const override { return rndmPtr->gamma(k0, r0()); }
  double pickRadiusTarg() const override { return rndmPtr->gamma(k0, r0()); }

private:

  // The power in the Gamma distribution.
  double& k0;

  // Return the average radius deduced from other parameters and
  // the total cross section.
  double r0() const {
    return sqrt(sigTot() / (M_PI * (2.0 * k0 + 4.0 * k0 * k0)));
  }

};

//==========================================================================

// ImpactParameterGenerator is able to generate a specific impact
// parameter together with a weight such that aweighted average over
// any quantity X(b) corresponds to the infinite integral over d^2b
// X(b). This base class gives a Gaussian profile, d^2b exp(-b^2/2w^2).

class ImpactParameterGenerator {

public:

  // The default constructor takes a general width (in femtometers) as
  // argument.
  ImpactParameterGenerator()
    : widthSave(0.0), collPtr(0), projPtr(0), targPtr(0),
      settingsPtr(0), rndmPtr(0) {}

  // Virtual destructor.
  virtual ~ImpactParameterGenerator() {}

  // Virtual init method.
  virtual bool init();
  void initPtr(Info & infoIn, SubCollisionModel & collIn,
    NucleusModel & projIn, NucleusModel & targIn);

  // Return a new impact parameter and set the corresponding weight provided.
  virtual Vec4 generate(double & weight) const;

  // Set the width (in femtometers).
  void width(double widthIn) { widthSave = widthIn; }

  // Get the width.
  double width() const { return widthSave; }

  // Update width based on the associated subcollision and nucleus models.
  void updateWidth();

private:

  // The width of a distribution.
  double widthSave;

protected:

  // Pointers from the controlling HeavyIons object.
  Info* infoPtr;
  SubCollisionModel* collPtr;
  NucleusModel* projPtr;
  NucleusModel* targPtr;
  Settings* settingsPtr;
  Rndm* rndmPtr;
  Logger* loggerPtr;

};

//==========================================================================

// A sub-collision model where each nucleon fluctuates independently
// according to a log-normal distribution. Nucleons in the projectile and
// target may fluctuate according to different parameters, which is relevant
// e.g. for hadron-ion collisions with generic hadron species.

class LogNormalSubCollisionModel : public FluctuatingSubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  LogNormalSubCollisionModel(int modeIn = 0)
    : FluctuatingSubCollisionModel(4, modeIn),
    kProj(parmSave[0]), kTarg(parmSave[1]),
    rProj(parmSave[2]), rTarg(parmSave[3]) {}

  // Virtual destructor.
  virtual ~LogNormalSubCollisionModel() {}

  //virtual SigEst getSig() const override;

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override {
    return { 0.01, 0.01, 0.10, 0.10,  1.00, 0.00 }; }
  vector<double> defParm() const override {
    return { 1.00, 1.00, 0.54, 0.54, 17.24, 0.33 }; }
  vector<double> maxParm() const override {
    return { 2.00, 2.00, 4.00, 4.00, 20.00, 2.00 }; }

protected:

  double pickRadiusProj() const override { return pickRadius(kProj, rProj); }
  double pickRadiusTarg() const override { return pickRadius(kTarg, rTarg); }

private:

  // The standard deviation of each log-normal distribution.
  double& kProj;
  double& kTarg;

  // The mean radius of each nucleon.
  double& rProj;
  double& rTarg;

  double pickRadius(double k0, double r0) const {
    double logSig = log(M_PI * pow2(r0)) + k0 * rndmPtr->gauss();
    return sqrt(exp(logSig) / M_PI);
  }
};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HISubCollisionModel_H
