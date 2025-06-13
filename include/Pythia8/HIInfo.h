// HIInfo.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the HIInfo, EventInfo and
// HIUserHooks classes.
//
// EventInfo: stores full nucleon-nucleon events with corresponding Info.
// HIInfo: info about a Heavy Ion run and its produced events.
// HIUserHooks: User hooks for HeavyIons models.

#ifndef Pythia8_HIInfo_H
#define Pythia8_HIInfo_H

#include "Pythia8/Pythia.h"
#include "Pythia8/HIBasics.h"
#include "Pythia8/HISubCollisionModel.h"

namespace Pythia8 {

//==========================================================================

// Class for collecting info about a Heavy Ion run and its produced
// events.

class HIInfo {

public:

  friend class HeavyIons;
  friend class Angantyr;

  // Constructor.
  HIInfo() = default;

  // The impact-parameter distance in the current event.
  double b() const { return bSave; }

  // The impact parameter angle.
  double phi() const { return phiSave; }

  // The summed elastic amplitude.
  double T() { return TSave; }

  // The average NN non-diffractive impact parameter to be used to
  // communicate to Pythia's MPI machinery.
  double avNDb() const { return avNDbSave; }

  // The number of attempted impact parameter points.
  long nAttempts() const { return NSave; }

  // The number of produced events.
  long nAccepted() const { return NAccSave; }

  // The total number of sub-collisions in the event.
  int nCollTot() const { return nCollSave[0]; }

  // The number of non-diffractive sub-collisions in the event
  int nCollND() const { return nCollSave[1]; }

  // The number of single diffractive projectile excitation
  // sub-collisions in the event
  int nCollSDP() const { return nCollSave[2]; }

  // The number of single diffractive target excitation
  // sub-collisions in the event
  int nCollSDT() const { return nCollSave[3]; }

  // The number of double diffractive sub-collisions in the event
  int nCollDD() const { return nCollSave[4]; }

  // The number of central diffractive sub-collisions in the event
  int nCollCD() const { return nCollSave[5]; }

  // The number of elastic sub-collisions in the event
  int nCollEl() const { return nCollSave[6]; }

  // The total number of interacting projectile nucleons in the event.
  int nPartProj() const { return nProjSave[0]; }

  // The number of absorptively wounded projectile nucleons in the event.
  int nAbsProj() const { return nProjSave[1]; }

  // The number of diffractively wounded projectile nucleons in the event.
  int nDiffProj() const { return nProjSave[2]; }

  // The number of elastically scattered target nucleons in the event.
  int nElProj() const { return nProjSave[3]; }

  // The total number of interacting target nucleons in the event.
  int nPartTarg() const { return nTargSave[0]; }

  // The number of absorptively wounded target nucleons in the event.
  int nAbsTarg() const { return nTargSave[1]; }

  // The number of diffractively wounded target nucleons in the event.
  int nDiffTarg() const { return nTargSave[2]; }

  // The number of elastically scattered target nucleons in the event.
  int nElTarg() const { return nTargSave[3]; }

  // The weight for this collision.
  double weight() const { return weightSave; }

  // The sum of weights of the produced events.
  double weightSum() const { return weightSumSave; }

  // The number of failed nucleon excitations in the current event.
  int nFail() const {
    return nFailSave;
  }

  // Register a failed nucleon excitation.
  void failedExcitation(const SubCollision& subColl) {
    subColl.failed = true;
    ++nFailSave;
  }

  // Reset Glauber statistics.
  void glauberReset();

  // The total cross section from the Glauber calculation.
  double glauberTot() const {
    return sigmaTotSave*FMSQ2MB;
  }

  // The error on the total cross section from the Glauber
  // calculation.
  double glauberTotErr() const {
    return sqrt(sigErr2TotSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The non-diffractive cross section from the Glauber calculation.
  double glauberND() const {
    return sigmaNDSave*FMSQ2MB;
  }

  // The error on the non-diffractive cross section from the Glauber
  // calculation.
  double glauberNDErr() const {
    return sqrt(sigErr2NDSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The total inelastic cross section from the Glauber calculation.
  double glauberINEL() const {
    return sigmaINELSave*FMSQ2MB;
  }

  // The error on the total inelastic cross section from the Glauber
  // calculation.
  double glauberINELErr() const {
    return sqrt(sigErr2INELSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The elastic cross section from the Glauber calculation.
  double glauberEL() const {
    return sigmaELSave*FMSQ2MB;
  }

  // The error on the elastic cross section from the Glauber
  // calculation.
  double glauberELErr() const {
    return sqrt(sigErr2ELSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The diffractive taget excitation cross section from the Glauber
  // calculation.
  double glauberDiffT() const {
    return sigmaDiffTSave*FMSQ2MB;
  }

  // The error on the diffractive taget excitation cross section from the
  // Glauber calculation.
  double glauberDiffTErr() const {
    return sqrt(sigErr2DiffTSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The diffractive projectile excitation cross section from the Glauber
  // calculation.
  double glauberDiffP() const {
    return sigmaDiffPSave*FMSQ2MB;
  }

  // The error on the diffractive projectile excitation cross section
  // from the Glauber calculation.
  double glauberDiffPErr() const {
    return sqrt(sigErr2DiffPSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The doubly diffractive excitation cross section from the Glauber
  // calculation.
  double glauberDDiff() const {
    return sigmaDDiffSave*FMSQ2MB;
  }

  // The error on the doubly diffractive excitation cross section from
  // the Glauber calculation.
  double glauberDDiffErr() const {
    return sqrt(sigErr2DDiffSave/max(1.0, double(NSave)))*FMSQ2MB;
  }

  // The elastic slope parameter.
  double glauberBSlope() const {
    return slopeSave/(sigmaTotSave*pow2(HBARC));
  }

  // The error on the elastic slope parameter.
  double glauberBSlopeErr() const {
    return sqrt((slopeErrSave/pow2(slopeSave) +
                 sigErr2TotSave/pow2(sigmaTotSave))/max(1.0, double(NSave)))*
                glauberBSlope();
  }

private:

  // Register a tried impact parameter point giving the total elastic
  // amplitude, the impact parameter and impact parameter generation
  // weight, together with the cross section scale.
  void addAttempt(double T, double bin, double phiin,
                  double bweight, double xSecScale);

  // Reweight event for whatever reason.
  void reweight(double w) {
    weightSave *= w;
  }

  // Select the primary process.
  void select(Info & in) {
    primInfo = in;
    primInfo.hiInfo = this;
  }

  // Accept an event and update statistics in info.
  void accept();

  // Reject an attmpted event.
  void reject() {}

  // Register Glauber statistics of the event.
  void glauberStatistics();

  // Collect running averages for Glauber cross section:
  static void runAvg(double & sig, double & sigErr2, double N,
    double w) {
    double delta = w - sig;
    sig += delta/N;
    sigErr2 += (delta*(w - sig) - sigErr2)/N;
  }

  // Id of the colliding nuclei.
  int idProjSave = 0, idTargSave = 0;

  // Impact parameter.
  double bSave = 0.0;
  double phiSave = 0.0;

  // Amplitude, attempts and weight.
  long NSave = 0, NAccSave = 0;
  double TSave = 0.0;
  // OBS: Cross sections should be accessed through the main info class.
  double sigmaTotSave = {}, sigmaNDSave = {},
         sigmaELSave = {}, sigmaINELSave = {},
         sigmaDiffPSave = {}, sigmaDiffTSave = {}, sigmaDDiffSave = {},
         slopeSave ={};
  double sigErr2TotSave = {}, sigErr2NDSave = {},
         sigErr2ELSave = {}, sigErr2INELSave = {},
         sigErr2DiffPSave = {}, sigErr2DiffTSave = {}, sigErr2DDiffSave = {},
         slopeErrSave ={};
  double avNDbSave = {};
  double weightSave = {}, weightSumSave = {}, xSecScaleSave = {};

  // Number of collisions and paricipants. See accessor functions for
  // indices.
  vector<int> nCollSave{}, nProjSave{}, nTargSave{};

  // Map of primary processes and the number of events and the sum of
  // weights.
  map<int,double> sumPrimW{}, sumPrimW2{};
  map<int,int> NPrim{};
  map<int,string> NamePrim{};

  // The info object of the primary process.
  Info primInfo{};

  // Number of failed nucleon excitations.
  int nFailSave = 0;

public:

  // Access to subcollision to be extracted by the user.
  const SubCollisionSet* subCollisionsPtr() { return subCollisionsPtrSave; }

private:

  // Set the subcollision pointer.
  void setSubCollisions(const SubCollisionSet* subCollisionsPtrIn) {
    subCollisionsPtrSave = subCollisionsPtrIn; }

  // Full information about the Glauber calculation, consisting of
  // all subcollisions.
  const SubCollisionSet* subCollisionsPtrSave{};

};

//==========================================================================

// This is the heavy ion user hooks class which in the future may be
// used inside a Pythia object to generate heavy ion collisons. For
// now it is used outside Pythia and requires access to a number of
// Pythia objects.

class HIUserHooks {

public:

  // The default constructor is empty.
  HIUserHooks(): idProjSave(0), idTargSave(0) {}

  // Virtual destructor.
  virtual ~HIUserHooks() {}

  // Initialize this user hook.
  virtual void init(int idProjIn, int idTargIn) {
    idProjSave = idProjIn;
    idTargSave = idTargIn;
  }

  // A user-supplied impact parameter generator.
  virtual bool hasImpactParameterGenerator() const { return false; }
  virtual shared_ptr<ImpactParameterGenerator> impactParameterGenerator()
    const { return nullptr; }

  // A user-supplied NucleusModel for the projectile and target.
  virtual bool hasProjectileModel() const { return false; }
  virtual shared_ptr<NucleusModel> projectileModel() const { return nullptr; }
  virtual bool hasTargetModel() const { return false; }
  virtual shared_ptr<NucleusModel> targetModel() const { return nullptr; }

  // A user-supplied SubCollisionModel for generating nucleon-nucleon
  // subcollisions.
  virtual bool hasSubCollisionModel() { return false; }
  virtual shared_ptr<SubCollisionModel> subCollisionModel() { return nullptr; }

  // A user-supplied ordering of events in (inverse) hardness.
  virtual bool hasEventOrdering() const { return false; }
  virtual double eventOrdering(const Event &, const Info &) { return -1; }

  // A user-supplied method for fixing up proton-neutron mismatch in
  // generated beams.
  virtual bool canFixIsoSpin() const { return false; }
  virtual bool fixIsoSpin(EventInfo &) { return false; }

  // A user-supplied method for shifting the event in impact parameter space.
  virtual bool canShiftEvent() const { return false; }
  virtual EventInfo & shiftEvent(EventInfo & ei) const { return ei; }

  // A user-supplied method of adding a diffractive excitation event
  // to another event, optionally connecting their colours.
  bool canAddNucleonExcitation() const { return false; }
  bool addNucleonExcitation(EventInfo &, EventInfo &, bool) const {
    return false; }

  // A user supplied wrapper around the Pythia::forceHadronLevel()
  virtual bool canForceHadronLevel() const { return false; }
  virtual bool forceHadronLevel(Pythia &) { return false; }

  // A user-supplied way of finding the remnants of an
  // non-diffrcative pp collision (on the target side if tside is
  // true) to be used to give momentum when adding.
  virtual bool canFindRecoilers() const { return false; }
  virtual vector<int>
  findRecoilers(const Event &, bool /* tside */, int /* beam */, int /* end */,
               const Vec4 & /* pdiff */, const Vec4 & /* pbeam */) const {
    return vector<int>();
  }

protected:

  // Information set in the init() function.
  // The PDG id of the projectile and target nuclei.
  int idProjSave, idTargSave;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HIInfo_H
