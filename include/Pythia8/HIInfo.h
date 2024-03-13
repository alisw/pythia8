// HIInfo.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the HIInfo, EventInfo and
// HIUserHooks classes, as well as the HIUnits namespace.
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
  HIInfo()
    : idProjSave(0), idTargSave(0), bSave(0.0), NSave(0), NAccSave(0),
      sigmaTotSave(0.0), sigmaNDSave(0.0), sigErr2TotSave(0.0),
      sigErr2NDSave(0.0), avNDbSave(0.0), weightSave(0.0), weightSumSave(0.0),
      nCollSave(10, 0), nProjSave(10, 0), nTargSave(10, 0), nFailSave(0),
      subCollisionsPtrSave(nullptr) {}

  // The impact-parameter distance in the current event.
  double b() const {
    return bSave;
  }

  // The impact parameter angle.
  double phi() const {
    return phiSave;
  }

  // The Monte Carlo integrated total cross section in the current run.
  double sigmaTot() const {
    return sigmaTotSave/millibarn;
  }

  // The estimated statistical error on sigmaTot().
  double sigmaTotErr() const {
    return sqrt(sigErr2TotSave/max(1.0, double(NSave)))/millibarn;
  }

  // The Monte Carlo integrated non-diffractive cross section in the
  // current run.
  double sigmaND() const {
    return sigmaNDSave/millibarn;
  }

  // The estimated statistical error on sigmaND().
  double sigmaNDErr() const {
    return sqrt(sigErr2NDSave/max(1.0, double(NSave)));
  }

  // The average NN non-diffractive impact parameter to be used to
  // communicate to Pythia's MPI machinery.
  double avNDb() const {
    return avNDbSave;
  }

  // The number of attempted impact parameter points.
  long nAttempts() const {
    return NSave;
  }

  // The number of produced events.
  long nAccepted() const {
    return NAccSave;
  }

  // The total number of separate sub-collisions.
  int nCollTot() const { return nCollSave[0]; }

  // The number of separate non-diffractive sub collisions in the
  // current event.
  int nCollND() const { return nCollSave[1]; }

  // The total number of non-diffractive sub collisions in the current event.
  int nCollNDTot() const { return nProjSave[1] + nTargSave[1] - nCollSave[1]; }

  // The number of separate single diffractive projectile excitation
  // sub collisions in the current event.
  int nCollSDP() const { return nCollSave[2]; }

  // The number of separate single diffractive target excitation sub
  // collisions in the current event.
  int nCollSDT() const { return nCollSave[3]; }

  // The number of separate double diffractive sub collisions in the
  // current event.
  int nCollDD() const { return nCollSave[4]; }

  // The number of separate central diffractive sub collisions in the
  // current event.
  int nCollCD() const { return nCollSave[5]; }

  // The number of separate elastic sub collisions.
  int nCollEL() const { return nCollSave[6]; }

  // The number of interacting projectile nucleons in the current event.
  int nPartProj() const { return nProjSave[0]; }

  // The number of absorptively wounded projectile nucleons in the
  // current event.
  int nAbsProj() const { return nProjSave[1]; }

  // The number of diffractively wounded projectile nucleons in the
  // current event.
  int nDiffProj() const { return nProjSave[2]; }

  // The number of elastically scattered projectile nucleons in the
  // current event.
  int nElProj() const { return nProjSave[3]; }

  // The number of interacting projectile nucleons in the current
  // event.
  int nPartTarg() const { return nTargSave[0]; }

  // The number of absorptively wounded projectile nucleons in the
  // current event.
  int nAbsTarg() const { return nTargSave[1]; }

  // The number of diffractively wounded projectile nucleons in the
  // current event.
  int nDiffTarg() const { return nTargSave[2]; }

  // The number of elastically scattered projectile nucleons in the
  // current event.
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
  void failedExcitation() {
    ++nFailSave;
  }

private:

  // Register a tried impact parameter point giving the total elastic
  // amplitude, the impact parameter and impact parameter generation weight.
  void addAttempt(double T, double bin, double phiin, double bweight);

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

  // Register a full sub collision.
  int addSubCollision(const SubCollision & c);

  // Register a participating projectile/target nucleon.
  int addProjectileNucleon(const Nucleon & n);
  int addTargetNucleon(const Nucleon & n);

  // Id of the colliding nuclei.
  int idProjSave, idTargSave;

  // Impact parameter.
  double bSave;
  double phiSave;

  // Cross section estimates.
  long NSave, NAccSave;
  double sigmaTotSave, sigmaNDSave, sigErr2TotSave, sigErr2NDSave;
  double avNDbSave;
  double weightSave, weightSumSave;

  // Number of collisions and paricipants. See accessor functions for
  // indices.
  vector<int> nCollSave, nProjSave, nTargSave;

  // Map of primary processes and the number of events and the sum of
  // weights.
  map<int,double> sumPrimW, sumPrimW2;
  map<int,int> NPrim;
  map<int,string> NamePrim;

  // The info object of the primary process.
  Info primInfo;

  // Number of failed nucleon excitations.
  int nFailSave;

public:

  // Access to subcollision to be extracted by the user.
  const SubCollisionSet* subCollisionsPtr() { return subCollisionsPtrSave; }

private:

  // Set the subcollision pointer.
  void setSubCollisions(const SubCollisionSet* subCollisionsPtrIn) {
    subCollisionsPtrSave = subCollisionsPtrIn; }

  // Full information about the Glauber calculation, consisting of
  // all subcollisions.
  const SubCollisionSet* subCollisionsPtrSave;

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
