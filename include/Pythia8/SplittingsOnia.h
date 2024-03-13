// SplittingsOnia.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Naomi Cooke, Philip Ilten, Leif Lonnblad, Steve Mrenna,
// Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for fragmentation functions for onia showers.

#ifndef Pythia8_SplittingsOnia_H
#define Pythia8_SplittingsOnia_H

#include "Pythia8/Basics.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/SigmaOnia.h"
#include "Pythia8/SimpleTimeShower.h"

namespace Pythia8 {

//==========================================================================

// A helper class used to set up the onia splittings.  This class
// derives from OniaSetup, where a shared implementation for similar
// configuration between SigmaOniaSetup and this class is
// implemented. This class extends the OniaSetup class for additional
// configuration that is not shared between SigmaOniaSetup and
// SplittingsOniaSetup.

class SplitOniaSetup : public OniaSetup {

public:

  // Constructor.
  SplitOniaSetup(Info* infoPtrIn, AlphaStrong* alphaSPtrIn, int flavourIn);

  // Initialise the splittings.
  void setup(vector<SplitOniaPtr> &splits, set<double> &thresholds,
    bool oniaIn = false);

  // Flag if only octet -> X splittings enabled.
  bool onlyOctet{false};

 private:

  // Additional stored validity and production flags
  // not inherited from OniaSetup.
  bool onia1S0{true}, valid1S0{true};

  // Additional stored pointers. Note, alpha_s cannot be taken from the Info
  // CoupSM pointer since this is not the TimeShower alpha_s.
  AlphaStrong* alphaSPtr{};

  // Stored vectors of settings.
  vector<int> states1S0, spins1S0;
  vector<string> meNames1S0;
  vector< vector<double> > mes1S0;
  vector<string> splitNames1S0, splitNames3S1, splitNames3PJ;
  vector< vector<bool> > splits1S0, splits3S1, splits3PJ;

};

//===========================================================================

// The SplitOnia encodes the generation of a particular
// emission of an onium state from a parton in a dipole.
//
// The basic generation in the shower assumes a function of the form
// dP(pT2, z) = alphaS/2pi dpT2/pT2 dz F(pT2, z)
// in an A -> B + C splitting, where B takes a fraction z of the
// momentum of A.
//
// A subclass of this base clase should encode not only the
// function F(pT2, z), but also an overestimate and the z-integral
// of this overestimate:
//
// G(pT2, z) >= F(pT2, z) and OIZ = \int_0^1 G(pT2, z) dz
//
// The ratio F/G for a generated point should then be used as a weight
// for the phase space point generated. The subclass should also be
// able to generate the z according to the overestimate and possible
// other internal variables that have been integrated out.

class SplitOnia {

public:

  // Constructor that needs to be called by the subclass.
  SplitOnia(int idAIn, int idBIn, int idCIn, double ldmeIn, Info* infoPtrIn,
    AlphaStrong* alphaSPtrIn) :
    idA(idAIn), idB(idBIn), idC(idCIn),
    mA(infoPtrIn->particleDataPtr->m0(idAIn)),
    mB(infoPtrIn->particleDataPtr->m0(idBIn)),
    mC(infoPtrIn->particleDataPtr->m0(idCIn)),
    m2A(pow2(mA)), m2B(pow2(mB)), m2C(pow2(mC)), ldme(ldmeIn),
    alphaMode(infoPtrIn->settingsPtr->mode("OniaShower:alphaScale")),
    loggerPtr(infoPtrIn->loggerPtr), alphaSPtr(alphaSPtrIn),
    rndmPtr(infoPtrIn->rndmPtr) {}

  // Virtual destructor.
  virtual ~SplitOnia() = default;

  // Calculate the splitting overestimate.
  virtual double overestimate(const TimeDipoleEnd &dip, double pT2Min,
    bool enh);

  // Return the weight of the splitting function and Jacobian at the
  // generated point divided by the overestimate.
  virtual double weight(const TimeDipoleEnd &dip) const = 0;

  // Generate internal degrees of freedom (in particluar z).
  // Default behaviour is to generate z from a flat distribution.
  virtual double generateZ(const TimeDipoleEnd &) {
    zGen = zMin + rndmPtr->flat()*(zMax - zMin); return zGen;}

  // Check if this emission is available for the current dipole end
  // with the given ID of the radiator.
  bool isActive(const TimeDipoleEnd &dip, int id, double m) {
    if (idA == 99 && dip.oniumType == 2) {
      idB = id; mA = mB = m; m2A = m2B = pow2(mA);}
    return ldme > 0 && (abs(id) == idA || (idA == 99 && dip.oniumType == 2))
      && dip.mDip > dip.mRec + mB + mC;}

  // Check if the emission is a colour octet state.
  bool isOctet() {return idC == 0;}

  // Check if the emission is enhanced.
  bool isEnhanced() {return enhance != 1;}

  // Set the enhancement factor.
  void setEnhance(const unordered_map<string, double> &enhanceFSRs) {
    enhance = 1;
    auto enhanceFSR = enhanceFSRs.find(enhanceName());
    if ( enhanceFSR != enhanceFSRs.end() && enhanceFSR->second > enhance )
      enhance = enhanceFSR->second;}

  // Return the enhancement weight.
  double enhanceWeight() {return enhance;}

  // Return the name to be used if this emission is enhanced.
  virtual string enhanceName() const {return "";}

  // Possibility to add more than one emitted particle.
  virtual vector<int> addEmitted(Event &, int, int, int, int,
    vector<TimeDipoleEnd> &) {return {};}

  // Update the dipole end with the generated values.
  virtual void updateDipole(TimeDipoleEnd &dip) {
    dip.z = zGen; dip.flavour = idC; dip.mFlavour = mC;
    dip.m2A = m2A; dip.m2B = m2B; dip.m2C = m2C;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
  }

  // Update the passed branching variables with the internal values.
  bool updateBranchVars(const TimeDipoleEnd* dip, Event& event,
    int &idRadIn, int &idEmtIn, int &colRadIn, int &acolRadIn, int &colEmtIn,
    int &acolEmtIn, int &appendEmtIn, double &pTorigIn, double &pTcorrIn,
    double &pzRadPlusEmtIn, double &pzRadIn, double &pzEmtIn,
    double &mRadIn, double &m2RadIn, double &mEmtIn);

protected:

  // Update the internal kinematics.
  virtual bool kinematics(const TimeDipoleEnd* dip, Event &event);

  // Set the z limits and return the difference.
  virtual void overestimate(const TimeDipoleEnd &, double) {}

  // Integrate the distribution used to generate z.
  // Default behaviour is to integrate a flat distribution.
  virtual double integrateZ() const {return zMax - zMin;}

  // Set the colour octet ID and ensure in particle database.
  void setOctetID(int state, double mSplit, Info* infoPtr);

  // The ID and mass (squared) of the radiating particle (A), the
  // scattered radiated particle (B, taking a fraction z of A's
  // momentum), and the emitted particle (C).
  int idA, idB, idC;
  double mA, mB, mC, m2A, m2B, m2C;

  // Enhancement, max value of alpha_S, and long-distance matrix element.
  double enhance{1}, ldme{-1};

  // The constant and variable overestimate prefactors.
  double cFac{0}, oFac{0};

  // The z limits used in the overestimated z-integral and generated z value.
  double zMin{0}, zMax{1}, zGen{0};

  // IDs, colors, and number of emitters to append.
  int idRad{0}, idEmt{0}, colRad{0}, acolRad{0}, colEmt{0}, acolEmt{0},
    appendEmt{1};

  // Transverse and longitudinal momentum, and masses.
  double pTorig{0}, pTcorr{0}, pzRadPlusEmt{0}, pzRad{0}, pzEmt{0}, mRad{0},
    m2Rad{0}, mEmt{0};

  // Mode to evaluate the final alpha_s scale.
  int alphaMode{1};

  // Return alpha_s at the requested scale.
  double alphaScale(double m2, double pT2, double s) const {
    if (alphaMode == 0) return alphaSPtr->alphaS(m2);
    else if (alphaMode == 2) return alphaSPtr->alphaS(s);
    else return alphaSPtr->alphaS(pT2);
  }

  // Logger pointer.
  Logger* loggerPtr{};

  // The alphaS object in current use.
  AlphaStrong* alphaSPtr{};

  // The random number generator in current use.
  Rndm* rndmPtr{};

};

//==========================================================================

// Splitting class for Q -> QQbar[1S0(1)] Q (Q = c or b).
// The splitting function is taken from equation 19 of Bra93a.

class Split2Q2QQbar1S01Q: public SplitOnia {

public:

  // Constructor.
  Split2Q2QQbar1S01Q(int flavourIn, int stateIn, double ldmeIn,
    Info *infoPtrIn, AlphaStrong* alphaSPtrIn) : SplitOnia(flavourIn,
    flavourIn, stateIn, ldmeIn, infoPtrIn, alphaSPtrIn) {}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return ehancement names.
  string enhanceName() const override {return "fsr:Q2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

};

//==========================================================================

// Splitting class for g -> QQbar[1S0(1)] g (Q = c or b).
// The splitting function is taken from equation 6 of Bra93.

class Split2g2QQbar1S01g: public SplitOnia {

public:

  // Constructor.
  Split2g2QQbar1S01g(int stateIn, double ldmeIn, Info *infoPtrIn,
    AlphaStrong* alphaSPtrIn) : SplitOnia(21, 21, stateIn, ldmeIn, infoPtrIn,
    alphaSPtrIn) {}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:G2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

};

//==========================================================================

// Splitting class for Q -> QQbar[3S1(1)] Q (Q = c or b).
// The splitting function is taken from equations 14 and 15 of Bra93a.

class Split2Q2QQbar3S11Q: public SplitOnia {

public:

  // Constructor.
  Split2Q2QQbar3S11Q(int flavourIn, int stateIn, double ldmeIn,
    Info *infoPtrIn, AlphaStrong* alphaSPtrIn) : SplitOnia(flavourIn,
    flavourIn, stateIn, ldmeIn, infoPtrIn, alphaSPtrIn) {}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return enhancement names.
  string enhanceName() const override {return "fsr:Q2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

};

//==========================================================================

// Splitting class for g -> QQbar[3S1(1)] g g (Q = c or b).
// The splitting function is taken from equations 3 and 8 of Bra95.

class Split2g2QQbar3S11gg : public SplitOnia {

public:

  // Constructor.
  Split2g2QQbar3S11gg(int stateIn, double ldmeIn, Info *infoPtrIn,
    AlphaStrong* alphaSPtrIn) : SplitOnia(21, 21, stateIn, ldmeIn,
    infoPtrIn, alphaSPtrIn) {}

  // Generate a 1/(z(1-z)) distribution.
  double generateZ(const TimeDipoleEnd &) override {
    double u(rndmPtr->flat());
    zGen = u < 0.5 ? zMin*pow(zMax/zMin, 2*u):
      1 - (1 - zMax)*pow((1 - zMin)/(1 - zMax), 2*u - 1);
    ygg  = pow(rndmPtr->flat(), 1.0/(1.0 - pg))*zGen;
    return zGen;
  }

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Update the dipole end with the generated values.
  void updateDipole(TimeDipoleEnd &dip) override {
    SplitOnia::updateDipole(dip); dip.m2gg = ygg*dip.pT2/(zGen*(1 - zGen));}

  // Update the internal kinematics.
  bool kinematics(const TimeDipoleEnd* dip, Event &event) override;

  // Possibility to add more than one emitted particle, here an extra gluon.
  vector<int> addEmitted(Event &, int, int, int, int,
    vector<TimeDipoleEnd> &) override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:G2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Integrate a 1/(z(1-z)) distribution.
  double integrateZ() const override {
    return log(zMax/zMin) + log((1 - zMin)/(1 - zMax));}

  // Additional kinematic variables.
  double ygg{0}, pg{0.5};

};

//==========================================================================

// Splitting class for Q -> QQbar[1P1(1)] Q (Q = c or b).
// The splitting function is taken from equations 16 - 21 of Yua94.

class Split2Q2QQbar1P11Q: public SplitOnia {

public:

  // Constructor.
  Split2Q2QQbar1P11Q(int flavourIn, int stateIn, double ldmeIn,
    Info *infoPtrIn, AlphaStrong* alphaSPtrIn) : SplitOnia(flavourIn,
    flavourIn, stateIn, ldmeIn, infoPtrIn, alphaSPtrIn) {}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:Q2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Defined as r = mQ/(mQ + mQ') and b = 1 - r, see Yua93 equation 7.
  double r{0.5}, b{0.5};

};

//==========================================================================

// Splitting class for Q -> QQbar[3PJ(1)] Q (Q = c or b).
// The splitting function is taken from equations 16, 23 - 27, 98,
// and 99 of Yua94.

class Split2Q2QQbar3PJ1Q: public SplitOnia {

public:

  // Constructor.
  Split2Q2QQbar3PJ1Q(int flavourIn, int stateIn, double ldmeIn, int spinIn,
    Info *infoPtrIn, AlphaStrong* alphaSPtrIn) : SplitOnia(flavourIn,
    flavourIn, stateIn, ldmeIn, infoPtrIn, alphaSPtrIn), spin(spinIn) {}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:Q2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Total angular momentum.
  int spin;

  // Defined as r = mQ/(mQ + mQ') and b = 1 - r, see Yua93 equation 7.
  double r{0.5}, b{0.5};

};

//==========================================================================

// Splitting class for g -> QQbar[3PJ(1)] g (Q = c or b).
// The splitting function is taken from equations 8 - 12 of Bra94.

class Split2g2QQbar3PJ1g: public SplitOnia {

public:

  // Constructor.
  Split2g2QQbar3PJ1g(int stateIn, double ldmeIn, int spinIn,
    Info *infoPtrIn, AlphaStrong* alphaSPtrIn, set<double> &thresholds) :
    SplitOnia(21, 21, stateIn, ldmeIn, infoPtrIn, alphaSPtrIn), spin(spinIn) {
    thresholds.insert({0.26*m2C, 3*m2C});
  }

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:G2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Total angular momentum.
  int spin;

};

//==========================================================================

// Splitting class for Q -> QQbar[3PJ(8)] Q (Q = c or b).
// The splitting function is taken from equations 54, 55, and 59 - 61 of
// Yua94.

class Split2Q2QQbar3PJ8Q: public SplitOnia {

public:

  // Constructor.
  Split2Q2QQbar3PJ8Q(int flavourIn, int stateIn, double ldmeIn, int spinIn,
    double mSplitIn, Info *infoPtrIn, AlphaStrong* alphaSPtrIn) : SplitOnia(
    flavourIn, flavourIn, stateIn, ldmeIn, infoPtrIn, alphaSPtrIn),
    spin(spinIn) {setOctetID(0, mSplitIn, infoPtrIn);}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:Q2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Total angular momentum.
  int spin;

  // Defined as r = mQ/(mQ + mQ') and b = 1 - r, see Yua93 equation 7.
  double r{0.5}, b{0.5};

};

//==========================================================================

// Splitting class for g -> QQbar[X(8)] (Q = c or b).

class Split2g2QQbarX8: public SplitOnia {

public:

  // Constructor.
  Split2g2QQbarX8(int stateIn, double ldmeIn, int spinIn, double mSplitIn,
    Info *infoPtrIn, AlphaStrong* alphaSPtrIn, set<double> &thresholds) :
    SplitOnia(21, stateIn, 0, ldmeIn, infoPtrIn, alphaSPtrIn), spin(spinIn) {
    setOctetID(0, mSplitIn, infoPtrIn);
    if (ldme > 0) thresholds.insert({m2B, m2B*(1.0 + delta)});}

  // Return the splitting overestimate.
  double overestimate(const TimeDipoleEnd &dip, double pT2Min,
    bool enh) override;

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override {
    return (dip.pT2 > m2B*(1.0 + delta) || dip.pT2 < m2B) ? 0 : 1;}

  // Update the internal branching variables.
  void updateDipole(TimeDipoleEnd &dip) override {
    dip.z = zGen; dip.flavour = idB; dip.mFlavour = mB;}

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:G2OX";}

private:

  // Update the internal kinematics.
  bool kinematics(const TimeDipoleEnd* dip, Event &event) override;

  // Total angular momentum and delta function width.
  int spin{0};
  double delta{0.01};

};

//==========================================================================

// Splitting class for QQbar[X(8)] -> QQbar[X(8)] + g (Q = c or b),
// treating the colour octet as a quark.

class Split2QQbarXq82QQbarX8g: public SplitOnia {

public:

  // Constructor.
  Split2QQbarXq82QQbarX8g(double colFacIn, Info *infoPtrIn,
    AlphaStrong* alphaSPtrIn) : SplitOnia(99, 99, 21, colFacIn, infoPtrIn,
      alphaSPtrIn) {};

  // Generate a 1/(1 - z) distribution.
  double generateZ(const TimeDipoleEnd &) override {
    zGen = 1 - (1 - zMax)*pow((1 - zMin)/(1 - zMax), rndmPtr->flat());
    return zGen;}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:O2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Update the internal kinematics.
  bool kinematics(const TimeDipoleEnd* dip, Event &event) override;

  // Integrate a 1/(1 - z) distribution.
  double integrateZ() const override {
    return log((1 - zMin)/(1 - zMax));}

};

//==========================================================================

// Splitting class for QQbar[X(8)] -> QQbar[X(8)] + g (Q = c or b),
// treating the colour octet as a heavy gluon.

class Split2QQbarXg82QQbarX8g: public Split2QQbarXq82QQbarX8g {

public:

  // Constructor.
  Split2QQbarXg82QQbarX8g(double colFacIn, Info *infoPtrIn,
    AlphaStrong* alphaSPtrIn) : Split2QQbarXq82QQbarX8g(colFacIn, infoPtrIn,
      alphaSPtrIn) {};

  // Generate a 1/(z*(1 - z)) distribution.
  double generateZ(const TimeDipoleEnd &) override {
    double u = rndmPtr->flat();
    zGen = u < 0.5 ?  zMin*pow(zMax/zMin, 2*u) :
      1 - (1 - zMax)*pow((1 - zMin)/(1 - zMax), 2*u - 1);
    return zGen;}

  // Return the splitting weight.
  double weight(const TimeDipoleEnd &dip) const override;

  // Return the enhancement names.
  string enhanceName() const override {return "fsr:O2OX";}

private:

  // Return the splitting overestimate.
  void overestimate(const TimeDipoleEnd &dip, double pT2Min) override;

  // Integrate a 1/(z*(1 - z)) distribution.
  double integrateZ() const override {
    return log(zMax/zMin) + log((1 - zMin)/(1 - zMax));}

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SplittingsOnia_H
