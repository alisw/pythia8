// SplittingsOnia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Naomi Cooke, Philip Ilten, Leif Lonnblad, Steve Mrenna,
// Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions for fragmentation functions for onia showers.

#include "Pythia8/SplittingsOnia.h"

namespace Pythia8 {

//==========================================================================

// Implementation of the SplitOniaSetup class.
// A helper class derived from OniaSetup used for setting up the
// parton shower splittings for onia.

//--------------------------------------------------------------------------

// The constructor.

SplitOniaSetup::SplitOniaSetup(Info* infoPtrIn, AlphaStrong* alphaSPtrIn,
  int flavourIn) :
  OniaSetup(infoPtrIn, flavourIn, "Shower"), alphaSPtr(alphaSPtrIn) {

  // Set the additional general switch settings.
  onia1S0 = settingsPtr->flag("OniaShower:all(1S0)");

  // Set the names of the additional long-distance matrix-element settings.
  meNames1S0.push_back(cat + ":O(1S0)[1S0(1)]");
  meNames1S0.push_back(cat + ":O(1S0)[3S1(8)]");
  meNames3PJ.push_back(cat + ":O(3PJ)[3P0(1)]");
  meNames3PJ.push_back(cat + ":O(3PJ)[3S1(8)]");

  // Set the names of the splitting settings.
  string hvq((flavour == 4) ? "c" : "b");
  splitNames1S0.push_back(cat + ":" + hvq + "2" + key + "(1S0)[1S0(1)]" + hvq);
  splitNames1S0.push_back(cat + ":g2" + key + "(1S0)[1S0(1)]g");
  splitNames1S0.push_back(cat + ":g2" + key + "(1S0)[3S1(8)]");
  splitNames3S1.push_back(cat + ":" + hvq + "2" + key + "(3S1)[3S1(1)]" + hvq);
  splitNames3S1.push_back(cat + ":g2" + key + "(3S1)[3S1(1)]gg");
  splitNames3S1.push_back(cat + ":g2" + key + "(3S1)[3S1(8)]");
  splitNames3PJ.push_back(cat + ":" + hvq + "2" + key + "(3PJ)[3PJ(1)]" + hvq);
  splitNames3PJ.push_back(cat + ":g2" + key + "(3PJ)[3PJ(1)]g");
  splitNames3PJ.push_back(cat + ":" + hvq + "2" + key + "(3PJ)[3S1(8)]" + hvq);
  splitNames3PJ.push_back(cat + ":g2" + key + "(3PJ)[3S1(8)]");

  // Initialise and check all settings.
  states1S0 = settingsPtr->mvec(cat + ":states(1S0)");
  initStates("(1S0)", states1S0, spins1S0, valid1S0);
  initSettings("(1S0)", states1S0.size(), meNames1S0, mes1S0, valid1S0);
  initSettings("(1S0)", states1S0.size(), splitNames1S0, splits1S0, valid1S0);
  states3S1 = settingsPtr->mvec(cat + ":states(3S1)");
  initStates("(3S1)", states3S1, spins3S1, valid3S1);
  initSettings("(3S1)", states3S1.size(), meNames3S1, mes3S1, valid3S1);
  initSettings("(3S1)", states3S1.size(), splitNames3S1, splits3S1, valid3S1);
  states3PJ = settingsPtr->mvec(cat + ":states(3PJ)");
  initStates("(3PJ)", states3PJ, spins3PJ, valid3PJ);
  initSettings("(3PJ)", states3PJ.size(), meNames3PJ, mes3PJ, valid3PJ);
  initSettings("(3PJ)", states3PJ.size(), splitNames3PJ, splits3PJ, valid3PJ);

}

//--------------------------------------------------------------------------

// Initialise the SplitOnias splitting kernels.

void SplitOniaSetup::setup(vector<SplitOniaPtr> &splits,
  set<double> &thresholds, bool oniaIn) {

  // Initialise the 1S0 splittings.
  int nSplits = splits.size();
  double ldmeFac = settingsPtr->parm("OniaShower:ldmeFac");
  if (valid1S0) {
    for (unsigned int i = 0; i < states1S0.size(); ++i) {
      bool flag = oniaIn || onia || onia1S0 || oniaFlavour;
      if (flag || splits1S0[0][i])
        splits.push_back(make_shared<Split2Q2QQbar1S01Q>
          (flavour, states1S0[i], ldmeFac*mes1S0[0][i], infoPtr, alphaSPtr));
      if (flag || splits1S0[1][i])
        splits.push_back(make_shared<Split2g2QQbar1S01g>
          (states1S0[i], ldmeFac*mes1S0[0][i], infoPtr, alphaSPtr));
      if (flag || splits1S0[2][i])
        splits.push_back(make_shared<Split2g2QQbarX8>
          (states1S0[i], ldmeFac*mes1S0[1][i], 0, mSplit,
            infoPtr, alphaSPtr, thresholds));
    }
  }

  // Initialise the 3S1 splittings.
  if (valid3S1) {
    for (unsigned int i = 0; i < states3S1.size(); ++i) {
      bool flag = oniaIn || onia || onia3S1 || oniaFlavour;
      if (flag || splits3S1[0][i])
        splits.push_back(make_shared<Split2Q2QQbar3S11Q>
          (flavour, states3S1[i], ldmeFac*mes3S1[0][i], infoPtr, alphaSPtr));
      if (flag || splits3S1[1][i])
        splits.push_back(make_shared<Split2g2QQbar3S11gg>
          (states3S1[i], ldmeFac*mes3S1[0][i], infoPtr, alphaSPtr));
      if (flag || splits3S1[2][i])
        splits.push_back(make_shared<Split2g2QQbarX8>
          (states3S1[i], ldmeFac*mes3S1[1][i], 0, mSplit,
            infoPtr, alphaSPtr, thresholds));
    }
  }

  // Initialise the 3PJ splittings.
  if (valid3PJ) {
    for (unsigned int i = 0; i < states3PJ.size(); ++i) {
      bool flag = oniaIn || onia || onia3PJ || oniaFlavour;
      if (flag || splits3PJ[0][i])
        splits.push_back(make_shared<Split2Q2QQbar3PJ1Q>
          (flavour, states3PJ[i], ldmeFac*mes3PJ[0][i], spins3PJ[i],
            infoPtr, alphaSPtr));
      if (flag || splits3PJ[1][i])
        splits.push_back(make_shared<Split2g2QQbar3PJ1g>
          (states3PJ[i], ldmeFac*mes3PJ[0][i], spins3PJ[i],
            infoPtr, alphaSPtr, thresholds));
      if (flag || splits3PJ[2][i])
        splits.push_back(make_shared<Split2Q2QQbar3PJ8Q>
          (flavour, states3PJ[i], ldmeFac*mes3PJ[1][i], spins3PJ[i], mSplit,
            infoPtr, alphaSPtr));
      if (flag || splits3PJ[3][i])
        splits.push_back(make_shared<Split2g2QQbarX8>
          (states3PJ[i], ldmeFac*mes3PJ[1][i], spins3PJ[i], mSplit,
            infoPtr, alphaSPtr, thresholds));
    }
  }

  // Initialise the colour octet splittings.
  onlyOctet = nSplits == (int)splits.size();
  int splitMode = settingsPtr->mode("OniaShower:octetSplit");
  if (splitMode == 2)
    splits.push_back(make_shared<Split2QQbarXq82QQbarX8g>
      (settingsPtr->parm("OniaShower:octetColFac"), infoPtr, alphaSPtr));
  else if (splitMode == 1)
    splits.push_back(make_shared<Split2QQbarXg82QQbarX8g>
      (1., infoPtr, alphaSPtr));

}

//==========================================================================

// Implementation of the SplitOnia class.

//--------------------------------------------------------------------------

// Set up the the overestimate. This calculates the splitting independent
// factor.

double SplitOnia::overestimate(const TimeDipoleEnd &dip, double pT2Min,
    bool enh) {

  zMin = 0.5 - sqrtpos( 0.25 - pT2Min / dip.m2DipCorr );
  zMax = 0.5 + sqrtpos( 0.25 - pT2Min / dip.m2DipCorr );
  if (zMax - zMin < 0) return 0;
  overestimate(dip, pT2Min);
  return cFac*oFac*integrateZ()*(enh ? enhance : 1);

}

//--------------------------------------------------------------------------

// Set branch variables. Call class specific version first, then set
// passed variables.

bool SplitOnia::updateBranchVars(const TimeDipoleEnd* dip, Event& event,
  int &idRadIn, int &idEmtIn, int &colRadIn, int &acolRadIn, int &colEmtIn,
  int &acolEmtIn, int &appendEmtIn, double &pTorigIn, double &pTcorrIn,
  double &pzRadPlusEmtIn, double &pzRadIn, double &pzEmtIn,
  double &mRadIn, double &m2RadIn, double &mEmtIn) {

  if (!kinematics(dip, event)) return false;
  idRadIn        = idRad;
  idEmtIn        = idEmt;
  colRadIn       = colRad;
  acolRadIn      = acolRad;
  colEmtIn       = colEmt;
  acolEmtIn      = acolEmt;
  appendEmtIn    = appendEmt;
  pTorigIn       = pTorig;
  pTcorrIn       = pTcorr;
  pzRadPlusEmtIn = pzRadPlusEmt;
  pzRadIn        = pzRad;
  pzEmtIn        = pzEmt;
  mRadIn         = mRad;
  m2RadIn        = m2Rad;
  mEmtIn         = mEmt;
  return true;

}

//--------------------------------------------------------------------------

// Set the kinematics. Suitable for 1 -> 2 onium emissions.

bool SplitOnia::kinematics(const TimeDipoleEnd* dip, Event& event) {

  idRad          = event[dip->iRadiator].id() > 0 ? idB : -idB;
  idEmt          = idC;
  colRad         = event[dip->iRadiator].col();
  acolRad        = event[dip->iRadiator].acol();
  appendEmt      = 1;
  pTorig         = sqrt(dip->pT2);
  double s       = dip->pT2/(dip->z*(1 - dip->z)) + dip->m2A;
  if (sqrt(s) + dip->mRec >= dip->mDip) return false;
  double m2df    = dip->m2Dip + s - dip->m2Rec;
  double pPlusBC = 0.5*(m2df + sqrt(pow2(m2df) - 4*s*dip->m2Dip))/dip->mDip;
  double pPlusB  = dip->z*pPlusBC;
  double pPlusC  = pPlusBC - pPlusB;
  double pT2corr = s*dip->z*(1 - dip->z) - dip->m2B*(1 - dip->z)
    - dip->m2C*dip->z;
  pTcorr         = sqrt(pT2corr);
  pzRad          = 0.5*(pPlusB - (dip->m2B + pT2corr)/pPlusB);
  pzEmt          = 0.5*(pPlusC - (dip->m2C + pT2corr)/pPlusC);
  pzRadPlusEmt   = pzRad + pzEmt;
  mRad           = sqrt(dip->m2B);
  m2Rad          = dip->m2B;
  mEmt           = sqrt(dip->m2C);
  return true;

}

//--------------------------------------------------------------------------

// Set the colour octet ID and ensure in particle database. Here, the state
// corresponds to 0 is 3S1, 1 is 1S0, and 2 is 3PJ.

void SplitOnia::setOctetID(int state, double mSplit, Info* infoPtr) {

  // Determine quark composition and quantum numbers of the physical state.
  int mod1(10), mod2(1), idHad(idC == 0 ? idB : idC);
  vector<int> digits;
  while (digits.size() < 7) {
    digits.push_back((idHad%mod1 - idHad%mod2) / mod2);
    mod1 *= 10;
    mod2 *= 10;
  }

  // Set the name of the physical and octet state.
  string stateName = "[3S1(8)]";
  if (state == 1) stateName = "[1S0(8)]";
  if (state == 2) stateName = "[3PJ(8)]";

  // Ensure the dummy particle for the colour-octet state is valid.
  ParticleData* particleDataPtr = infoPtr->particleDataPtr;
  int idOct = 9900000 + digits[1]*10000 + state*1000 + digits[5]*100
    + digits[4]*10 + digits[0];
  double m0     = particleDataPtr->m0(idHad) + abs(mSplit);
  double mWidth = 0.0;
  if (!particleDataPtr->isParticle(idOct)) {
    string nameOct    = particleDataPtr->name(idHad) + stateName;
    int    spinType   = state == 1 ? 1 : 3;
    int    chargeType = particleDataPtr->chargeType(idHad);
    int    colType    = 2;
    particleDataPtr->addParticle(idOct, nameOct, spinType, chargeType, colType,
      m0, mWidth, m0, m0);
    ParticleDataEntryPtr entry = particleDataPtr->particleDataEntryPtr(idOct);
    if (entry->id() != 0) entry->addChannel(1, 1.0, 0, idHad, 21);
  } else if (mSplit > 0 && abs(particleDataPtr->m0(idOct) - m0) > 1E-5) {
    particleDataPtr->m0(idOct, m0);
    particleDataPtr->mWidth(idOct, mWidth);
    particleDataPtr->mMin(idOct, m0);
    particleDataPtr->mMax(idOct, m0);
  } else if (particleDataPtr->m0(idOct) <= particleDataPtr->m0(idHad)) {
    loggerPtr->ERROR_MSG("mass of intermediate colour-octet state "
      "increased to be greater than the physical state");
    particleDataPtr->m0(idOct, m0);
    particleDataPtr->mWidth(idOct, mWidth);
    particleDataPtr->mMin(idOct, m0);
    particleDataPtr->mMax(idOct, m0);
  }
  if (idC == 0) {idB = idOct; mB = m0; m2B = pow2(m0);}
  else {idC = idOct; mC = m0; m2C = pow2(m0);}

}

//==========================================================================

// Implementation of the Split2Q2QQbar1S01Q class.
// Splitting class for Q -> QQbar[1S0(1)] Q (Q = c or b).
// The splitting function is taken from equation 19 of Bra93a.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2Q2QQbar1S01Q::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of the z-dependence.
  double zFac = 2.5;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/(8*m2A);

  // Set the overestimate factor.
  oFac = jsFac*zFac;

  // Set the constant prefactor. From Bra93a this is 8/(27*pi*mQ) but
  // the Sudakov generation gives us a factor of 2*pi. LDME units are GeV^3.
  cFac = ldme*16/(27*mA);

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2Q2QQbar1S01Q::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Bra93a.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2Q(m2A), m2O(m2C), s(pT2/(z*(1 - z)) + m2Q);

  // Check kinematic limits.
  if (s <= m2O/z + m2Q/(1 - z)) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than (mO + mQ)^2 - m2Q, i.e. approximately 8m2Q.
  double jFac = s - m2Q;

  // The s-dependent factors, which combined with the Jacobian is
  // always less than (1/8m2Q).
  double sFac = alphaScale(m2O, pT2, s)/pow2(s - m2Q);

  // The z-dependent part of the fragmentation function, including
  // a 1/(s - m2Q)^2 contribution to be dimensionless.
  double zFac = (pow2(s) - 2*m2Q*s - 15*pow2(m2Q)
    - z*(s - m2Q)*(s - pow2(mA + mC)) + 4*s*(s - m2Q)*z*(1 - z)/(2 - z)
    - 4*m2Q*(s - m2Q)*(1 - 3*z)*z/(2 - z) + 4*(pow2((s - m2Q)*z))*(1 - z)
    /(pow2(2 - z)))/pow2(s - m2Q);

  // Return the weight.
  return jFac*sFac*zFac/oFac;

}

//==========================================================================

// Implementation of the Split2g2QQbar1S01g class.
// Splitting class for g -> QQbar[1S0(1)] g (Q = c or b).
// The splitting function is taken from equation 6 of Bra93.
// Note, twice the mass of the heavy quark, 2mQ, has been replaced
// with the onia mass, mC.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2g2QQbar1S01g::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of the z-dependence.
  double zFac = 4;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/m2C;

  // Set the overestimate factor and return the total overestimate.
  oFac = jsFac*zFac;

  // Set the constant prefactor. The Bra93 prefactor is
  // 1/(6*pi*mQ). From the Sudakov generation we have a factor of
  // 2*pi, then we have mQ = mO/2, and a factor of 1/2 for the
  // colour. LDME units are GeV^3.
  cFac = ldme/(3*mC);

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2g2QQbar1S01g::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Bra93.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the onium state and system.
  double m2O(m2C), s(pT2/(z*(1 - z)));

  // Check kinematic limits.
  if (s <= m2O/z) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than mO^2.
  double jFac = s;

  // Calculate equation 6 from Bra93.
  double zFac = pow2(s) + pow2(m2O) - 2*z*(s + m2O)*s + 2*pow2(z*s);
  // Include the factor of 1/(s - 4mQ^2)^2 here to make dimensionless.
  zFac /= pow2(s - m2O);

  // The s-dependent part of the fragmentation function.
  double sFac = alphaScale(m2O, pT2, s)/pow2(s);

  // Return the weight.
  return jFac*sFac*zFac/oFac;

}

//==========================================================================

// Implementation of the Split2Q2QQbar3S11Q class.
// Splitting class for Q -> QQbar[3S1(1)] Q (Q = c or b).
// The splitting function is taken from equations 14 and 15 of Bra93a.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2Q2QQbar3S11Q::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of the z-dependence.
  double zFac = 2.5;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/(8*m2A);

  // Set the overestimate factor and return the total overestimate.
  oFac = jsFac*zFac;

  // Set constant prefactor. From Bra93a this is 8/(27*pi*mQ) but the
  // Sudakov generation gives us a factor of 2*pi. LDME units are GeV^3.
  cFac = ldme*16/(27*mA);

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2Q2QQbar3S11Q::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Bra93a.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2Q(m2A), m2O(m2C), s(pT2/(z*(1 - z)) + m2Q);

  // Check kinematic limits.
  if (s <= m2O/z + m2Q/(1 - z)) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than (mO + mQ)^2 - m2Q, i.e. approximately 8m2Q.
  double jFac = s - m2Q;

  // The s-dependent factors, which combined with the Jacobian is
  // always less than (1/8m2Q).
  double sFac = alphaScale(m2O, pT2, s)/pow2(s - m2Q);

  // The z-dependent part of the fragmentation function, including
  // a 1/(s - m2Q)^2 contribution to be dimensionless.
  double zFac = (pow2(s) - 2*m2Q*s - 47*pow2(m2Q)
    - z*(s - m2Q)*(s - pow2(mA + mC)) + 4*s*(s - m2Q)*z*(1 - z)/(2 - z)
    - 4*m2Q*(s - m2Q)*(8 - 7*z - 5*pow2(z))/(2 - z)
    + 12*(pow2((s - m2Q)*z))*(1 - z)/(pow2(2 - z)))/pow2(s - m2Q);

  // Return the weight.
  return jFac*sFac*zFac/oFac;

}


//==========================================================================

// Implementation of the Split2g2QQbar3S11gg class.
// Splitting class for g -> QQbar[3S1(1)] g g (Q = c or b).
// The splitting function is taken from equations 3 and 8 of Bra95.
// Note, twice the mass of the heavy quark, 2mQ, has been replaced
// with the onia mass, mC.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2g2QQbar3S11gg::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of dimensionless z-dependent parts,
  // including the Jacobian.
  double jzFac = 2.5;

  // Set the overestimate factor and return the total overestimate.
  oFac = pow2(alphaSPtr->alphaS(pT2Min))*jzFac;

  // The constant prefactor. Note Braaten has 5/5184*pi*mc3, but we
  // have mQ = mO/2, 2*pi for the Sudakov generation, and 1/2 for the
  // color. LDME units are GeV^3.
  cFac = ldme*5/(5184*pow3(mC/2));

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2g2QQbar3S11gg::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Bra95.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the onium state and system.
  double m2O(m2C), s(pT2/(z*(1 - z)));

  // Check kinematic limits.
  double m2gg = ygg*s;
  if (s <= m2O/z + m2gg/(1 - z)) return 0;

  // The Jacobian going from dr dy to dpT2/pT2 dygg.
  double jFac = (z*(1 - z))/(2*pT2/m2O);

  // The s-dependent factors including a factor of 2 to make zFac below one.
  double sFac = 2*alphaSPtr->alphaS(pT2)*alphaScale(m2O, pT2, s)
    *z*pow(1 - z, 2 - pg)*pow(ygg, pg);

  // The scaling variables used in the splitting function.
  double r = z*(1 - z)/(pT2/m2O);
  double y = (1 + r - ygg)/2;
  if (2*y >= 1 + r) return 0;
  if (2*z*y <= r + z*z) return 0;

  // Equations 4 - 6 of Bra95.
  double f0 = pow2(r)*(1 + r)*(3 + 12*r + 13*pow2(r))
    - 16*pow2(r)*(1 + r)*(1 + 3*r)*y
    - 2*r*(3 - 9*r - 21*pow2(r) + 7*pow3(r))*pow2(y)
    + 8*r*(4 + 3*r + 3*pow2(r))*pow3(y)
    - 4*r*(9 - 3*r - 4*pow2(r))*pow4(y)
    - 16*(1 + 3*r + 3*pow2(r))*pow5(y)
    + 8*(6 + 7*r)*pow6(y) - 32*pow7(y);
  double f1 = -2*r*(1 + 5*r + 19*pow2(r) + 7*pow3(r))*y
    + 96*pow2(r)*(1 + r)*pow2(y)
    + 8*(1 - 5*r - 22*pow2(r) - 2*pow3(r))*pow3(y)
    + 16*r*(7 + 3*r)*pow4(y)
    - 8*(5 + 7*r)*pow5(y) + 32*pow6(y);
  double f2 = r*(1 + 5*r + 19*pow2(r) + 7*pow3(r))
    - 48*pow2(r)*(1 + r)*y
    - 4*(1 - 5*r - 22*pow2(r) - 2*pow3(r))*pow2(y)
    - 8*r*(7 + 3*r)*pow3(y) + 4*(5 + 7*r)*pow4(y) - 16*pow5(y);

  // Equations 7 - 9 of Bra95.
  double g0 = pow3(r)*(1 - r)*(3 + 24*r + 13*pow2(r))
    - 4*pow3(r)*(7 - 3*r - 12*pow2(r))*y
    - 2*pow3(r)*(17 + 22*r - 7*pow2(r))*pow2(y)
    + 4*pow2(r)*(13 + 5*r - 6*pow2(r))*pow3(y)
    - 8*r*(1 + 2*r + 5*pow2(r) + 2*pow3(r))*pow4(y)
    - 8*r*(3 - 11*r - 6*pow2(r))*pow5(y)
    + 8*(1 - 2*r - 5*pow2(r))*pow6(y);
  double g1 = -2*pow2(r)*(1 + r)*(1 - r)*(1 + 7*r)*y
    + 8*pow2(r)*(1 + 3*r)*(1 - 4*r)*pow2(y)
    + 4*r*(1 + 10*r + 57*pow2(r) + 4*pow3(r))*pow3(y)
    - 8*r*(1 + 29*r + 6*pow2(r))*pow4(y)
    - 8*(1 - 8*r - 5*pow2(r))*pow5(y);
  double g2 = pow2(r)*(1 + r)*(1 - r)*(1 + 7*r)
    - 4*pow2(r)*(1 + 3*r)*(1 - 4*r)*y
    - 2*r*(1 + 10*r + 57*pow2(r) + 4*pow3(r))*pow2(y)
    + 4*r*(1 + 29*r + 6*pow2(r))*pow3(y)
    + 4*(1 - 8*r - 5*pow2(r))*pow4(y);

  // Prefactors.
  double preFac = 1.0/(pow2(1 - y)*pow2(y - r)*pow2(pow2(y) - r));
  double gFac   = (1 + r - 2*y)/(2*(y - r)*sqrt(pow2(y) -r))
    *log(pow2(y - r + sqrt(pow2(y) - r))/(r*(1 + r - 2*y)));

  // The final splitting function.
  double zFac = preFac*(f0 + z*f1 + z*z*f2 + gFac*(g0 + z*g1 + z*z*g2));

  // Return the weight.
  double wt = jFac*sFac*zFac/oFac;
  if (wt > 0 && wt < 0.5) wt = wt > 0.5*rndmPtr->flat() ? 0.5 : 0;
  return wt;

}

//--------------------------------------------------------------------------

// Update the kinematics.

bool Split2g2QQbar3S11gg::kinematics(const TimeDipoleEnd* dip,
  Event &event) {

  idRad          = idB;
  idEmt          = idC;
  colRad         = event[dip->iRadiator].col();
  acolRad        = event[dip->iRadiator].acol();
  colEmt         = 0;
  acolEmt        = 0;
  appendEmt      = 2;
  pTorig         = sqrt(dip->pT2);
  double s       = dip->pT2/(dip->z*(1 - dip->z)) + m2A;
  if (sqrt(s) + dip->mRec >= dip->mDip) return false;
  double m2gg    = dip->m2gg;
  double m2df    = dip->m2Dip + s - dip->m2Rec;
  double pPlusBC = 0.5*(m2df + sqrt(pow2(m2df) - 4*s*dip->m2Dip))/dip->mDip;
  double pPlusB  = dip->z*pPlusBC;
  double pPlusC  = pPlusBC - pPlusB;
  double pT2corr = s*dip->z*(1 - dip->z) - m2gg*(1 - dip->z) - m2C*dip->z;
  pTcorr         = sqrt(pT2corr);
  pzRad          = 0.5*(pPlusB - (m2gg + pT2corr)/pPlusB);
  pzEmt          = 0.5*(pPlusC - (m2C + pT2corr)/pPlusC);
  pzRadPlusEmt   = pzRad + pzEmt;
  mRad           = sqrt(m2gg);
  m2Rad          = m2gg;
  mEmt           = mC;
  return true;

}

//--------------------------------------------------------------------------

// Possibility to add more than one emitted particle, here an extra gluon.

vector<int> Split2g2QQbar3S11gg::addEmitted(Event & e, int iRad, int,
  int iRec, int iDipSel, vector<TimeDipoleEnd> & dipEnds ) {

  // The probability for two gluons with a combined small invariant
  // mass to affect the event differently from a single gluon is
  // negligible.
  double eg = e[iRad].mCalc()/2;
  if (eg <= 0.1) return {};
  int colType = dipEnds[iDipSel].colType;

  // The momenta of the two gluons.
  Vec4 g1(0.0, 0.0,  eg, eg);
  Vec4 g2(0.0, 0.0, -eg, eg);

  // Random rotation.
  RotBstMatrix R;
  R.rot(acos(2.0*rndmPtr->flat() - 1.0), 2*M_PI*rndmPtr->flat());
  R.bst(e[iRad].p());

  // Add in a new colour line.
  int iEm1 = e.copy(iRad, 51);
  int iEm2 = e.copy(iEm1);
  int colNew = e.nextColTag();

  // Transform and set the momenta.
  e[iEm1].p(R*g1);
  e[iEm2].p(R*g2);
  e[iEm1].m(0.0);
  e[iEm2].m(0.0);
  e[iEm1].mothers(iRad);
  e[iEm2].mothers(iRad);
  e[iRad].daughters(iEm1, iEm2);
  e[iRad].status(-57);

  // Set the colors.
  if ( colType > 0 ) {
    e[iEm1].col(colNew);
    e[iEm2].acol(colNew);
  } else {
    e[iEm1].acol(colNew);
    e[iEm2].col(colNew);
  }

  // Now add/modify dipoles.
  vector<TimeDipoleEnd> newDips;
  for (auto &dip : dipEnds) {
    if (dip.colType && dip.iRadiator == iRad) {
      dip.iRadiator = iEm1;
      if (dip.iRecoiler == iRec && dip.colType == colType) {
        newDips.push_back(dip);
        dip.iRecoiler = iEm2;
        dip.isrType = 0;
        newDips.back().iRadiator = iEm2;
      }
    } else if (dip.colType*colType < 0 &&
      dip.iRadiator == iRec && dip.iRecoiler == iRad) {
      newDips.push_back(dip);
      dip.iRecoiler = iEm2;
      newDips.back().iRadiator = iEm2;
      newDips.back().iRecoiler = iEm1;
      newDips.back().colType = -colType;
    } else if (dip.iRecoiler == iRad) dip.iRecoiler = iEm1;
    // TO-DO: weak dipoles should be added here as well.
  }

  // Add the dipoles and return.
  for (auto &dip : newDips) dipEnds.push_back(dip);
  return {iEm1, iEm2};

}

//==========================================================================

// Implementation of the Split2Q2QQbar1P11Q class.
// Splitting class for Q -> QQbar[1P1(1)] Q (Q = c or b).
// The splitting function is taken from equations 16 - 12 of Yua94.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2Q2QQbar1P11Q::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of the z-dependence.
  double zFac = 180;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/(8*m2A);

  // Set the overestimate factor and return the total overestimate.
  oFac = jsFac*zFac;

  // Set the constant prefactor. From Yua94 this is 32/81*r*b^3, but
  // the Sudakov generation gives us a factor of 2*pi. This splitting
  // was derived for a bc state, so we have an additional symmetry
  // factor of 2 for a cc state. LDME units are given as GeV^3, so we
  // additionally divide by mQ^3.
  cFac = ldme*4*M_PI*32/81*r*pow3(b)/pow3(mA);
}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2Q2QQbar1P11Q::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Yua94.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2Q(m2A), m2O(m2C), s(pT2/(z*(1 - z)) + m2Q);

  // Check kinematic limits.
  if (s <= m2O/z + m2Q/(1 - z)) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than (mO + mQ)^2 - m2Q, i.e. approximately 8m2Q.
  double jFac = s - m2Q;

  // Calculate equations 18 - 21 from Yua94.
  vector<double> zFac(4, 0);
  zFac[0] = 64*pow2(r)*pow3(b)*pow4(1 - b*z);
  zFac[1] = 8*r*b*pow3(1 - b*z)*(3 - 2*r - 2*pow2(r)
    - 2*b*(2 + 4*r - pow2(r))*z + pow2(b)*(1 - 2*r)*pow2(z));
  zFac[2] = -pow2(1 - b*z)*(2*(1 - 2*r + 4*pow2(r))
    - (3 - 42*r + 64*pow2(r) - 16*pow3(r))*z
    - 2*r*b*(23 - 14*r - 4*pow2(r))*pow2(z)
    + pow2(b)*(1 + 12*r)*(1 - 2*r)*pow3(z));
  zFac[3] = (1 - z)*(1 - 2*(1 - 2*r)*z + (3 - 2*r + 2*pow2(r))*pow2(z)
    - 2*b*(2 + r - 2*pow2(r))*pow3(z) + pow3(b)*(2 + pow2(r))*pow4(z));

  // The s-dependent part of the fragmentation function.
  // This is equation 23 of Yua94. In the limit s -> 0, this approaches
  // 1/s^2.
  double sFac = 0;
  for (int i = 0; i < 4; ++i)
    sFac += zFac[i]*pow(m2O, 4 - i)/pow(s - pow2(b)*m2O, 5 - i);
  sFac *= alphaScale(m2O, pT2, s)/pow4(1 - b*z);

  // Return the weight.
  return jFac*sFac/oFac;

}

//==========================================================================

// Implementation of the Split2Q2QQbar3PJ1Q class.
// Splitting class for Q -> QQbar[3PJ(1)] Q (Q = c or b).
// The splitting function is taken from equations 16, 23 - 27, 98,
// and 99 of Yua94.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2Q2QQbar3PJ1Q::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of the z-dependence.
  double zFac = 180;
  if (spin == 1) zFac = 120;
  if (spin == 2) zFac = 40;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/(8*m2A);

  // Set the overestimate factor and return the total overestimate.
  oFac = jsFac*zFac;

  // Set the constant prefactor. From Yua94 this is 32/243*r*b^3, but
  // the Sudakov generation gives us a factor of 2*pi. This splitting
  // was derived for a bc state, so we have an additional symmetry
  // factor of 2 for a cc state. LDME units are given as GeV^3, so we
  // additionally divide by mQ^3.
  cFac = ldme*4*M_PI*32/243*r*pow3(b)/pow3(mA);
}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2Q2QQbar3PJ1Q::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Yua94.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2Q(m2A), m2O(m2C), s(pT2/(z*(1 - z)) + m2Q);

  // Check kinematic limits.
  if (s <= m2O/z + m2Q/(1 - z)) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than (mO + mQ)^2 - m2Q, i.e. approximately 8m2Q.
  double jFac = s - m2Q;

  // Calculate equations 24 - 35 from Yua94.
  vector<double> zFac(4, 0);
  if (spin == 0) {
    zFac[0] = 64*pow2(r)*pow3(b)*pow4(1 - b*z);
    zFac[1] = 8*r*b*pow3(1 - b*z)*(1 - 18*r + 14*pow2(r)
      - 2*b*(1 - 2*r + 7*pow2(r))*z + pow2(b)*(1 + 2*r)*pow2(z));
    zFac[2] = -pow2(1 - b*z)*(2*(1 - 4*r)*(1 + 6*r - 4*pow2(r))
      - (5 + 14*r - 8*pow2(r) + 80*pow3(r) - 64*pow4(r))*z
      + 2*b*(2 + 9*r + 18*pow2(r) - 28*pow3(r) - 16*pow4(r))*pow2(z)
      - pow2(b)*(1 + 6*r + 16*pow2(r) - 32*pow3(r))*pow3(z));
    zFac[3] = (1 - z)*pow2(1 - 4*r - (1 - 4*r)*(1 - 2*r)*z
      - r*b*(3 - 4*r)*pow2(z));
  } else if (spin == 1) {
    zFac[0] = 192*pow2(r)*pow3(b)*pow4(1 - b*z);
    zFac[1] = 24*r*b*pow3(1 - b*z)*(2*(1 - r - pow2(r))
      - b*(3 + 10*r - 2*pow2(r))*z + pow2(b)*pow2(z));
    zFac[2] = -6*pow2(1 - b*z)*(2*(1 + 2*r) - (5 - 2*r + 6*pow2(r))*z
      + 2*b*(2 - 3*r - 4*pow2(r))*pow2(z)
      - pow2(b)*(1 - 2*r + 2*pow2(r))*pow3(z));
    zFac[3] = 6*(1 - z)*(1 - 2*(1 - 2*r)*z + (1 - 4*r)*(1 - 2*r)*pow2(z)
      + 2*r*b*(1 - 2*r)*pow3(z) + pow2(r)*pow2(b)*pow4(z));
  } else if (spin == 2) {
    zFac[0] = 320*pow2(r)*pow3(b)*pow4(1 - b*z);
    zFac[1] = 8*r*pow2(b)*pow3(1 - b*z)*(2*(4 + 13*r)
      - (1 + 70*r - 26*pow2(r))*z - b*(7 + 8*r)*pow2(z));
    zFac[2] = -4*pow2(b)*pow2(1 - b*z)*(4*(1 + 4*r)
      - (7 + 12*r - 32*pow2(r))*z
      + 2*(1 + 13*r - 26*pow2(r) + 8*pow3(r))*pow2(z)
      + (1 - 30*r - 5*pow2(r) + 4*pow3(r))*pow3(z));
    zFac[3] = 4*pow2(b)*(1 - z)*(2 - 4*(1 - 2*r)*z
      + (5 - 8*r + 12*pow2(r))*pow2(z)
      - 2*(1 - 2*r)*(3 + 2*pow2(r))*pow3(z)
      + (3 - 12*r + 12*pow2(r) + 2*pow4(r))*pow4(z));
  }

  // The s-dependent part of the fragmentation function.
  // This is equation 23 of Yua94. In the limit s -> 0, this approaches
  // 1/s^2.
  double sFac = 0;
  for (int i = 0; i < 4; ++i)
    sFac += zFac[i]*pow(m2O, 4 - i)/pow(s - pow2(b)*m2O, 5 - i);
  sFac *= alphaScale(m2O, pT2, s)/pow4(1 - b*z);

  // Return the weight.
  return jFac*sFac/oFac;

}

//==========================================================================

// Implementation of the Split2g2QQbar3PJ1g class.
// Splitting class for g -> QQbar[3PJ(1)] g (Q = c or b).
// The splitting function is taken from equations 8 - 12 of Bra94.
// Note, twice the mass of the heavy quark, 2mQ, has been replaced
// with the onia mass, mC.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2g2QQbar3PJ1g::overestimate(
  const TimeDipoleEnd &dip, double pT2Min) {

  // Numerical overestimate of the z-dependence. The overestimate is
  // pT2 dependent, thresholds are used to make the sampling more
  // efficient.
  double zFac(1.5*(2*spin + 1));
  if (dip.pT2 < 3*m2C)    zFac *= 25;
  if (dip.pT2 < 0.26*m2C) zFac *= 250;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/m2C;

  // Set the overestimate factor and return the total overestimate.
  oFac = jsFac*zFac;

  // Set the constant prefactor. From Bra94 this is mQ^2/27, but the Sudakov
  // generation gives us a factor of 2*pi, and we have a color factor
  // of 1/2. LDME units are given as GeV^3, so we additionally divide
  // by mQ^3, where we take mQ = mO/2.
  cFac = ldme*M_PI/(27*mC/2);

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2g2QQbar3PJ1g::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match Bra94.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the onium state and system.
  double m2O(m2C), s(pT2/(z*(1 - z)));

  // Check kinematic limits (extra cut-off for pT2 divergence).
  if (s <= m2O/z || pT2 < 0.3) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than mO^2.
  double jFac = s;

  // Calculate equations 10 - 12 from Bra94.
  double zFac = 0;
  if (spin == 0) {
    zFac = pow2(s - 3*m2O)*(pow2(s - m2O) - 2*(1 - z)*(z*s - m2O)*s);
  } else if (spin == 1) {
    zFac = 6*pow2(s)*(pow2(s - m2O) - 2*(1 - z)*(z*s - m2O)*(s - 2*m2O));
  } else if (spin == 2) {
    zFac = 2*(pow2(s - m2O)*(pow2(s) + 6*pow2(m2O))
      - 2*(1 - z)*(z*s - m2O)*s*(pow2(s) - 6*s*m2O + 6*pow2(m2O)));
  }
  // Include the factor of 1/(s - 4mQ^2)^4 here to make dimensionless.
  zFac /= pow4(s - m2O);

  // The s-dependent part of the fragmentation function.
  double sFac = alphaScale(m2O, pT2, s)/pow2(s);

  // Return the weight.
  return jFac*sFac*zFac/oFac;

}

//==========================================================================

// Implementation of the Split2Q2QQbar3PJ8Q class.
// Splitting class for Q -> QQbar[3PJ(8)] Q (Q = c or b).
// The splitting function is taken from equations 54, 55, and 59 - 61
// of Yua94.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2Q2QQbar3PJ8Q::overestimate(const TimeDipoleEnd &, double pT2Min) {

  // Numerical overestimate of the z-dependence.
  double zFac = 30;

  // Overestimate of s-dependent prefactor and alpha_S (jFac*sFac in weight).
  double jsFac = alphaSPtr->alphaS(pT2Min)/(8*m2A);

  // Set the overestimate factor and return the total overestimate.
  oFac = jsFac*zFac;

  // Set the constant prefactor. From Yua94 this is 1/81*r*b^3, but
  // the Sudakov generation gives us a factor of 2*pi. This splitting
  // was derived for a bc state, so we have an additional symmetry
  // factor of 2 for a cc state. LDME units are given as GeV^3, so we
  // additionally divide by mQ^3. There is also a factor of 2J + 1,
  // see equation 98.
  cFac = (2*spin + 1)*ldme*4*M_PI/81*r*pow3(b)/pow3(mA);

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2Q2QQbar3PJ8Q::weight(const TimeDipoleEnd &dip) const {

  // Redefine z to match the Braaten paper.
  double pT2(dip.pT2), z(1 - zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2Q(m2A), m2O(m2C), s(pT2/(z*(1 - z)) + m2Q);

  // Check kinematic limits.
  if (s <= m2O/z + m2Q/(1 - z)) return 0;

  // The Jacobian going from ds to dpT2/pT2. This expression is
  // always larger than (mO+mQ)^2 - m2Q, i.e. approximately 8m2Q.
  double jFac = s - m2Q;

  // Calculate equations 59 - 61 from Yua94.
  vector<double> zFac(3, 0);
  zFac[0] = -12*r*b*pow2(1 - b*z);
  zFac[1] = -(1 -b*z)*(2*(1 + 2*r) - (1 + 12*r - 4*pow2(r))*z
    -b*(1 + 2*r)*pow2(z));
  zFac[2] = (1 - z)*(1 + 2*r*z + (2 + pow2(r))*pow2(z));

  // The s-dependent part of the fragmentation function.
  // This is equation 55 of Yua94. In the limit s -> 0, this approaches
  // 1/s^2.
  double sFac = 0;
  for (int i = 0; i < 3; ++i)
    sFac += zFac[i]*pow(m2O, 3 - i)/pow(s - pow2(b)*m2O, 4 - i);
  sFac *= alphaScale(m2O, pT2, s)/pow2(1 - b*z);

  // Return the weight.
  return jFac*sFac/oFac;

}

//==========================================================================

// Implementation of the Split2g2QQbarX8 class.
// Splitting class for g -> QQbar[X(8)] (Q = c or b).

//--------------------------------------------------------------------------

// Return the splitting overestimate.

double Split2g2QQbarX8::overestimate(
  const TimeDipoleEnd &dip, double, bool enh) {

  // Define the constant prefactor. From Bra94 this is pi/24, but the
  // Sudakov generation gives us a factor of 2*pi. LDME units are
  // given as GeV^3, so we additionally divide by mQ^3. There is also
  // a 2J + 1 factor for 3PJ states, see equation 23.
  cFac = pow2(M_PI)*(2*spin + 1)*ldme/(12*pow3(mB/2));

  if (dip.pT2 > m2B*(1.0 + delta)) return 1e-20;
  if (dip.pT2 < m2B) return 0;
  double alpha = alphaScale(m2B, dip.pT2, dip.pT2)/(2*M_PI);
  return (-log1p(-alpha*cFac)/(alpha*log1p(delta)))*(enh ? enhance: 1.0);

}

//--------------------------------------------------------------------------

// Update the internal branching variables.

bool Split2g2QQbarX8::kinematics(const TimeDipoleEnd* dip, Event &event) {

  idRad        = idB;
  idEmt        = idC;
  colRad       = event[dip->iRadiator].col();
  acolRad      = event[dip->iRadiator].acol();
  colEmt       = 0;
  acolEmt      = 0;
  appendEmt    = 0;
  pTorig       = sqrt(dip->pT2);
  double s     = m2B;
  double arg   = pow2(dip->m2Dip - s - dip->m2Rec) - 4*s*dip->m2Rec;
  if (arg < 0) return false;
  pzRadPlusEmt = 0.5*sqrtpos(arg)/dip->mDip;
  pTcorr       = 0;
  pzRad        = pzRadPlusEmt;
  pzEmt        = 0;
  mRad         = mB;
  m2Rad        = pow2(mRad);
  mEmt         = 0;
  return true;

}

//==========================================================================

// Implementation of the Split2QQbarXq82QQbarX8g class.
// Splitting class for QQbar[X(8)] -> QQbar[X(8)] + g (Q = c or b),
// treating the colour octet as a quark. The splitting function is
// taken from equation 2.4 of Ska20.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2QQbarXq82QQbarX8g::overestimate(const TimeDipoleEnd &, double) {

  // Set the overestimate factor and constant prefactor.
  oFac = 2;
  cFac = 2./3.*ldme;

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2QQbarXq82QQbarX8g::weight(const TimeDipoleEnd &dip) const {

  // Define pT2 and z.
  double pT2(dip.pT2), z(zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2Q(m2A/4), m2O(m2A), s(pT2/(z*(1 - z)));

  // Check kinematic limits.
  if (s <= m2O/z) return 0;

  // The z-dependent part of the fragmentation function.
  double zFac = (1 + pow2(z))/(1 - z) - 2*m2Q/s;

  // Return the weight (include 1/(1 - z) overestimate).
  return zFac/oFac*(1 - z);

}

//--------------------------------------------------------------------------

// Set the kinematics for colour octet splittings.

bool Split2QQbarXq82QQbarX8g::kinematics(
  const TimeDipoleEnd* dip, Event& event) {
  bool kin = SplitOnia::kinematics(dip, event);
  idRad   = event[dip->iRadiator].id();
  colRad  = event[dip->iRadiator].col();;
  acolRad = event.nextColTag();
  colEmt  = acolRad;
  acolEmt = event[dip->iRadiator].acol();
  return kin;

}

//==========================================================================

// Implementation of the Split2QQbarXq82QQbarX8g class.
// Splitting class for QQbar[X(8)] -> QQbar[X(8)] + g (Q = c or b),
// treating the colour octet as a gluon. The splitting function is
// taken from equation 2.4 of Ska20.

//--------------------------------------------------------------------------

// Return the splitting overestimate.

void Split2QQbarXg82QQbarX8g::overestimate(const TimeDipoleEnd &, double) {

  // Set the overestimate factor and constant prefactor.
  oFac = 2;
  cFac = 3./2.*ldme;

}

//--------------------------------------------------------------------------

// Return the splitting weight.

double Split2QQbarXg82QQbarX8g::weight(const TimeDipoleEnd &dip) const {

  // Define pT2 and z.
  double pT2(dip.pT2), z(zGen);

  // Mass squared of the radiating quark, onium state, and system.
  double m2O(m2A), s(pT2/(z*(1 - z)));

  // Check kinematic limits.
  if (s <= m2O/z) return 0;

  // The z-dependent part of the fragmentation function.
  double zFac = 2*z/(1 - z) - 2*m2O/s + 4./3.*((1 - z)/z + z*(1 - z));

  // Return the weight (include 1/(z*(1-z)) overestimate).
  return zFac/oFac*z*(1 - z);

}

//==========================================================================

} // end namespace Pythia8
