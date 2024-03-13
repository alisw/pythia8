// StringFragmentation.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the StringEnd and
// StringFragmentation classes.

#include "Pythia8/StringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The StringEnd class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.

// Avoid unphysical solutions to equation system.
const double StringEnd::TINY = 1e-6;

// Assume two (eX, eY) regions are related if pT2 differs by less.
const double StringEnd::PT2SAME = 0.01;

// Fictitious typical mass and pT used to look ahead for approximately where
// the next hadron will be produced, to quantify string close-packing there.
const double StringEnd::MEANMMIN = 0.2;
const double StringEnd::MEANM    = 0.8;
const double StringEnd::MEANPT   = 0.4;

//--------------------------------------------------------------------------

// Set up initial endpoint values from input.

void StringEnd::setUp(bool fromPosIn, int iEndIn, int idOldIn, int iMaxIn,
  double pxIn, double pyIn, double GammaIn, double xPosIn, double xNegIn,
  int colIn) {

  // Simple transcription from input.
  fromPos    = fromPosIn;
  iEnd       = iEndIn;
  iMax       = iMaxIn;
  flavSelNow = *flavSelPtr;
  flavOld    = FlavContainer(idOldIn);
  pxOld      = pxIn;
  pyOld      = pyIn;
  GammaOld   = GammaIn;
  iPosOld    = (fromPos) ? 0 : iMax;
  iNegOld    = (fromPos) ? iMax : 0;
  xPosOld    = xPosIn;
  xNegOld    = xNegIn;
  colOld     = colIn;

}

//--------------------------------------------------------------------------

// Fragment off one hadron from the string system, in flavour and pT.

void StringEnd::newHadron(double kappaRatio, bool forbidPopcornNow,
  bool allowPop, double strangeFac, double probQQmod) {

  // Avoid contradiction in allowPop and forbidPopcornNow.
  if (!allowPop) {
    flavOld.nPop = 0;
    forbidPopcornNow = false;
  }

  // In case we are using the thermal model or Gaussian with
  // mT2 suppression we have to pick the pT first.
  if (thermalModel || mT2suppression) {

    // Pick its transverse momentum.
    pair<double, double> pxy = pTSelPtr->pxy(flavNew.id, kappaRatio);
    pxNew = pxy.first;
    pyNew = pxy.second;
    pxHad = pxOld + pxNew;
    pyHad = pyOld + pyNew;
    double pT2Had = pow2(pxHad) + pow2(pyHad);

    // Pick new flavour and form a new hadron.
    // For forbidPopcornNow == true it must be a baryon.
    do {
      flavNew = flavSelPtr->pick( flavOld, sqrt(pT2Had), kappaRatio);
      idHad   = flavSelPtr->getHadronID( flavOld, flavNew);
    } while (idHad == 0 || (forbidPopcornNow && (abs(idHad)/1000)%10 == 0));

    // Get its mass and thereby define its transverse mass.
    mHad   = flavSelPtr->getHadronMassWin(idHad);
    mT2Had = pow2(mHad) + pow2(pxHad) + pow2(pyHad);
  }

  // In case of the Gaussian without mT2 suppression we pick
  // the new flavour first to make the width flavour dependent.
  else {

    // Pick new flavour and form a new hadron.
    // For forbidPopcornNow == true it must be a baryon.
    do {
      // Reinitialise probabilities if close-packing.
      if ( (closePacking && (probQQmod < 1. || kappaRatio > 1.))
        || strangeFac > 0. ) {
        flavSelNow.init(kappaRatio, strangeFac, probQQmod);
        flavNew = flavSelNow.pick( flavOld, -1.0, 0.0, allowPop );
        idHad   = flavSelNow.combine( flavOld, flavNew);
      } else {
        flavNew = flavSelPtr->pick( flavOld );
        idHad   = flavSelPtr->combine( flavOld, flavNew);
      }
    } while (idHad == 0 || (forbidPopcornNow && (abs(idHad)/1000)%10 == 0)
      || ((abs(flavNew.id) > 10) && !allowPop) );

    // Pick its transverse momentum.
    pair<double, double> pxy = pTSelPtr->pxy(flavNew.id, kappaRatio);
    pxNew = pxy.first;
    pyNew = pxy.second;
    pxHad = pxOld + pxNew;
    pyHad = pyOld + pyNew;

    // Pick its mass and thereby define its transverse mass.
    mHad   = particleDataPtr->mSel(idHad);
    mT2Had = pow2(mHad) + pow2(pxHad) + pow2(pyHad);
  }

}

//--------------------------------------------------------------------------

// Make (anti)baryon from stepping over pearl (anti)quark.

void StringEnd::pearlHadron(StringSystem& system, int idPearlIn,
  Vec4 pPearlIn) {

  // Find transverse momentum contribution of pearl in current region.
  StringRegion& region = system.region( iPosOld, iNegOld);
  region.project( pPearlIn);
  double pxPearl   = region.px();
  double pyPearl   = region.py();

  // Define baryon including pearl.
  flavOld.id = flavSelPtr->makeDiquark(idPearlIn, -flavOld.id);
  idHad = 0;
  do {
    idHad = flavSelPtr->combine( flavOld, flavNew);
  } while (idHad == 0 || abs(idHad) > 10000);
  mHad   = particleDataPtr->mSel(idHad);
  mT2Had = pow2(mHad) + pow2(pxHad + pxPearl) + pow2(pyHad + pyPearl);
}

//--------------------------------------------------------------------------

// Fragment off one hadron from the string system, in momentum space,
// by taking steps from positive end.

Vec4 StringEnd::kinematicsHadron( StringSystem& system,
  StringVertex& newVertex, bool useInputZ, double zHadIn,
  bool pearlIn, Vec4 pPearlIn) {

  // Ensure pearl momentum is not used if pearlIn = false.
  if (pearlIn) pPearlIn = {0., 0., 0., 0.};

  // Pick fragmentation step z and calculate new Gamma.
  if (useInputZ) zHad = zHadIn;
  else zHad = zSelPtr->zFrag( flavOld.id, flavNew.id, mT2Had);
  GammaNew = (1. - zHad) * (GammaOld + mT2Had / zHad);

  // Set up references that are direction-neutral;
  // ...Dir for direction of iteration and ...Inv for its inverse.
  int&    iDirOld = (fromPos) ? iPosOld : iNegOld;
  int&    iInvOld = (fromPos) ? iNegOld : iPosOld;
  int&    iDirNew = (fromPos) ? iPosNew : iNegNew;
  int&    iInvNew = (fromPos) ? iNegNew : iPosNew;
  double& xDirOld = (fromPos) ? xPosOld : xNegOld;
  double& xInvOld = (fromPos) ? xNegOld : xPosOld;
  double& xDirNew = (fromPos) ? xPosNew : xNegNew;
  double& xInvNew = (fromPos) ? xNegNew : xPosNew;
  double& xDirHad = (fromPos) ? xPosHad : xNegHad;
  double& xInvHad = (fromPos) ? xNegHad : xPosHad;

  // Start search for new breakup in the old region.
  iDirNew = iDirOld;
  iInvNew = iInvOld;
  Vec4 pTNew;

  // Each step corresponds to trying a new string region.
  for (int iStep = 0; ; ++iStep) {

    // Referance to current string region.
    StringRegion& region = system.region( iPosNew, iNegNew);
    colNew = fromPos ? region.colPos : region.colNeg;

    // Now begin special section for rapid processing of low region.
    if (iStep == 0 && iPosOld + iNegOld == iMax && !pearlIn) {

      // A first step within a low region is easy.
      if (mT2Had < zHad * xDirOld * (1. - xInvOld) * region.w2) {

        // Translate into x coordinates.
        xDirHad = zHad * xDirOld;
        xInvHad = mT2Had / (xDirHad * region.w2);
        xDirNew = xDirOld - xDirHad;
        xInvNew = xInvOld + xInvHad;

        // Store breakup vertex information from the fragmentation process.
        newVertex.store( fromPos, iPosNew, iNegNew, xPosNew, xNegNew);

        // Find and return four-momentum of the produced particle.
        return region.pHad( xPosHad, xNegHad, pxHad, pyHad);

      // A first step out of a low region also OK, if there are more regions.
      // Negative energy signals failure, i.e. in last region.
      } else {
        --iInvNew;
        if (iInvNew < 0) return Vec4(0., 0., 0., -1.);

        // Momentum taken by stepping out of region. Continue to next region.
        xInvHad = 1. - xInvOld;
        xDirHad = 0.;
        pSoFar  = region.pHad( xPosHad, xNegHad, pxOld, pyOld);
        continue;
      }

    // Else, for first step, take into account starting pT.
    } else if (iStep == 0) {
      pSoFar  = region.pHad( 0., 0., pxOld, pyOld);
      pSoFar += pPearlIn;
      pTNew   = region.pHad( 0., 0., pxNew, pyNew);
    }


    // Now begin normal treatment of nontrivial regions.
    // Set up four-vectors in a region not visited before.
    if (!region.isSetUp) region.setUp(
      system.regionLowPos(iPosNew).pPos,
      system.regionLowNeg(iNegNew).pNeg,
      system.regionLowPos(iPosNew).colPos,
      system.regionLowNeg(iNegNew).colNeg, true);

    // If new region is vanishingly small, continue immediately to next.
    // Negative energy signals failure to do this, i.e. moved too low.
    if (region.isEmpty) {
      xDirHad = (iDirNew == iDirOld) ? xDirOld : 1.;
      xInvHad = 0.;
      pSoFar += region.pHad( xPosHad, xNegHad, 0., 0.);
      ++iDirNew;
      if (iDirNew + iInvNew > iMax) return Vec4(0., 0., 0., -1.);
      continue;
    }

    // Reexpress pTNew w.r.t. base vectors in new region, if possible.
    // Recall minus sign from normalization e_x^2 = e_y^2 = -1.
    double pxNewTemp = -pTNew * region.eX;
    double pyNewTemp = -pTNew * region.eY;
    if (abs( pxNewTemp * pxNewTemp + pyNewTemp * pyNewTemp
      - pxNew * pxNew - pyNew * pyNew) < PT2SAME) {
      pxNew = pxNewTemp;
      pyNew = pyNewTemp;
    }

    // Four-momentum taken so far, including new pT.
    Vec4 pTemp = pSoFar + region.pHad( 0., 0., pxNew, pyNew);

    // Derive coefficients for m2 expression.
    // cM2 * x+ + cM3 * x- + cM4 * x+ * x- = m^2 - cM1;
    double cM1 = pTemp.m2Calc();
    double cM2 = 2. * (pTemp * region.pPos);
    double cM3 = 2. * (pTemp * region.pNeg);
    double cM4 = region.w2;
    if (!fromPos) swap( cM2, cM3);

    // Derive coefficients for Gamma expression.
    // cGam2 * x+ + cGam3 * x- + cGam4 * x+ * x- = Gamma_new - cGam1;
    double cGam1 = 0.;
    double cGam2 = 0.;
    double cGam3 = 0.;
    double cGam4 = 0.;
    for (int iInv = iInvNew; iInv <= iMax - iDirNew; ++iInv) {
      double xInv = 1.;
      if (iInv == iInvNew) xInv = (iInvNew == iInvOld) ? xInvOld : 0.;
      for (int iDir = iDirNew; iDir <= iMax - iInv; ++iDir) {
        double xDir = (iDir == iDirOld) ? xDirOld : 1.;
        int iPos = (fromPos) ? iDir : iInv;
        int iNeg = (fromPos) ? iInv : iDir;
        StringRegion& regionGam =  system.region( iPos, iNeg);
        if (!regionGam.isSetUp) regionGam.setUp(
          system.regionLowPos(iPos).pPos,
          system.regionLowNeg(iNeg).pNeg,
          system.regionLowPos(iPos).colPos,
          system.regionLowNeg(iNeg).colNeg, true);
        double w2 = regionGam.w2;
        cGam1 += xDir * xInv * w2;
        if (iDir == iDirNew) cGam2 -= xInv * w2;
        if (iInv == iInvNew) cGam3 += xDir * w2;
        if (iDir == iDirNew && iInv == iInvNew) cGam4 -= w2;
      }
    }

    // Solve (m2, Gamma) equation system => r2 * x-^2 + r1 * x- + r0 = 0.
    double cM0   = pow2(mHad) - cM1;
    double cGam0 = GammaNew - cGam1;
    double r2    = cM3 * cGam4 - cM4 * cGam3;
    double r1    = cM4 * cGam0 - cM0 * cGam4 + cM3 * cGam2 - cM2 * cGam3;
    double r0    = cM2 * cGam0 - cM0 * cGam2;
    double root  = sqrtpos( r1*r1 - 4. * r2 * r0 );
    if (abs(r2) < TINY || root < TINY) return Vec4(0., 0., 0., -1.);
    xInvHad      = 0.5 * (root / abs(r2) - r1 / r2);
    if (abs(cM2 + cM4 * xInvHad) < TINY) return Vec4(0., 0., 0., -1.);
    xDirHad      = (cM0 - cM3 * xInvHad) / (cM2 + cM4 * xInvHad);

    // Define position of new trial vertex.
    xDirNew = (iDirNew == iDirOld) ? xDirOld - xDirHad : 1. - xDirHad;
    xInvNew = (iInvNew == iInvOld) ? xInvOld + xInvHad : xInvHad;

    // Step up to new region if new x- > 1.
    if (xInvNew > 1.) {
      xInvHad = (iInvNew == iInvOld) ? 1. - xInvOld : 1.;
      xDirHad = 0.;
      pSoFar += region.pHad( xPosHad, xNegHad, 0., 0.);
      --iInvNew;
      if (iInvNew < 0) return Vec4(0., 0., 0., -1.);
      continue;

    // Step down to new region if new x+ < 0.
    } else if (xDirNew < 0.) {
      xDirHad = (iDirNew == iDirOld) ? xDirOld : 1.;
      xInvHad = 0.;
      pSoFar += region.pHad( xPosHad, xNegHad, 0., 0.);
      ++iDirNew;
      if (iDirNew + iInvNew > iMax) return Vec4(0., 0., 0., -1.);
      continue;
    }

    // Store breakup vertex information from the fragmentation process.
    newVertex.store( fromPos, iPosNew, iNegNew, xPosNew, xNegNew);

    // Else we have found the correct region, and can set the new
    // colour index and return four-momentum.
    colNew = fromPos ? region.colPos : region.colNeg;
    if (useInputZ) return Vec4( 0., 0., 0., 0.);
    else return pSoFar + region.pHad( xPosHad, xNegHad, pxNew, pyNew);

  // End of "infinite" loop of stepping to new region.
  }

}

//--------------------------------------------------------------------------

// Generate momentum for some possible next hadron, based on mean values
// to get an estimate for rapidity and pT.

Vec4 StringEnd::kinematicsHadronTmp( StringSystem system, Vec4 pRem,
  double phi, double mult) {

  // Now estimate the energy the next hadron will take.
  double mRem     = pRem.mCalc();
  double meanM    = (mRem > 0.0) ? max( MEANMMIN, min( MEANM, mRem) ) : MEANM;
  double meanMT2  = pow2(meanM) + pow2(MEANPT);
  double GammaNow = (1.0 + aLund) / bLund;
  // Modify Gamma value in case of earlier fails.
  if (mult > 0.0) GammaNow *= mult;
  double zPlus, zMinus;
  if (GammaOld > 1e-10) {
    double tmp = ( GammaNow + meanMT2 - GammaOld ) / GammaOld;
    zPlus      = (-0.5 * tmp + sqrt(0.25 * pow2(tmp) + meanMT2 / GammaOld));
    zMinus     = (-0.5 * tmp - sqrt(0.25 * pow2(tmp) + meanMT2 / GammaOld));
  // Special case of first hadron.
  } else {
    zPlus  = GammaNow / (GammaNow + meanMT2);
    zMinus = -1.0;
  }
  bool zPlusOk    = (zPlus < 1.0) && (zPlus > 0.0);
  bool zMinusOk   = (zMinus < 1.0) && (zMinus > 0.0);
  // Negative energy signals failure.
  if ( (!zPlusOk) && (!zMinusOk) ) return Vec4(0., 0., 0., -1.);
  double zHadTmp  = (zPlusOk ? zPlus : zMinus);
  double pxHadTmp = cos(phi) * MEANPT;
  double pyHadTmp = sin(phi) * MEANPT;

  // First make a copy of all variables to not overwrite anything.
  int    iPosOldTmp = iPosOld, iNegOldTmp = iNegOld;
  int    iPosNewTmp = iPosNew, iNegNewTmp = iNegNew;
  double xPosOldTmp = xPosOld, xNegOldTmp = xNegOld;
  double xPosNewTmp = xPosNew, xNegNewTmp = xNegNew;
  double xPosHadTmp = xPosHad, xNegHadTmp = xNegHad;
  double pxNewTmp   = pxNew,   pxOldTmp   = pxOld;
  double pyNewTmp   = pyNew,   pyOldTmp   = pyOld;
  Vec4   pSoFarTmp  = pSoFar;

  // Set up references that are direction-neutral;
  // ...Dir for direction of iteration and ...Inv for its inverse.
  int&    iDirOld = (fromPos) ? iPosOldTmp : iNegOldTmp;
  int&    iInvOld = (fromPos) ? iNegOldTmp : iPosOldTmp;
  int&    iDirNew = (fromPos) ? iPosNewTmp : iNegNewTmp;
  int&    iInvNew = (fromPos) ? iNegNewTmp : iPosNewTmp;
  double& xDirOld = (fromPos) ? xPosOldTmp : xNegOldTmp;
  double& xInvOld = (fromPos) ? xNegOldTmp : xPosOldTmp;
  double& xDirNew = (fromPos) ? xPosNewTmp : xNegNewTmp;
  double& xInvNew = (fromPos) ? xNegNewTmp : xPosNewTmp;
  double& xDirHad = (fromPos) ? xPosHadTmp : xNegHadTmp;
  double& xInvHad = (fromPos) ? xNegHadTmp : xPosHadTmp;

  // Start search for new breakup in the old region.
  iDirNew = iDirOld;
  iInvNew = iInvOld;
  Vec4 pTNew;

  // Each step corresponds to trying a new string region.
  for (int iStep = 0; ; ++iStep) {

    // Referance to current string region.
    StringRegion region = system.region( iPosNewTmp, iNegNewTmp);

    // Now begin special section for rapid processing of low region.
    if (iStep == 0 && iPosOldTmp + iNegOldTmp == iMax) {

      // A first step within a low region is easy. Make sure we use this
      // region in case it's the last one.
      if ( (meanMT2 < zHadTmp * xDirOld * (1. - xInvOld) * region.w2)
        || (iInvNew < 1) ) {
        if (iInvNew < 1)
          zHadTmp = meanMT2 / (xDirOld * (1. - xInvOld) * region.w2);

        // Translate into x coordinates.
        xDirHad = zHadTmp * xDirOld;
        xInvHad = meanMT2 / (xDirHad * region.w2);
        xDirNew = xDirOld - xDirHad;
        xInvNew = xInvOld + xInvHad;

        // Find and return four-momentum of the produced particle.
        return region.pHad( xPosHadTmp, xNegHadTmp, pxHadTmp, pyHadTmp);

      // A first step out of a low region also OK, if there are more regions.
      // Negative energy signals failure, i.e. in last region.
      } else {
        --iInvNew;
        // Should be covered by the above check.
        if (iInvNew < 0) return Vec4(0., 0., 0., -1.);

        // Momentum taken by stepping out of region. Continue to next region.
        xInvHad   = 1. - xInvOld;
        xDirHad   = 0.;
        pSoFarTmp = region.pHad( xPosHadTmp, xNegHadTmp, pxOldTmp, pyOldTmp);
        continue;
      }

    // Else, for first step, take into account starting pT.
    } else if (iStep == 0) {
      pSoFarTmp = region.pHad( 0., 0., pxOldTmp, pyOldTmp);
      pTNew     = region.pHad( 0., 0., pxNewTmp, pyNewTmp);
    }

    // Now begin normal treatment of nontrivial regions.
    // Set up four-vectors in a region not visited before.
    if (!region.isSetUp) region.setUp(
      system.regionLowPos(iPosNewTmp).pPos,
      system.regionLowNeg(iNegNewTmp).pNeg,
      system.regionLowPos(iPosNewTmp).colPos,
      system.regionLowNeg(iNegNewTmp).colNeg, true);

    // If new region is vanishingly small, continue immediately to next.
    // Negative energy signals failure to do this, i.e. moved too low.
    if (region.isEmpty) {
      xDirHad = (iDirNew == iDirOld) ? xDirOld : 1.;
      xInvHad = 0.;
      pSoFarTmp += region.pHad( xPosHadTmp, xNegHadTmp, 0., 0.);
      ++iDirNew;
      if (iDirNew + iInvNew > iMax) return Vec4(0., 0., 0., -1.);
      continue;
    }

    // Reexpress pTNew w.r.t. base vectors in new region, if possible.
    // Recall minus sign from normalization e_x^2 = e_y^2 = -1.
    double pxNewRegNow = -pTNew * region.eX;
    double pyNewRegNow = -pTNew * region.eY;
    if (abs( pxNewRegNow * pxNewRegNow + pyNewRegNow * pyNewRegNow
      - pxNewTmp * pxNewTmp - pyNewTmp * pyNewTmp) < PT2SAME) {
      pxNewTmp = pxNewRegNow;
      pyNewTmp = pyNewRegNow;
    }

    // Four-momentum taken so far, including new pT.
    Vec4 pTemp = pSoFarTmp + region.pHad( 0., 0., pxNewTmp, pyNewTmp);

    // Derive coefficients for m2 expression.
    // cM2 * x+ + cM3 * x- + cM4 * x+ * x- = m^2 - cM1;
    double cM1 = pTemp.m2Calc();
    double cM2 = 2. * (pTemp * region.pPos);
    double cM3 = 2. * (pTemp * region.pNeg);
    double cM4 = region.w2;
    if (!fromPos) swap( cM2, cM3);

    // Derive coefficients for Gamma expression.
    // cGam2 * x+ + cGam3 * x- + cGam4 * x+ * x- = Gamma_new - cGam1;
    double cGam1 = 0.;
    double cGam2 = 0.;
    double cGam3 = 0.;
    double cGam4 = 0.;
    for (int iInv = iInvNew; iInv <= iMax - iDirNew; ++iInv) {
      double xInv = 1.;
      if (iInv == iInvNew) xInv = (iInvNew == iInvOld) ? xInvOld : 0.;
      for (int iDir = iDirNew; iDir <= iMax - iInv; ++iDir) {
        double xDir = (iDir == iDirOld) ? xDirOld : 1.;
        int iPos = (fromPos) ? iDir : iInv;
        int iNeg = (fromPos) ? iInv : iDir;
        StringRegion regionGam =  system.region( iPos, iNeg);
        if (!regionGam.isSetUp) regionGam.setUp(
          system.regionLowPos(iPos).pPos,
          system.regionLowNeg(iNeg).pNeg,
          system.regionLowPos(iPos).colPos,
          system.regionLowNeg(iNeg).colNeg, true);
        double w2 = regionGam.w2;
        cGam1 += xDir * xInv * w2;
        if (iDir == iDirNew) cGam2 -= xInv * w2;
        if (iInv == iInvNew) cGam3 += xDir * w2;
        if (iDir == iDirNew && iInv == iInvNew) cGam4 -= w2;
      }
    }

    // Solve (m2, Gamma) equation system => r2 * x-^2 + r1 * x- + r0 = 0.
    double cM0   = pow2(meanM) - cM1;
    double cGam0 = GammaNow - cGam1;
    double r2    = cM3 * cGam4 - cM4 * cGam3;
    double r1    = cM4 * cGam0 - cM0 * cGam4 + cM3 * cGam2 - cM2 * cGam3;
    double r0    = cM2 * cGam0 - cM0 * cGam2;
    double root  = sqrtpos( r1*r1 - 4. * r2 * r0 );
    if (abs(r2) < TINY || root < TINY) return Vec4(0., 0., 0., -1.);
    xInvHad      = 0.5 * (root / abs(r2) - r1 / r2);
    xDirHad      = (cM0 - cM3 * xInvHad) / (cM2 + cM4 * xInvHad);

    // Define position of new trial vertex.
    xDirNew = (iDirNew == iDirOld) ? xDirOld - xDirHad : 1. - xDirHad;
    xInvNew = (iInvNew == iInvOld) ? xInvOld + xInvHad : xInvHad;

    // Step up to new region if new x- > 1.
    if (xInvNew > 1.) {
      xInvHad = (iInvNew == iInvOld) ? 1. - xInvOld : 1.;
      xDirHad = 0.;
      pSoFarTmp += region.pHad( xPosHadTmp, xNegHadTmp, 0., 0.);
      --iInvNew;
      if (iInvNew < 0) return Vec4(0., 0., 0., -1.);
      continue;

    // Step down to new region if new x+ < 0.
    } else if (xDirNew < 0.) {
      xDirHad = (iDirNew == iDirOld) ? xDirOld : 1.;
      xInvHad = 0.;
      pSoFarTmp += region.pHad( xPosHadTmp, xNegHadTmp, 0., 0.);
      ++iDirNew;
      if (iDirNew + iInvNew > iMax) return Vec4(0., 0., 0., -1.);
      continue;
    }

    // Else we have found the correct region, and can return four-momentum.
    return pSoFarTmp + region.pHad( xPosHadTmp, xNegHadTmp, pxNewTmp,
      pyNewTmp);

  // End of "infinite" loop of stepping to new region.
  }

}

//--------------------------------------------------------------------------

// Update string end information after a hadron has been removed.

void StringEnd::update() {

  flavOld.anti(flavNew);
  iPosOld  = iPosNew;
  iNegOld  = iNegNew;
  pxOld    = -pxNew;
  pyOld    = -pyNew;
  GammaOld = GammaNew;
  xPosOld  = xPosNew;
  xNegOld  = xNegNew;
  colOld   = colNew;

}

//--------------------------------------------------------------------------

// Store old string end information.

void StringEnd::storePrev() {

  flavOldPrev  = flavOld;
  iPosOldPrev  = iPosOld;
  iNegOldPrev  = iNegOld;
  pxOldPrev    = pxOld;
  pyOldPrev    = pyOld;
  GammaOldPrev = GammaOld;
  xPosOldPrev  = xPosOld;
  xNegOldPrev  = xNegOld;
  colOldPrev   = colOld;

}

//--------------------------------------------------------------------------

// Update string end information to previous string break.

void StringEnd::updateToPrev() {

  flavOld  = flavOldPrev;
  iPosOld  = iPosOldPrev;
  iNegOld  = iNegOldPrev;
  pxOld    = pxOldPrev;
  pyOld    = pyOldPrev;
  GammaOld = GammaOldPrev;
  xPosOld  = xPosOldPrev;
  xNegOld  = xNegOldPrev;
  colOld   = colOldPrev;

}

//==========================================================================

// The StringFragmentation class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to (flavour-, energy) join the two string ends.
const int    StringFragmentation::NTRYFLAV      = 10;
const int    StringFragmentation::NTRYJOIN      = 30;

// The last few times gradually increase the stop mass to make it easier.
const int    StringFragmentation::NSTOPMASS     = 15;
const double StringFragmentation::FACSTOPMASS   = 1.05;

// For closed string, pick a Gamma by taking a step with fictitious mass.
const double StringFragmentation::CLOSEDM2MAX   = 25.;
const double StringFragmentation::CLOSEDM2FRAC  = 0.1;

// Do not allow too large argument to exp function.
const double StringFragmentation::EXPMAX        = 50.;

// Matching criterion that p+ and p- not the same (can happen in gg loop).
const double StringFragmentation::MATCHPOSNEG   = 1e-4;

// Fail and try again when two legs combined to diquark (3 loops).
const int    StringFragmentation::NTRYJNMATCH   = 20;
const double StringFragmentation::EEXTRAJNMATCH = 0.5;
const double StringFragmentation::MDIQUARKMIN   = -2.0;

// Protect against numerical precision giving zero or negative m2.
const double StringFragmentation::M2MINJRF      = 1e-4;
const double StringFragmentation::EMINJRF       = 1e-2;

// Iterate junction rest frame equation until convergence or too many tries.
const double StringFragmentation::CONVJRFEQ     = 1e-12;
const int    StringFragmentation::NTRYJRFEQ     = 40;

// Retry smearing of breakup vertex if too big shifts.
const int    StringFragmentation::NTRYSMEAR     = 100;

// Check that breakup vertex does not have negative tau^2 or t within roundoff.
const double StringFragmentation::CHECKPOS     = 1e-10;

//--------------------------------------------------------------------------

// Initialize and save pointers.

void StringFragmentation::init(StringFlav* flavSelPtrIn,
  StringPT* pTSelPtrIn, StringZ* zSelPtrIn, FragModPtr fragModPtrIn) {

  // Save pointers.
  flavSelPtr      = flavSelPtrIn;
  pTSelPtr        = pTSelPtrIn;
  zSelPtr         = zSelPtrIn;
  flavRopePtr     = fragModPtrIn;

  // Local copy of flavSel for possible dynamical modifications.
  flavSelNow = *flavSelPtr;

  // Initialize the StringFragmentation class.
  stopMass          = zSelPtr->stopMass();
  stopNewFlav       = zSelPtr->stopNewFlav();
  stopSmear         = zSelPtr->stopSmear();
  pNormJunction     = parm("StringFragmentation:pNormJunction");
  pMaxJunction      = 5 * pNormJunction;
  eBothLeftJunction = parm("StringFragmentation:eBothLeftJunction");
  eMaxLeftJunction  = parm("StringFragmentation:eMaxLeftJunction");
  eMinLeftJunction  = parm("StringFragmentation:eMinLeftJunction");

  // Calculation and definition of hadron space-time production vertices.
  hadronVertex    = mode("HadronVertex:mode");
  setVertices     = flag("Fragmentation:setVertices")
                 || flag("HadronLevel:Rescatter");
  kappaVtx        = parm("HadronVertex:kappa");
  smearOn         = flag("HadronVertex:smearOn");
  xySmear         = parm("HadronVertex:xySmear");
  maxSmear        = parm("HadronVertex:maxSmear");
  constantTau     = flag("HadronVertex:constantTau");
  maxTau          = parm("HadronVertex:maxTau");

  // Tracing of colours for primary hadrons.
  traceColours    = flag("StringFragmentation:TraceColours");

  // Joining of nearby partons along the string.
  mJoin           = parm("FragmentationSystems:mJoin");

  // Initialize the b parameter of the z spectrum, used when joining jets.
  bLund           = zSelPtr->bAreaLund();

  // Charm and bottom quark masses used for space-time offset.
  mc              = particleDataPtr->m0(4);
  mb              = particleDataPtr->m0(5);

  // MPI pT0, used for calculating effective number of strings.
  pT20            = pow2(parm("MultipartonInteractions:pT0Ref"));

  // Initialize the hadrons instance of an event record.
  hadrons.init( "(string fragmentation)", particleDataPtr);

  // Send on pointers to the two StringEnd instances.
  posEnd.init( particleDataPtr, flavSelPtr, pTSelPtr, zSelPtr, *settingsPtr);
  negEnd.init( particleDataPtr, flavSelPtr, pTSelPtr, zSelPtr, *settingsPtr);

  // Optionally allow for closepacking and reduced diquark formation.
  closePacking             = flag("ClosePacking:doClosePacking");
  closePackingTension      = parm("ClosePacking:tension");
  closePackingTensionRatio = parm("ClosePacking:tensionRatio");
  closePackingPT20         = pow2(parm("ClosePacking:PT0"));
  qqFacP                   = parm("ClosePacking:qqFacP");
  qqFacQ                   = parm("ClosePacking:qqFacQ");

  // Optionally hard baryon in beam remnant handling.
  dampPopcorn     = parm("BeamRemnants:dampPopcorn");
  hardRemn        = flag("BeamRemnants:hardRemnantBaryon");
  aRemn           = parm("BeamRemnants:aRemnantBaryon");
  bRemn           = parm("BeamRemnants:bRemnantBaryon");

  // Optionally allow pearl fragmentation.
  gluonPearl   = flag("StringFragmentation:pearlFragmentation");
  pearlFac     = parm("StringFragmentation:pearlProbFactor");

  // Optionally allow strangeness enhancement around the junction.
  strangeJunc  = flag("StringFragmentation:strangeJunctions");
  strangeParm  = parm("StringFragmentation:strangeJuncFactor");
  strangePearl = parm("StringFragmentation:strangePearlFactor");

}

//--------------------------------------------------------------------------

// Perform the fragmentation.

bool StringFragmentation::fragment( int iSub, const ColConfig& colConfig,
  Event& event) {

  // Find partons and their total four-momentum.
  iParton            = colConfig[iSub].iParton;
  iPos               = iParton[0];
  if (iPos < 0) iPos = iParton[1];
  int idPos          = event[iPos].id();
  iNeg               = iParton.back();
  int idNeg          = event[iNeg].id();
  pSum               = colConfig[iSub].pSum;

  // Set default values.
  kappaRatio = 1.;
  probQQmod = 1.;

  // Rapidity pairs [yCol, yAcol] of string piece ends.
  vector< vector< pair<double,double> > > rapPairs = colConfig.rapPairs;

  // Reset the local event record and vertex arrays.
  hadrons.clear();
  stringVertices.clear();
  legMinVertices.clear();
  legMidVertices.clear();
  posEnd.hadSoFar = 0;
  negEnd.hadSoFar = 0;
  nExtraJoin = 0;

  // Keep track of the string breaks of this string.
  StringBreaks stringBreaks;
  stringBreaks.setEnds(iPos, iNeg);

  // For closed gluon string: pick first breakup region.
  isClosed = colConfig[iSub].isClosed;
  if (isClosed) iParton = findFirstRegion(iSub, colConfig, event);

  // For junction topology: fragment off two of the string legs.
  // Then iParton overwritten to remaining leg + leftover diquark.
  pJunctionHadrons = 0.;
  hasJunction = colConfig[iSub].hasJunction;
  pearlFrag = false;
  eCutoff = 0.;
  if (hasJunction && !fragmentToJunction(event, rapPairs)) return false;
  int junctionHadrons = hadrons.size();
  bool usedPearl = (pearlFrag) ? false : true;
  if (!pearlFrag) pPearl = {0., 0., 0., 0. };
  if (hasJunction) {
    idPos  = event[ iParton[0] ].id();
    idNeg  = event[ iParton.back() ].id();
    pSum  -= pJunctionHadrons;
  }

  // Set up kinematics of string evolution ( = motion).
  system.setUp(iParton, event);
  stopMassNow = stopMass;
  bool extraStrange = false;
  bool allowPop = true;

  // Fallback loop, when joining in the middle fails.  Bailout if stuck.
  for ( int iTry = 0; ; ++iTry) {
    if (iTry > NTRYJOIN) {
      loggerPtr->ERROR_MSG("stuck in joining");
      if (nExtraJoin > 0) event.popBack(nExtraJoin);
      return false;
    }

    // After several failed tries join some (extra) nearby partons.
    if (iTry == NTRYJOIN / 3) {
      nExtraJoin += extraJoin( 2., event);
      system.setUp(iParton, event);
    }
    if (iTry == 2 * NTRYJOIN / 3) {
      nExtraJoin += extraJoin( 4., event);
      system.setUp(iParton, event);
    }

    // After several failed tries gradually allow larger stop mass.
    if (iTry > NTRYJOIN - NSTOPMASS) stopMassNow
      *= (max( abs(posEnd.flavOld.id), abs(negEnd.flavOld.id)) < 4)
      ? FACSTOPMASS : FACSTOPMASS * FACSTOPMASS;

    // Set up flavours of two string ends, and reset other info.
    setStartEnds(idPos, idNeg, system);
    pRem = pSum;

    // Begin fragmentation loop, interleaved from the two ends.
    bool fromPos;

    // Variables used to help identifying baryons from junction splittings.
    bool usedPosJun = false, usedNegJun = false;

    // Keep track of the momentum of hadrons taken from left and right.
    Vec4 hadMomPos, hadMomNeg;

    // Inform the UserHooks about the string to he hadronised.
    if ( userHooksPtr && userHooksPtr->canChangeFragPar() )
      userHooksPtr->setStringEnds(&posEnd, &negEnd, iParton);

    for ( ; ; ) {

      // Take a step either from the positive or the negative end.
      fromPos = (!usedPearl) ? false : (rndmPtr->flat() < 0.5);
      StringEnd& nowEnd = (fromPos) ? posEnd : negEnd;

      // Check how many nearby string pieces there are for the next hadron.
      kappaRatio = 1.;
      if (closePacking) kappaEffRatio(system, nowEnd, fromPos, iParton,
        rapPairs, pRem.mCalc(), event);

      // The FlavourRope treatment changes the fragmentation parameters.
      if (flavRopePtr) {
        if (!flavRopePtr->doChangeFragPar(flavSelPtr, zSelPtr, pTSelPtr,
          (fromPos ? hadMomPos.m2Calc() : hadMomNeg.m2Calc()), iParton,
          (fromPos ? idPos : idNeg)) )
          loggerPtr->ERROR_MSG(
            "FlavourRope failed to change fragmentation parameters");
      }

      // Possibility for a user to change the fragmentation parameters.
      if ( (userHooksPtr != 0) && userHooksPtr->canChangeFragPar() ) {
         if ( !userHooksPtr->doChangeFragPar( flavSelPtr, zSelPtr, pTSelPtr,
           (fromPos ? idPos : idNeg),
           (fromPos ? hadMomPos.m2Calc() : hadMomNeg.m2Calc()),
           iParton, &nowEnd) )
           loggerPtr->ERROR_MSG("failed to change hadronisation parameters");
      }

      // Check whether to use special hard diquark handling in beam remnant.
      bool forbidPopcornNow = false;
      if (dampPopcorn < 1. && !hasJunction && (nowEnd.hadSoFar == 0)) {
        int iNow = (fromPos) ? iPos : iNeg;
        if (event[iNow].isDiquark()) {
          int iMother = iNow;
          while (event[iMother].statusAbs() / 10 == 7)
            iMother = event[iMother].mother1();
          if (event[iMother].statusAbs() == 63
            && rndmPtr->flat() > dampPopcorn) forbidPopcornNow = true;
        }
      }

      // Optionally allow for strangeness enhancement around the junction.
      double strangeFac = 0.;
      if (strangeJunc) {
        if (!fromPos && hasJunction && !usedNegJun && !pearlFrag)
          strangeFac = strangeParm;
        else if (extraStrange) strangeFac = strangePearl;
      }

      if (!allowPop) nowEnd.flavOld.nPop = 0;
      // Construct trial hadron and check that energy remains.
      nowEnd.newHadron(
        kappaRatio, forbidPopcornNow, allowPop, strangeFac, probQQmod);
      if ( energyUsedUp(fromPos) ) {
        // Break out of fragmentation loop if no pearl or pearl is processed.
        if (usedPearl) break;
        // Unprocessed pearls require no diquarks next to the pearl.
        // Redo previous step if was a diquark or want strangeness enhancement.
        if ( int(hadrons.size()) > 0 &&
          (abs(nowEnd.flavOld.id) > 10 || (strangeJunc && !extraStrange))) {
          stringVertices.pop_back();
          hadMomNeg -= hadrons.back().p();
          pRem += hadrons.back().p();
          hadrons.popBack();
          nowEnd.hadSoFar -= 1;
          nowEnd.updateToPrev();
          if (strangeJunc) extraStrange = true;
          allowPop = false;
          continue;
        }
        // Redo current step if is a diquark or want strangeness enhancement.
        else if (abs(nowEnd.flavNew.id) > 10
          || (strangeJunc && !extraStrange)) {
          extraStrange = true;
          allowPop = false;
          continue;
        }
        // Make the pearl baryon and use in finalTwo joining.
        else {
          // Cannot process a system with diquark endpoint close to pearl.
          // Try combining the pearl with the other end if possible.
          if (abs(nowEnd.flavOld.id) > 10) {
            // Check if other string end is also a diquark.
            if ((fromPos && abs(negEnd.flavOld.id) > 10)
              || (!fromPos && abs(posEnd.flavOld.id) > 10)) {
              loggerPtr->ERROR_MSG(
                "cannot combine pearl and diquark endpoint");
              event.popBack(nExtraJoin);
              return false;
            // Swap ends to the non-diquark end if possible.
            } else {
              fromPos = !fromPos;
              nowEnd = (fromPos) ? posEnd : negEnd;
            }
          }
          nowEnd.pearlHadron(system, idPearl, pPearl);
          break;
        }
      }

      // Optionally allow a hard baryon fragmentation in beam remnant.
      bool useInputZ = false;
      double zUse    = 0.5;
      if (forbidPopcornNow && hardRemn) {
        useInputZ = true;
        zUse      = zSelPtr->zLund( aRemn, bRemn);
      }

      // Construct kinematics of the new hadron.
      Vec4 pHad = nowEnd.kinematicsHadron(system, newVertex,
        useInputZ, zUse);
      int statusHad = (fromPos) ? 83 : 84;
      if (!usedPearl) statusHad = 85;

      // Check if stepping over the pearl. If not, make the hadron.
      if (usedPearl || pRem.e() - pHad.e() > eCutoff) {
        nowEnd.hadSoFar += 1;
        extraStrange = false;
      }

      // If steppng over the pearl, make the junction baryon.
      else {

        // Unprocessed pearls require no diquarks next to the pearl.
        // Redo previous step if was a diquark or want strangeness enhancement.
        if ( int(hadrons.size()) > 0 &&
          (abs(nowEnd.flavOld.id) > 10 || (strangeJunc && !extraStrange))) {
          stringVertices.pop_back();
          hadMomNeg -= hadrons.back().p();
          pRem += hadrons.back().p();
          hadrons.popBack();
          nowEnd.hadSoFar -= 1;
          nowEnd.updateToPrev();
          if (strangeJunc) extraStrange = true;
          allowPop = false;
          continue;
        }

        // Redo current step if is a diquark or want strangeness enhancement.
        if (abs(nowEnd.flavNew.id) > 10 || strangeJunc) {
          if (strangeJunc) strangeFac = strangePearl;
          nowEnd.newHadron(kappaRatio, false, false, strangeFac, true);
          pHad = nowEnd.kinematicsHadron(system, newVertex, false);
        }
        // Check again if stepped over pearl. If not, make the hadron.
        if (pRem.e() - pHad.e() > eCutoff) nowEnd.hadSoFar += 1;
        // Construct hadron including pearl quark.
        // Once junction baryon is made, go to regular q-qbar fragmentation.
        else {
          // Cannot process a system with diquark endpoint close to pearl.
          if (int(hadrons.size()) < 1 && abs(nowEnd.flavOld.id) > 10) {
            loggerPtr->ERROR_MSG("cannot combine pearl and diquark endpoint");
            event.popBack(nExtraJoin);
            return false;
          }
          // Construct pearl baryon.
          nowEnd.pearlHadron(system, idPearl, pPearl);
          pHad = nowEnd.kinematicsHadron(system, newVertex,
            false, zUse, true, pPearl);
          statusHad = 88;
          usedPearl = true;
          allowPop = true;
        }
      }

      // Assign status code 87 for regular junction baryons.
      // Temporarily use 871 and 872 to distinguish pos and neg end.
      if (abs(nowEnd.idHad) > 1000 && abs(nowEnd.idHad) < 10000
          && statusHad != 88) {
        if (fromPos && !usedPosJun && event[iPos].statusAbs() == 74)  {
          statusHad = 871;
          usedPosJun = true;
        } else if (!fromPos && !usedNegJun
            && (event[iNeg].statusAbs() == 74 || hasJunction))  {
          statusHad = 872;
          usedNegJun = true;
        }
      }

      // Possibility for a user to veto the hadron production.
      if ( (userHooksPtr != 0) && userHooksPtr->canChangeFragPar() ) {
        // Provide full particle info for veto decision.
        if ( userHooksPtr->doVetoFragmentation( Particle( nowEnd.idHad,
          statusHad, iPos, iNeg, 0, 0, 0, 0, pHad, nowEnd.mHad), &nowEnd ) )
          continue;
      }

      // Bookkeeping of momentum taken away.
      if (fromPos) hadMomPos += pHad;
      else         hadMomNeg += pHad;

      // Append new vertex.
      stringVertices.push_back( newVertex);

      // Append produced hadron and string break.
      int colHadOld = nowEnd.colOld;
      int colHadNew = nowEnd.colNew;
      if ( !nowEnd.fromPos ) swap(colHadOld, colHadNew);
      if ( !usedPearl ) {
        swap(colHadOld, colHadNew);
        if ( particleDataPtr->chargeType( nowEnd.idHad ) != 0 )
          nowEnd.idHad = -nowEnd.idHad;
      }
      hadrons.append( nowEnd.idHad, statusHad, iPos, iNeg,
        0, 0, colHadOld, colHadNew, pHad, nowEnd.mHad);
      stringBreaks.add(nowEnd.flavNew.id);
      if (pHad.e() < 0.) break;

      // Update string end and remaining momentum.
      nowEnd.storePrev();
      nowEnd.update();
      pRem -= pHad;

    // End of fragmentation loop.
    }

    // Check how many nearby string pieces there are for the last hadron.
    StringEnd& nowEnd = (fromPos) ? posEnd : negEnd;
    kappaRatio = 1.;
    if (closePacking) kappaEffRatio( system, nowEnd, fromPos, iParton,
      rapPairs, pRem.mCalc(), event);

    // When done, join in the middle. If this works, then really done.
    if ( finalTwo(fromPos, event, usedPosJun, usedNegJun, usedPearl, pPearl) )
      break;

    // Else remove produced particles (except from first two junction legs)
    // and start all over.
    int newHadrons = hadrons.size() - junctionHadrons;
    hadrons.popBack(newHadrons);
    stringVertices.clear();
    posEnd.hadSoFar = 0;
    negEnd.hadSoFar = 0;
    if (pearlFrag) usedPearl = false;
  }

  // Junctions & extra joins: remove fictitious end/pearl, restore
  // original partons.
  if (nExtraJoin > 0) {
    event.popBack(nExtraJoin);
    iParton = colConfig[iSub].iParton;
  }

  // Store the hadrons in the normal event record, ordered from one end.
  store(event);

  // Store the string breaks in the event record.
  event.addStringBreaks(stringBreaks);

  // Store hadron production space-time vertices.
  bool saneVertices = (setVertices) ? setHadronVertices( event) : true;

  // Done.
  return saneVertices;

}

//--------------------------------------------------------------------------

// Find region where to put first string break for closed gluon loop.

vector<int> StringFragmentation::findFirstRegion(int iSub,
  const ColConfig& colConfig, const Event& event) const {

  // Partons and their total four-momentum.
  vector<int> iPartonIn = colConfig[iSub].iParton;

  // Evaluate mass-squared for all adjacent gluon pairs.
  vector<double> m2Pair;
  double m2Sum = 0.;
  int size = iPartonIn.size();
  for (int i = 0; i < size; ++i) {
    double m2Now = 0.5 * event[ iPartonIn[i] ].p()
      * event[ iPartonIn[(i + 1)%size] ].p();
    m2Pair.push_back(m2Now);
    m2Sum += m2Now;
  }

  // Pick breakup region with probability proportional to mass-squared.
  double m2Reg = m2Sum * rndmPtr->flat();
  int iReg = -1;
  do m2Reg -= m2Pair[++iReg];
  while (m2Reg > 0. && iReg < size - 1);

  // Create reordered parton list, with breakup string region duplicated.
  vector<int> iPartonOut;
  for (int i = 0; i < size + 2; ++i)
    if (size > 0) iPartonOut.push_back( iPartonIn[(i + iReg)%size] );

  // Done.
  return iPartonOut;

}

//--------------------------------------------------------------------------

// Set flavours and momentum position for initial string endpoints.

void StringFragmentation::setStartEnds(int idPosIn, int idNegIn,
const StringSystem& systemNow, int legNow) {

  // Variables characterizing string endpoints: defaults for open string.
  int    idPos       = idPosIn;
  int    idNeg       = idNegIn;
  double px          = 0.;
  double py          = 0.;
  double Gamma       = 0.;
  double xPosFromPos = 1.;
  double xNegFromPos = 0.;
  double xPosFromNeg = 0.;
  double xNegFromNeg = 1.;

  // For closed gluon loop need to pick an initial flavour.
  if (isClosed) {
    do {
      int idTry = flavSelPtr->pickLightQ();
      FlavContainer flavTry(idTry, 1);
      flavTry = flavSelPtr->pick( flavTry);
      flavTry = flavSelPtr->pick( flavTry);
      idPos   = flavTry.id;
      idNeg   = -idPos;
    } while (idPos == 0);

    // Also need pT and breakup vertex position in region.
    pair<double, double> pxy = pTSelPtr->pxy(idPos);
    px = pxy.first;
    py = pxy.second;
    double m2Region = systemNow.regionLowPos(0).w2;
    double m2Temp   = min( CLOSEDM2MAX, CLOSEDM2FRAC * m2Region);
    do {
      double zTemp = zSelPtr->zFrag( idPos, idNeg, m2Temp);
      xPosFromPos  = 1. - zTemp;
      xNegFromPos  = m2Temp / (zTemp * m2Region);
    } while (xNegFromPos > 1.);
    Gamma = xPosFromPos * xNegFromPos * m2Region;
    xPosFromNeg = xPosFromPos;
    xNegFromNeg = xNegFromPos;
  }

  // Initialize two string endpoints.
  posEnd.setUp( true, iPos, idPos, systemNow.iMax,  px,  py,
                Gamma, xPosFromPos, xNegFromPos,
                systemNow.regionLowPos(0).colPos);
  negEnd.setUp( false, iNeg, idNeg, systemNow.iMax, -px, -py,
                Gamma, xPosFromNeg, xNegFromNeg,
                systemNow.regionLowNeg(0).colPos);
  // Store breakup vertex information from the first and last points.
  if (setVertices) {
    if (legNow == legMin) legMinVertices.push_back(
      StringVertex( true, 0, systemNow.iMax, xPosFromPos, xNegFromPos) );
    else if (legNow == legMid) legMidVertices.push_back(
      StringVertex( true, 0, systemNow.iMax, xPosFromPos, xNegFromPos) );
    else {
      stringVertices.push_back(
       StringVertex (true, 0, systemNow.iMax, xPosFromPos, xNegFromPos) );
      stringVertices.push_back(
       StringVertex( false, systemNow.iMax, 0, xPosFromNeg, xNegFromNeg) );
    }
  }
  // For closed gluon loop can allow popcorn on one side but not both.
  if (isClosed) {
    flavSelPtr->assignPopQ(posEnd.flavOld);
    flavSelPtr->assignPopQ(negEnd.flavOld);
    if (rndmPtr->flat() < 0.5) posEnd.flavOld.nPop = 0;
    else                       negEnd.flavOld.nPop = 0;
    posEnd.flavOld.rank = 1;
    negEnd.flavOld.rank = 1;
  }

  // Done.

}

//--------------------------------------------------------------------------

// Check remaining energy-momentum whether it is OK to continue.

bool StringFragmentation::energyUsedUp(bool fromPos) {

  // If remaining negative energy then abort right away.
  if (pRem.e() < 0.) return true;

  // Calculate W2_minimum and done if remaining W2 is below it.
  double wMin = stopMassNow
    + particleDataPtr->constituentMass(posEnd.flavOld.id)
    + particleDataPtr->constituentMass(negEnd.flavOld.id);
  if (fromPos) wMin += stopNewFlav
    * particleDataPtr->constituentMass(posEnd.flavNew.id);
  else         wMin += stopNewFlav
    * particleDataPtr->constituentMass(negEnd.flavNew.id);
  wMin *= 1. + (2. * rndmPtr->flat() - 1.) * stopSmear;
  w2Rem = pRem.m2Calc();
  if (w2Rem < pow2(wMin)) return true;

  // Else still enough energy left to continue iteration.
  return false;

}

//--------------------------------------------------------------------------

// Store hadron production vertices in the event record.

bool StringFragmentation::setHadronVertices( Event& event) {

  // Order breakup points from one end to the other.
  bool saneVertices = true;
  int vertexSize = stringVertices.size();
  vector<StringVertex> orderedVertices;
  for (int i = 0; i < vertexSize; ++i) if (stringVertices[i].fromPos)
    orderedVertices.push_back( stringVertices[i] );
  for (int i = vertexSize - 1; i >= 0; --i) if (!stringVertices[i].fromPos)
    orderedVertices.push_back( stringVertices[i] );

  // Obtain space-time picture for breakup points.
  vector<Vec4> longitudinal;
  int finalVertexPos = 0;
  vector<int> iPartonIn = (hasJunction) ? iPartonMax : iParton;
  int id1 = event[ iPartonIn.front() ].idAbs();
  int id2 = (hasJunction) ? idDiquark : event[ iPartonIn.back() ].idAbs();
  int iHadJunc = (hasJunction)
    ?  legMinVertices.size() + legMidVertices.size() - 2 : 0;

   // Loop over ordered vertices to obtain longitudinal space-time locations.
  for (int i = 0; i < vertexSize; ++i) {
    int iPosIn    = orderedVertices[i].iRegPos;
    int iNegIn    = orderedVertices[i].iRegNeg;
    double xPosIn = orderedVertices[i].xRegPos;
    double xNegIn = orderedVertices[i].xRegNeg;

    // Special case for complicated solutions in finalTwo: calculated after.
    if (iPosIn == -1 && iNegIn == -1) {
      longitudinal.push_back( Vec4( 0., 0., 0., 0.) );
      finalVertexPos = i;

    // The normal cases.
    } else {
      StringRegion currentRegion = system.region( iPosIn, iNegIn);
      Vec4 gluonOffset = currentRegion.gluonOffset( iPartonIn, event, iPosIn,
        iNegIn) / kappaVtx;
      Vec4 noOffset = (xPosIn * currentRegion.pPos +
        xNegIn * currentRegion.pNeg) / kappaVtx;
      Vec4 pRegion = (currentRegion.pPos + currentRegion.pNeg) / kappaVtx;
      Vec4 fromBreaks = noOffset + gluonOffset;

      // Correction added to the space-time location of breakup points
      // if negative squared invariant time.
      if (fromBreaks.m2Calc() < 0.) {
        double cPlus = (-pRegion * fromBreaks + sqrt( pow2(pRegion
          * fromBreaks) - pRegion.m2Calc() * fromBreaks.m2Calc()))
          / pRegion.m2Calc();
        Vec4 pCorrection = noOffset + cPlus * pRegion;
        Vec4 fracCorrection;
        bool betterFrac = false;
        if (xPosIn < 0. || xPosIn > 1. || xNegIn < 0. || xNegIn > 1.) {
          xPosIn = max( 0., min( 1., xPosIn) );
          xNegIn = max( 0., min( 1., xNegIn) );
          fracCorrection = (xPosIn * currentRegion.pPos
            + xNegIn * currentRegion.pNeg) / kappaVtx;
          betterFrac = abs(noOffset.e() - pCorrection.e())
                     > abs(noOffset.e() - fracCorrection.e());
        }
        noOffset = (betterFrac) ? fracCorrection : pCorrection;
        fromBreaks = noOffset + gluonOffset;
      }

      // Store vertex and check positivity.
      longitudinal.push_back(fromBreaks);
      if (fromBreaks.m2Calc() < -CHECKPOS * max(1., pow2(fromBreaks.e())))
        loggerPtr->WARNING_MSG("negative tau^2 from breaks");
    }
  }

  // Breakup longitudinal space--time location for the special finalTwo case.
  if (finalVertexPos != 0) {
    double xPosIn = orderedVertices[finalVertexPos].xRegPos;
    double xNegIn = orderedVertices[finalVertexPos].xRegNeg;
    Vec4 v1 = longitudinal[finalVertexPos - 1];
    Vec4 v2 = longitudinal[finalVertexPos + 1];
    Vec4 vl = v1 - v2;
    Vec4 vk = pPosFinalReg + pNegFinalReg;
    double r = 0.;
    if (vl.m2Calc() > 0.) {
      Vec4 va = ( 0.5 * (v1 + v2) + xPosIn * pPosFinalReg +
        xNegIn * pNegFinalReg ) / kappaVtx;
      r = (va * vk - sqrt(pow2(va * vk) - vk.m2Calc() * va.m2Calc()) )
        / vk.m2Calc();
    } else r = 0.5 * sqrt(-vl.m2Calc() / vk.m2Calc());
    longitudinal[finalVertexPos] = ( 0.5 * (v1 + v2) + (xPosIn - r)
      * pPosFinalReg + (xNegIn - r) * pNegFinalReg ) / kappaVtx;
  }

  // Longitudinal offset of breakup points for massive quarks.
  for (int i = 0; i < vertexSize; ++i) {
    int iPosIn = orderedVertices[i].iRegPos;
    int iNegIn = orderedVertices[i].iRegNeg;
    if (iPosIn != -1) {
      StringRegion currentRegion = system.region( iPosIn, iNegIn);
      if ( currentRegion.massiveOffset( iPosIn, iNegIn, system.iMax,
        id1, id2, mc, mb) ) {
        Vec4 v1, v2;

        // Initial endpoint correction.
        if (i == 0 && (id1 == 4 || id1 == 5)) {
          int iPosIn2 = orderedVertices[i + 1].iRegPos;
          int iNegIn2 = orderedVertices[i + 1].iRegNeg;
          v2 = longitudinal[i + 1];
          double mHad =  event[event.size() + iHadJunc - hadrons.size()].m();
          double pPosMass = particleDataPtr->m0(id1);
          if (iPosIn == iPosIn2 && iNegIn == iNegIn2) {
            v1 = longitudinal[i];
            longitudinal[i] = v1 + (pPosMass / mHad) * (v2 - v1);
            if (longitudinal[i].m2Calc()
              < -CHECKPOS * max(1., pow2(longitudinal[i].e())))
              loggerPtr->WARNING_MSG(
                "negative tau^2 for endpoint massive correction");
          } else {
            StringRegion region2 = system.region( iPosIn2, iNegIn2);
            Vec4 gluonOffset = currentRegion.gluonOffset( iPartonIn, event,
              iPosIn, iNegIn);
            v1 = (region2.pPos + gluonOffset) / kappaVtx;
            longitudinal[i] = v1 + (pPosMass / mHad) * (v2 - v1);
            if (longitudinal[i].m2Calc()
              < -CHECKPOS * max(1., pow2(longitudinal[i].e())))
              loggerPtr->WARNING_MSG(
                "negative tau^2 for endpoint massive correction");
            continue;
          }
        }

        // Final endpoint correction.
        if (i == vertexSize - 1 && (id2 == 4 || id2 == 5) && !hasJunction) {
          int iPosIn2 = orderedVertices[i - 1].iRegPos;
          int iNegIn2 = orderedVertices[i - 1].iRegNeg;
          double mHad =  event[i - 1 + event.size() + iHadJunc
            - hadrons.size()].m();
          double pNegMass = particleDataPtr->m0(id2);
          if (iPosIn == iPosIn2 && iNegIn == iNegIn2) {
            v1 = longitudinal[i];
            v2 = longitudinal[i - 1] + currentRegion.massOffset / kappaVtx;
            longitudinal[i] = v1 + (pNegMass / mHad) * (v2 - v1);
            if (longitudinal[i].m2Calc()
              < -CHECKPOS * max(1., pow2(longitudinal[i].e())))
              loggerPtr->WARNING_MSG(
                "negative tau^2 for endpoint massive correction");

          } else {
            StringRegion region2 = system.region( iPosIn2, iNegIn2);
            Vec4 gluonOffset = currentRegion.gluonOffset( iPartonIn, event,
              iPosIn, iNegIn);
            v1 = (region2.pNeg + gluonOffset) / kappaVtx;
            v2 = longitudinal[i - 1];
            longitudinal[i] = v1 + (pNegMass / mHad) * (v2 - v1);
            if (longitudinal[i].m2Calc()
              < -CHECKPOS * max(1., pow2(longitudinal[i].e())))
              loggerPtr->WARNING_MSG(
                "negative tau^2 for endpoint massive correction");
            continue;
          }
        }

        // Add mass offset for all breakup points.
        Vec4 massOffset = currentRegion.massOffset / kappaVtx;
        Vec4 position = longitudinal[i] - massOffset;
        if (position.m2Calc() < 0.) {
          double cMinus = 0.;
          if (position.m2Calc() > -CHECKPOS * max(1., pow2(position.e())))
            position.e( position.pAbs() );
          else {
            if (massOffset.m2Calc() > 1e-6)
              cMinus = (longitudinal[i] * massOffset
                - sqrt(pow2(longitudinal[i] * massOffset)
                - longitudinal[i].m2Calc() * massOffset.m2Calc()))
                / massOffset.m2Calc();
            else cMinus = 0.5 * longitudinal[i].m2Calc()
                   / (longitudinal[i] * massOffset);
            position = longitudinal[i] - cMinus * massOffset;
          }
        }
        longitudinal[i] = position;
      }
    }
  }

  // Begin smearing in transverse space. Endpoint vertices unchanged.
  vector<Vec4> spaceTime;
  for (int i = 0; i < vertexSize; ++i) {
    Vec4& longi = longitudinal[i];
    Vec4 positionTot = longi;

    // Find two spacelike transverse four-vector directions.
    if (smearOn && (isClosed || (i > 0 && i < vertexSize - 1))) {
      Vec4 eX, eY;
      int iPosIn = orderedVertices[i].iRegPos;
      int iNegIn = orderedVertices[i].iRegNeg;
      if (iPosIn == -1 && iNegIn == -1) {
        eX = eXFinalReg;
        eY = eYFinalReg;
      } else{
        StringRegion currentRegion = system.region(iPosIn, iNegIn);
        eX = currentRegion.eX;
        eY = currentRegion.eY;
      }

      // Loop over tries; give up if struck.
      double longiLen = sqrt(longi.pAbs2() + pow2(longi.e()) + pow2(xySmear));
      for (int iTry = 0; ; ++iTry) {
        if (iTry == NTRYSMEAR) {
          loggerPtr->WARNING_MSG("failed to smear vertex (normal string)");
          positionTot = longi;
          break;
        }

        // Smearing calculated randomly following a Gaussian.
        double transX = rndmPtr -> gauss();
        double transY = rndmPtr -> gauss();
        Vec4 transPos = xySmear * (transX * eX + transY * eY) / sqrt(2.);
        positionTot = transPos + longi;

        // Keep proper or actual time constant when including the smearing.
        if (constantTau)
          positionTot.e( sqrt(longi.m2Calc() + positionTot.pAbs2()) );
        if ( sqrt(transPos.pAbs2() + pow2(positionTot.e() - longi.e()))
          < maxSmear * longiLen) break;
      }
    }
    spaceTime.push_back(positionTot);
  }

  // Space-time location of the breakup points in two initial junction legs.
  vector<Vec4> legMinSpaceTime, legMidSpaceTime;
  if (hasJunction && !pearlFrag) {
    int hadSoFar = 0;
    // Loop over the two lowest-energy legs.
    for (int legLoop = 0; legLoop < 2; ++legLoop) {
      vector<StringVertex> legVertices = (legLoop == 0) ? legMinVertices
        : legMidVertices;
      StringSystem& systemNow =  (legLoop == 0) ? systemMin : systemMid;
      vector<int> iPartonNow = (legLoop == 0) ? iPartonMinLeg : iPartonMidLeg;
      vector<Vec4> longitudinalPos;
      if (int(legVertices.size()) < 2) continue;

      // Loop over all breakup points of the specific leg.
      for (int i = 0; i < int(legVertices.size()); i++) {

        // Obtain longitudinal space-time location of breakup vertices.
        int iPosIn = legVertices[i].iRegPos;
        int iNegIn = legVertices[i].iRegNeg;
        double xPosIn = legVertices[i].xRegPos;
        double xNegIn = legVertices[i].xRegNeg;
        StringRegion currentRegion = systemNow.region( iPosIn, iNegIn);
        Vec4 gluonOffset = currentRegion.gluonOffsetJRF( iPartonNow, event,
          iPosIn, iNegIn, MtoJRF) / kappaVtx;
        Vec4 noOffset = (xPosIn * currentRegion.pPos
          + xNegIn * currentRegion.pNeg) / kappaVtx;
        Vec4 pRegion = (currentRegion.pPos + currentRegion.pNeg) / kappaVtx;
        Vec4 fromBreaks = noOffset + gluonOffset;

        // Correction added to the space-time location of breakup points
        // if negative squared invariant time.
        if (fromBreaks.m2Calc() < 0.) {
          double cPlus = (-pRegion * fromBreaks + sqrt( pow2(pRegion
            * fromBreaks) - pRegion.m2Calc() * fromBreaks.m2Calc()))
            / pRegion.m2Calc();
          Vec4 pCorrection = noOffset + cPlus * pRegion;
          Vec4 fracCorrection;
          bool betterFrac = false;
          if (xPosIn < 0. || xPosIn > 1. || xNegIn < 0. || xNegIn > 1.) {
            xPosIn = max( 0., min( 1., xPosIn) );
            xNegIn = max( 0., min( 1., xNegIn) );
            fracCorrection = (xPosIn * currentRegion.pPos
              + xNegIn * currentRegion.pNeg) / kappaVtx;
            betterFrac = abs(noOffset.e() - pCorrection.e())
                       > abs(noOffset.e() - fracCorrection.e());
          }
          noOffset = (betterFrac) ? fracCorrection : pCorrection;
          fromBreaks = noOffset + gluonOffset;
        }

        // Store vertex and check positivity.
        longitudinalPos.push_back(fromBreaks);
        if (fromBreaks.m2Calc() < -CHECKPOS * max(1., pow2(fromBreaks.e())))
          loggerPtr->WARNING_MSG("negative tau^2 from breaks");
      }

      // Longitudinal offset of breakup points for massive quarks.
      for (int i = 0; i < int(legVertices.size()); ++i) {
        int iPosIn = legVertices[i].iRegPos;
        int iNegIn = legVertices[i].iRegNeg;
        StringRegion currentRegion = systemNow.region( iPosIn, iNegIn);
        int id = event[ iPartonNow.front() ].idAbs();

        if ( currentRegion.massiveOffset( iPosIn, iNegIn, systemNow.iMax,
          id, 2, mc, mb) ) {
          Vec4 v1, v2;
          // Initial endpoint correction.
          if (i == 0 && (id == 4 || id == 5)) {
            int iPosIn2 = legVertices[i + 1].iRegPos;
            int iNegIn2 = legVertices[i + 1].iRegNeg;
            v2 = longitudinalPos[i + 1];
            double mHad =  event[hadSoFar + event.size() - hadrons.size()].m();
            double pPosMass = particleDataPtr->m0(id);
            if (iPosIn == iPosIn2 && iNegIn == iNegIn2) {
              v1 = longitudinalPos[i];
              longitudinalPos[i] = v1 + (pPosMass / mHad) * (v2 - v1);
              if (longitudinalPos[i].m2Calc()
                < -CHECKPOS * max(1., pow2(longitudinalPos[i].e())))
                loggerPtr->WARNING_MSG(
                  "negative tau^2 for endpoint massive correction");
            } else {
              StringRegion region2 = systemNow.region( iPosIn2, iNegIn2);
              Vec4 gluonOffset =  currentRegion.gluonOffsetJRF( iPartonNow,
              event, iPosIn2, iNegIn2, MtoJRF);
              v1 = (region2.pPos + gluonOffset) / kappaVtx;
              longitudinalPos[i] = v1 + (pPosMass / mHad) * (v2 - v1);
              if (longitudinalPos[i].m2Calc()
                < -CHECKPOS * max(1., pow2(longitudinalPos[i].e())))
                loggerPtr->WARNING_MSG(
                  "negative tau^2 for endpoint massive correction");
              continue;
            }
          }

          // Add mass offset for all breakup points.
          Vec4 massOffset = currentRegion.massOffset / kappaVtx;
          Vec4 position = longitudinalPos[i] - massOffset;

          // Correction for non-physical  situations.
          if (position.m2Calc() < 0.) {
            double cMinus = 0.;
            if (position.m2Calc() > -CHECKPOS * max(1., pow2(position.e())))
              position.e( position.pAbs() );
            else {
              if (massOffset.m2Calc() > 1e-6)
                cMinus = (longitudinalPos[i] * massOffset
                  - sqrt(pow2(longitudinalPos[i] * massOffset)
                  - longitudinalPos[i].m2Calc() * massOffset.m2Calc()))
                  / massOffset.m2Calc();
              else cMinus = 0.5 * longitudinalPos[i].m2Calc()
                     / (longitudinalPos[i] * massOffset);
              position = longitudinalPos[i] - cMinus * massOffset;
            }
          }
          longitudinalPos[i] = position;
        }
      }

      // Begin smearing in transverse space.
      for (int i = 0; i < int(legVertices.size()); ++i) {
        Vec4& longi = longitudinalPos[i];
        Vec4 positionTot = longi;

        // Find two spacelike transverse four-vector directions.
        if (smearOn && i > 0) {
          int iPosIn = legVertices[i].iRegPos;
          int iNegIn = legVertices[i].iRegNeg;
          StringRegion currentRegion = systemNow.region( iPosIn, iNegIn);
          Vec4 eX = currentRegion.eX;
          Vec4 eY = currentRegion.eY;

          // Loop over tries; give up if struck.
          double longiLen = sqrt(longi.pAbs2() + pow2(longi.e())
            + pow2(xySmear));
          for (int iTry = 0; ; ++iTry) {
            if (iTry == NTRYSMEAR) {
              loggerPtr->WARNING_MSG(
                "failed to smear vertex (junction string)");
              positionTot = longi;
              break;
            }

            // Smearing calculated randomly following a Gaussian.
            double transX = rndmPtr->gauss();
            double transY = rndmPtr->gauss();
            Vec4 transPos = xySmear * (transX * eX + transY * eY) / sqrt(2.);
            positionTot = transPos + longi;

            // Keep proper or actual time constant when including the smearing.
            if (constantTau)
              positionTot.e( sqrt(longi.m2Calc() + positionTot.pAbs2()) );
            if ( sqrt(transPos.pAbs2() + pow2(positionTot.e() - longi.e()))
              < maxSmear * longiLen) break;
          }
        }

        // Boost from the rest frame of the junction to the original frame.
        positionTot.rotbst(MfromJRF);

        // Recalculate time to compensate for numerical precision loss
        // in iterative calculation of MfromJRF. Store final result.
        if (positionTot.m2Calc() < 0.)
          positionTot.e( positionTot.pAbs() );
        if (legLoop == 0) legMinSpaceTime.push_back(positionTot);
        else legMidSpaceTime.push_back(positionTot);

      // End of loops over all breakup point of two lowest-energy legs.
      }
      hadSoFar = hadSoFar + legVertices.size() - 1;
    }
  }

  // Calculate hadron production points of the two initial junction legs.
  if (hasJunction) {
    int hadSoFar = 0;

    // Loop over the two lowest-energy legs.
    for (int legLoop = 0; legLoop < 2; legLoop++) {
      vector<Vec4>& finalLocation = (legLoop == 0) ? legMinSpaceTime
        : legMidSpaceTime;
      vector<int> iPartonNow = (legLoop == 0) ? iPartonMinLeg : iPartonMidLeg;
      int id =  event[iPartonNow.front()].idAbs();

      // Calculate hadron production points from breakup vertices
      // using one of the three definitions.
      for (int i = 0; i < int(finalLocation.size()) - 1; ++i) {
        Vec4 middlePoint =  0.5 * (finalLocation[i] + finalLocation[i + 1]);
        if (abs(middlePoint.mCalc()) > maxTau) saneVertices = false;
        int iHad = i + hadSoFar + event.size() - hadrons.size();
        Vec4 pHad = event[iHad].p();
        Vec4 prodPoints = Vec4( 0., 0., 0., 0.);

        // Reduced oscillation period if system contains massive quarks.
        double mHad = event[iHad].m();
        double redOsc = (i == 0 && (id == 4 || id == 5))
           ? 1. - pow2(particleDataPtr->m0(id) / mHad) : 1.;

        // Hadron production points calculation depending on definition.
        if (hadronVertex == 0) prodPoints = middlePoint;
        else if (hadronVertex == 1)
          prodPoints = middlePoint + 0.5 * redOsc * pHad / kappaVtx;
        else {
          prodPoints = middlePoint - 0.5 * redOsc * pHad / kappaVtx;
          if (prodPoints.m2Calc() < 0. || prodPoints.e() < 0.) {
            double midpProd = redOsc * middlePoint * pHad;
            double tau0fac = 2. * (midpProd - sqrt( pow2(midpProd)
              - middlePoint.m2Calc() * pow2(redOsc * mHad)))
              / pow2(redOsc * mHad);
            prodPoints = middlePoint - 0.5 * tau0fac * redOsc * pHad
              / kappaVtx;
          }
        }
        event[iHad].vProd( event[iHad].vProd() + prodPoints * FM2MM );
      }

      // End of the two legs loop. Number of hadrons with stored vertices.
      if (finalLocation.size() > 0) hadSoFar += finalLocation.size() - 1;
    }
  }

  // Normal string system or last leg: calculate hadron production points
  // from breakup vertices using one of the three definitions.
  for (int i = 0; i < int(spaceTime.size()) - 1; ++i) {
    Vec4 middlePoint = 0.5 * (spaceTime[i] + spaceTime[i + 1]);
    if (abs(middlePoint.mCalc()) > maxTau) saneVertices = false;
    int iHad = i + iHadJunc + event.size() - hadrons.size();
    Vec4 pHad = event[iHad].p();
    Vec4 prodPoints = Vec4( 0., 0., 0., 0.);

    // Reduced oscillation period if system contains massive quarks.
    double redOsc = 1.;
    double mHad = event[iHad].m();
    if (i == 0 && (id1 == 4 || id1 == 5))
      redOsc = 1. - pow2( particleDataPtr->m0(id1) / mHad);
    if (i == int(spaceTime.size()) - 1 && (id2 == 4 || id2 == 5))
      redOsc = 1. - pow2( particleDataPtr->m0(id2) / mHad);

    // Hadron production points calculation depending on definition.
    if (hadronVertex == 0) prodPoints = middlePoint;
    else if (hadronVertex == 1)
      prodPoints = middlePoint + 0.5 * redOsc * pHad / kappaVtx;
    else {
      prodPoints = middlePoint - 0.5 * redOsc * pHad / kappaVtx;
      if (prodPoints.m2Calc() < 0. || prodPoints.e() < 0.) {
        double midpProd = redOsc * middlePoint * pHad;
        double tau0fac = 2. * ( midpProd - sqrt( pow2(midpProd)
          - middlePoint.m2Calc() * pow2(redOsc * mHad))) / pow2(redOsc * mHad);
        prodPoints = middlePoint - 0.5 * tau0fac *  redOsc * pHad / kappaVtx;
      }
    }
    event[iHad].vProd( event[iHad].vProd() + prodPoints * FM2MM );
  }

  // Done.
  if (!saneVertices)
    loggerPtr->ERROR_MSG("too large |tau| so make new try");
  return saneVertices;

}

//--------------------------------------------------------------------------

// Produce the final two partons to complete the system.

bool StringFragmentation::finalTwo(bool fromPos, const Event& event,
  bool usedPosJun, bool usedNegJun, bool usedPearlIn, Vec4 pPearlIn) {

  if (usedPearlIn) pPearlIn = {0., 0., 0., 0.};

  // Check whether we went too far in p+-.
  if (pRem.e() < 0.  || w2Rem < 0. || (hadrons.size() > 0
    && hadrons.back().e() < 0.) ) return false;
  if ( posEnd.iPosOld > negEnd.iPosOld || negEnd.iNegOld > posEnd.iNegOld)
    return false;
  if ( posEnd.iPosOld == negEnd.iPosOld && posEnd.xPosOld < negEnd.xPosOld)
    return false;
  if ( posEnd.iNegOld == negEnd.iNegOld && posEnd.xNegOld > negEnd.xNegOld)
    return false;

  // Impossible to join two diquarks.
  FlavContainer flav1 = (fromPos) ? posEnd.flavNew.anti() : posEnd.flavOld;
  FlavContainer flav2 = (fromPos) ? negEnd.flavOld : negEnd.flavNew.anti();
  if (flav1.isDiquark() && flav2.isDiquark()) return false;

  // Construct preliminary hadron pT as if no region changes.
  double pHadPrelim[2] = { 0.0, 0.0 };
  if (fromPos) {
    pHadPrelim[0] = negEnd.pxOld-posEnd.pxNew;
    pHadPrelim[1] = negEnd.pyOld-posEnd.pyNew;
  } else {
    pHadPrelim[0] = posEnd.pxOld-negEnd.pxNew;
    pHadPrelim[1] = posEnd.pyOld-negEnd.pyNew;
  }
  double pThadPrelim = sqrt( pow2(pHadPrelim[0]) + pow2(pHadPrelim[1]) );

  // Construct the final hadron from the leftover flavours. Break if stuck.
  int idHad;
  for (int iTry = 0; iTry < NTRYFLAV; ++iTry) {
    idHad = flavSelPtr->getHadronID(
      flav1, flav2, pThadPrelim, kappaRatio, true);
    if (idHad != 0) break;
  }
  if (idHad == 0) return false;

  // Store the final particle and its new pT, and construct its mass.
  if (fromPos) {
    negEnd.idHad = idHad;
    negEnd.pxNew = -posEnd.pxNew;
    negEnd.pyNew = -posEnd.pyNew;
    negEnd.mHad  = flavSelPtr->getHadronMassWin(idHad);
  } else {
    posEnd.idHad = idHad;
    posEnd.pxNew = -negEnd.pxNew;
    posEnd.pyNew = -negEnd.pyNew;
    posEnd.mHad  = flavSelPtr->getHadronMassWin(idHad);
  }

  // String region in which to do the joining.
  StringRegion region = finalRegion();
  if (region.isEmpty) return false;

  // Project remaining momentum along longitudinal and transverse directions.
  region.project( pRem);
  double pxRem   = region.px() - posEnd.pxOld - negEnd.pxOld;
  double pyRem   = region.py() - posEnd.pyOld - negEnd.pyOld;
  double xPosRem = region.xPos();
  double xNegRem = region.xNeg();

  // Share extra pT kick evenly between final two hadrons.
  posEnd.pxOld += 0.5 * pxRem;
  posEnd.pyOld += 0.5 * pyRem;
  negEnd.pxOld += 0.5 * pxRem;
  negEnd.pyOld += 0.5 * pyRem;
  // Project the pearl momentum along longitudinal and transverse directions.
  double xPosPearl = 0;
  double xNegPearl = 0;
  if (!usedPearlIn) {
    region.project( pPearlIn);
    double pxPearl   = region.px();
    double pyPearl   = region.py();
    xPosPearl = region.xPos();
    xNegPearl = region.xNeg();
    // Give extra pT from the pearl to the pearl baryon.
    if (fromPos) {
      negEnd.pxOld += pxPearl;
      negEnd.pyOld += pyPearl;
    } else {
      posEnd.pxOld += pxPearl;
      posEnd.pyOld += pyPearl;
    }
  }

  // Construct new pT and mT of the final two particles.
  posEnd.pxHad  = posEnd.pxOld + posEnd.pxNew;
  posEnd.pyHad  = posEnd.pyOld + posEnd.pyNew;
  posEnd.mT2Had = pow2(posEnd.mHad) + pow2(posEnd.pxHad)
    + pow2(posEnd.pyHad);
  negEnd.pxHad  = negEnd.pxOld + negEnd.pxNew;
  negEnd.pyHad  = negEnd.pyOld + negEnd.pyNew;
  negEnd.mT2Had = pow2(negEnd.mHad) + pow2(negEnd.pxHad)
    + pow2(negEnd.pyHad);

  // Construct remaining system transverse mass, pxHad^2 + pyHad^2 + w2Rem,
  // including p+ and p- pearl contribution in w2Rem.
  Vec4 pRemPearl = pRem + pPearlIn;
  double wT2Rem = pRemPearl.m2Calc() + pow2( posEnd.pxHad + negEnd.pxHad)
    + pow2( posEnd.pyHad + negEnd.pyHad);

  // Check that kinematics possible.
  if ( sqrt(wT2Rem) < sqrt(posEnd.mT2Had) + sqrt(negEnd.mT2Had) )
    return false;
  double lambda2 = pow2( wT2Rem - posEnd.mT2Had - negEnd.mT2Had)
    - 4. * posEnd.mT2Had * negEnd.mT2Had;
  if (lambda2 <= 0.) return false;

  // Construct kinematics, as viewed in the transverse rest frame.
  double lambda = sqrt(lambda2);
  double probReverse = 1. / (1. + exp( min( EXPMAX, bLund * lambda) ) );
  double xpzPos = 0.5 * lambda/ wT2Rem;
  if (probReverse > rndmPtr->flat()) xpzPos = -xpzPos;
  double xmDiff = (posEnd.mT2Had - negEnd.mT2Had) / wT2Rem;
  double xePos  = 0.5 * (1. + xmDiff );
  double xeNeg  = 0.5 * (1. - xmDiff );

  // Translate this into kinematics in the string frame.
  Vec4 pHadPos = region.pHad( (xePos + xpzPos) * (xPosRem + xPosPearl),
    (xePos - xpzPos) * (xNegRem + xNegPearl), posEnd.pxHad, posEnd.pyHad);
  Vec4 pHadNeg = region.pHad( (xeNeg - xpzPos) * (xPosRem + xPosPearl),
    (xeNeg + xpzPos) * (xNegRem + xNegPearl), negEnd.pxHad, negEnd.pyHad);

  // Project pHadPos momentum fraction on the positive region
  // to obtain breakup vertices with respect to that region.
  if (setVertices) {
    if (usedPearlIn) {
      StringRegion posRegion = system.region( posEnd.iPosOld, posEnd.iNegOld);
      posRegion.project(pHadPos);
      double xFromPosPos = posEnd.xPosOld - posRegion.xPos();
      double xFromPosNeg = posEnd.xNegOld + posRegion.xNeg();

      // Same, but projecting pHadNeg fractions on the negative region.
      StringRegion negRegion = system.region( negEnd.iPosOld, negEnd.iNegOld);
      negRegion.project(pHadNeg);
      double xFromNegPos = negEnd.xPosOld + negRegion.xPos();
      double xFromNegNeg = negEnd.xNegOld - negRegion.xNeg();

      // Store energy-momentum coordinates for the final breakup vertex.
      // If projections give valid results, store them as breakup fractions.
      if (xFromPosPos > 0. && xFromPosPos < 1. && xFromPosNeg > 0.
        && xFromPosNeg < 1.) newVertex.store( fromPos, posEnd.iPosOld,
        posEnd.iNegOld, xFromPosPos, xFromPosNeg);
      else if (xFromNegPos > 0. && xFromNegPos < 1. && xFromNegNeg > 0.
        && xFromNegNeg < 1.) newVertex.store( fromPos, negEnd.iPosOld,
        negEnd.iNegOld, xFromNegPos, xFromNegNeg);

      // If above procedures do not work, calculate a new zHad and use
      // the kinematicsHadron method, first from the positive end.
      else {
        double gammaPosOld = posEnd.GammaOld;
        double gammaNegOld = negEnd.GammaOld;
        double zNewReg = 0.;
        if (posEnd.hadSoFar == 0) zNewReg = wT2Rem / (wT2Rem + gammaNegOld);
        else zNewReg = 0.5 * ( sqrt( pow2(wT2Rem + gammaNegOld - gammaPosOld)
          + 4. * wT2Rem * gammaPosOld) - (wT2Rem + gammaNegOld - gammaPosOld) )
          / gammaPosOld;
        double zHad = (xePos + xpzPos) * zNewReg;
        Vec4 proof = posEnd.kinematicsHadron(system, newVertex, true, zHad);

        // Try negative-end kinematicsHadron method if positive-end one failed.
        if (proof.e() < -1e-8) {
          if (negEnd.hadSoFar == 0) zNewReg = wT2Rem / (wT2Rem + gammaPosOld);
          else zNewReg = 0.5 * ( sqrt( pow2(wT2Rem + gammaPosOld - gammaNegOld)
              + 4. * wT2Rem * gammaNegOld)
            - (wT2Rem + gammaPosOld - gammaNegOld) ) / gammaNegOld;
          zHad = (xeNeg + xpzPos) * zNewReg;
          proof = negEnd.kinematicsHadron( system, newVertex, true, zHad);

          // As last resort use the final region created by the
          // finalTwo method.
          if (proof.e() < -1e-8) {
            pPosFinalReg = region.pPos;
            pNegFinalReg = region.pNeg;
            eXFinalReg = region.eX;
            eYFinalReg = region.eY;
            newVertex.store( true, -1, -1, 1. - (xePos + xpzPos) * xPosRem,
              (xePos - xpzPos) * xNegRem);
          }
        }
      }
    }
    // If dealing with a pearl, use the final region created by the
    // finalTwo method.
    else {
      pPosFinalReg = region.pPos;
      pNegFinalReg = region.pNeg;
      eXFinalReg = region.eX;
      eYFinalReg = region.eY;
      newVertex.store(
        true, -1, -1, 1. - (xePos + xpzPos) * (xPosRem + xPosPearl),
        (xePos - xpzPos) * (xNegRem + xNegPearl));
    }
  }

  // Update status codes for junction baryons.
  int statusHadPos = 83;
  int statusHadNeg = 84;

  // Assign code 87 for regular junction baryons.
  // Use 871 and 872 to distinguish between fromPos or fromNeg, reassign later.
  // Check status of iInvDir for possible popcorn meson shifting baryon number.
  if (fromPos) {
    if (!usedPearlIn) statusHadPos = 88;
    else if (abs(posEnd.idHad) > 1000 && abs(posEnd.idHad) < 10000
        && !usedPosJun && event[iPos].statusAbs() == 74) {
      statusHadPos = 871;
      usedPosJun = true;
    }
    if (abs(idHad) > 1000 && abs(idHad) < 10000) {
      if ((!usedPosJun && event[iPos].statusAbs() == 74)
          || (!usedNegJun && (event[iNeg].statusAbs() == 74 || hasJunction))) {
        statusHadNeg = 872;
      }
    }
  } else {
    if (!usedPearlIn) statusHadNeg = 88;
    else if (abs(negEnd.idHad) > 1000 && abs(negEnd.idHad) < 10000
      && !usedNegJun && (event[iNeg].statusAbs() == 74 || hasJunction)) {
        statusHadNeg = 872;
        usedNegJun = true;
    }
    if (abs(idHad) > 1000 && abs(idHad) < 10000) {
      if ((!usedPosJun && event[iPos].statusAbs() == 74)
          || (!usedNegJun && (event[iNeg].statusAbs() == 74 || hasJunction))) {
        statusHadPos = 871;
      }
    }
  }

  int colMid = (fromPos? negEnd.colOld: posEnd.colOld);
  // Possibility for a user to veto the hadron production.
  if ( (userHooksPtr != 0) && userHooksPtr->canChangeFragPar() ) {
    // Provide full particle info for veto decision.
    if ( userHooksPtr->doVetoFragmentation(
      Particle( posEnd.idHad, statusHadPos, posEnd.iEnd, negEnd.iEnd,
        0, 0, posEnd.colOld, colMid, pHadPos, posEnd.mHad),
      Particle( negEnd.idHad, statusHadNeg, posEnd.iEnd, negEnd.iEnd,
        0, 0, colMid, negEnd.colOld, pHadNeg, negEnd.mHad),
      &posEnd, &negEnd ) ) return false;
  }

  // Add produced particles to the event record.
  hadrons.append( posEnd.idHad, statusHadPos, posEnd.iEnd, negEnd.iEnd,
    0, 0, posEnd.colOld, colMid, pHadPos, posEnd.mHad);
  hadrons.append( negEnd.idHad, statusHadNeg, posEnd.iEnd, negEnd.iEnd,
    0, 0, colMid, negEnd.colOld, pHadNeg, negEnd.mHad);

  // Append the new vertex.
  stringVertices.push_back( newVertex);

  // It worked.
  return true;

}

//--------------------------------------------------------------------------

//  Construct a special joining region for the final two hadrons.

StringRegion StringFragmentation::finalRegion() {

  // Simple case when both string ends are in the same region.
  if (posEnd.iPosOld == negEnd.iPosOld && posEnd.iNegOld == negEnd.iNegOld)
    return system.region( posEnd.iPosOld, posEnd.iNegOld);

  // Start out with empty region. (Empty used for error returns.)
  StringRegion region;

  // Add up all remaining p+.
  Vec4 pPosJoin;
  int colPos = system.regionLowPos(posEnd.iPosOld).colPos;
  int colNeg = system.regionLowNeg(negEnd.iNegOld).colNeg;
  if ( posEnd.iPosOld == negEnd.iPosOld) {
    double xPosJoin = posEnd.xPosOld - negEnd.xPosOld;
    if (xPosJoin < 0.) return region;
    pPosJoin = system.regionLowPos(posEnd.iPosOld).pHad( xPosJoin, 0., 0., 0.);
  } else {
    for (int iPosNow = posEnd.iPosOld; iPosNow <= negEnd.iPosOld; ++iPosNow) {
      if (iPosNow == posEnd.iPosOld) pPosJoin
        += system.regionLowPos(iPosNow).pHad( posEnd.xPosOld, 0., 0., 0.);
      else if (iPosNow == negEnd.iPosOld) pPosJoin
        += system.regionLowPos(iPosNow).pHad( 1. - negEnd.xPosOld, 0., 0., 0.);
      else pPosJoin += system.regionLowPos(iPosNow).pHad( 1., 0., 0., 0.);
    }
  }

  // Add up all remaining p-.
  Vec4 pNegJoin;
  if ( negEnd.iNegOld == posEnd.iNegOld) {
    double xNegJoin = negEnd.xNegOld - posEnd.xNegOld;
    if (xNegJoin < 0.) return region;
    pNegJoin = system.regionLowNeg(negEnd.iNegOld).pHad( 0., xNegJoin, 0., 0.);
  } else {
    for (int iNegNow = negEnd.iNegOld; iNegNow <= posEnd.iNegOld; ++iNegNow) {
      if (iNegNow == negEnd.iNegOld) pNegJoin
        += system.regionLowNeg(iNegNow).pHad( 0., negEnd.xNegOld, 0., 0.);
      else if (iNegNow == posEnd.iNegOld) pNegJoin
        += system.regionLowNeg(iNegNow).pHad( 0., 1. - posEnd.xNegOld, 0., 0.);
      else pNegJoin += system.regionLowNeg(iNegNow).pHad( 0., 1., 0., 0.);
    }
  }

  // For a closed gluon loop pPosJoin == pNegJoin and the above does not work.
  // So reshuffle; "perfect" for g g systems, OK in general.
  Vec4 pTest = pPosJoin - pNegJoin;
  if ( abs(pTest.px()) + abs(pTest.py()) + abs(pTest.pz()) + abs(pTest.e())
    < MATCHPOSNEG * (pPosJoin.e() + pNegJoin.e()) ) {
    Vec4 delta
      = system.regionLowPos(posEnd.iPosOld + 1).pHad( 1., 0., 0., 0.)
      - system.regionLowNeg(negEnd.iNegOld + 1).pHad( 0., 1., 0., 0.);
    // If reshuffle did not help then pick random axis to break tie.
    // (Needed for low-mass q-g-qbar with q-qbar perfectly parallel.)
    if ( abs(delta.px()) + abs(delta.py()) + abs(delta.pz()) + abs(delta.e())
      < MATCHPOSNEG * (pPosJoin.e() + pNegJoin.e()) ) {
      double cthe = 2. * rndmPtr->flat() - 1.;
      double sthe = sqrtpos(1. - cthe * cthe);
      double phi  = 2. * M_PI * rndmPtr->flat();
      delta = 0.5 * min( pPosJoin.e(), pNegJoin.e())
        * Vec4( sthe * sin(phi), sthe * cos(phi), cthe, 0.);
      loggerPtr->WARNING_MSG("random axis needed to break tie");
    }
    pPosJoin -= delta;
    pNegJoin += delta;
  }

  // Construct a new region from remaining p+ and p-.
  region.setUp( pPosJoin, pNegJoin, colPos, colNeg);
  if (region.isEmpty) return region;

  // Project the existing pTold vectors onto the new directions.
  Vec4 pTposOld = system.region( posEnd.iPosOld, posEnd.iNegOld).pHad(
    0., 0.,  posEnd.pxOld, posEnd.pyOld);
  region.project( pTposOld);
  posEnd.pxOld = region.px();
  posEnd.pyOld = region.py();
  Vec4 pTnegOld = system.region( negEnd.iPosOld, negEnd.iNegOld).pHad(
      0., 0.,  negEnd.pxOld, negEnd.pyOld);
  region.project( pTnegOld);
  negEnd.pxOld = region.px();
  negEnd.pyOld = region.py();

  // Done.
  return region;

}

//--------------------------------------------------------------------------

// Store the hadrons in the normal event record, ordered from one end.

void StringFragmentation::store(Event& event) {

  // Starting position.
  int iFirst = event.size();

  // Remove colour indices to avoid confusion by default.
  if ( !traceColours )
    for (int i = 0; i < hadrons.size(); ++i) {
      hadrons[i].col(0);
      hadrons[i].acol(0);
    }

  // Copy straight over from first two junction legs.
  if (hasJunction) {
    for (int i = 0; i < hadrons.size(); ++i)
      if (hadrons[i].status() == 85 || hadrons[i].status() == 86)
        event.append( hadrons[i] );
  }

  // Loop downwards, copying all from the positive end.
  for (int i = 0; i < hadrons.size(); ++i) {
    if (hadrons[i].status() == 83) {
      event.append( hadrons[i] );
    } else if (hadrons[i].status() == 871) {
      hadrons[i].status(87);
      event.append( hadrons[i] );
    }
  }

  // Loop upwards, copying all from the negative end.
  for (int i = hadrons.size() - 1; i >= 0 ; --i) {
    if (hadrons[i].status() == 84) {
      event.append( hadrons[i] );
    } else if (hadrons[i].status() == 872) {
      hadrons[i].status(87);
      event.append( hadrons[i] );
    } else if (hadrons[i].status() == 88)
      event.append( hadrons[i] );
  }

  int iLast = event.size() - 1;

  // Set decay vertex when this is displaced.
  if (event[posEnd.iEnd].hasVertex()) {
    Vec4 vDec = event[posEnd.iEnd].vDec();
    for (int i = iFirst; i <= iLast; ++i) event[i].vProd( vDec );
  }

  // Set lifetime of hadrons.
  for (int i = iFirst; i <= iLast; ++i)
    event[i].tau( event[i].tau0() * rndmPtr->exp() );

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i)
    if (iParton[i] >= 0) {
      event[ iParton[i] ].statusNeg();
      event[ iParton[i] ].daughters(iFirst, iLast);
    }

}

//--------------------------------------------------------------------------

// Fragment off two of the string legs in to a junction.

bool StringFragmentation::fragmentToJunction(Event& event,
  vector< vector< pair<double,double> > >& rapPairs) {

  // Identify range of partons on the three legs.
  // (Each leg begins with an iParton[i] = -(10 + 10*junctionNumber + leg),
  // and partons then appear ordered from the junction outwards.)
  int legBeg[3] = {0, 0, 0};
  fill(begin(legEnd), end(legEnd), 0);
  int leg = -1;
  // PS (4/10/2011) Protect against invalid systems
  if (iParton[0] > 0) {
    loggerPtr->ERROR_MSG("iParton[0] not a valid junctionNumber");
    return false;
  }
  for (int i = 0; i < int(iParton.size()); ++i) {
    if (iParton[i] < 0) {
      if (leg == 2) {
        loggerPtr->ERROR_MSG("unprocessed multi-junction system");
        return false;
      }
      legBeg[++leg] = i + 1;
    }
    else legEnd[leg] = i;
  }

  // Initialise variables for JRF finding procedure.
  listJRF.clear();
  weightJRF.clear();
  pSumJRF = 0.;
  weightSum = 0.;
  vPearl = 0.;
  for (leg = 0; leg < 3; ++leg) {
    iLeg[leg] = legBeg[leg];
    int j = iParton[iLeg[leg]];
    pLeg[leg] = event[j].p();
    m2Leg[leg] = max( event[j].m2Calc(), 0.);
    idLeg[leg] = event[j].id();
    endpoint[leg] = (iLeg[leg] == legEnd[leg]);
  }

  // Do iterative JRF finding procedure, stepping through partons on each leg.
  lastJRF = false;
  do {

    // Check for non-standard junction configurations.
    if (perturbedJRF(event)) continue;
    if (collinearPair(event)) continue;

    // Find the 120 degree JRF for the current set of partons.
    Vec4 vToJRF = junctionRestFrame(pLeg[0], pLeg[1], pLeg[2], false);

    // Check if a valid velocity has been found. Typically errors arise from
    // soft small mass 4-vectors. These are usually endpoints that are not
    // given a constituent mass, or clustered gluons.
    bool juncCoM = false;
    if( isnan(vToJRF.e()) ) {
      vToJRF = pLeg[0] + pLeg[1] + pLeg[2];
      vToJRF /= vToJRF.e();
      loggerPtr->WARNING_MSG("JRF finding failure, use CoM instead");
      if (endpoint[0] || endpoint[1] || endpoint[2]) lastJRF = true;
      else juncCoM = true;
    }

    // Store information on JRF and update legs.
    listJRF.push_back( vToJRF );
    updateLegs(event, listJRF.back(), juncCoM );

  } while (!lastJRF);


  // Allow for pearl-on-a-string gluon-approximation for simple cases.
  if (gluonPearl) {
    if (vPearl > 0 && abs(event[ iParton[legEnd[ legPearl ]] ].id()) < 10) {
      if (pearlOnAString(event, legPearl)) {
        pearlFrag = true;
        idPearl = event[ iParton[legEnd[ legPearl ]] ].id();
      }
    }
  }

  // Set up pearl-on-a-string gluon-approximation.
  if (pearlFrag) {

    // Find the largest energy leg in frame of fragmentation.
    iPartonMax.clear();
    double eLeg[2] = {0., 0.};
    int j = (legPearl == 2) ? 0: legPearl + 1;
    int k = (j == 2) ? 0: j + 1;
    for (int i = legEnd[j]; i >= legBeg[j]; --i)
      eLeg[0] += event[iParton[i]].e();
    for (int i = legEnd[k]; i >= legBeg[k]; --i)
      eLeg[1] += event[iParton[i]].e();

    int legBig, legSmall;
    if (eLeg[0] < eLeg[1]) {
      eCutoff = eLeg[1];
      legBig = k;
      legSmall = j;
    } else {
      eCutoff = eLeg[0];
      legBig = j;
      legSmall = k;
    }

    // Store partons on legBig in reverse order in original frame.
    for (int i = legEnd[legBig]; i >= legBeg[legBig]; --i)
      iPartonMax.push_back( iParton[i] );

    // Create massless pseudo-gluon to represent the leg as pearl on string.
    int iBaryon  = event.append( 21, 78, 0, 0, 0, 0, 0, 0,
      gPearl, 0.);
    iPartonMax.push_back( iBaryon );

    // Store partons on legSmall from junction outwards.
    // Reverse colour and endpoint id.
    for (int i = legBeg[legSmall]; i <= legEnd[legSmall]; ++i) {
      int iNew = event.append( event[ iParton[i] ] );
      int aCol = event[iParton[i]].acol();
      event[iNew].acol( event[iParton[i]].col() );
      event[iNew].col( aCol );
      if (i == legEnd[legSmall]) event[iNew].id( -event[iParton[i]].id() );
      iPartonMax.push_back( iNew );
    }

    nExtraJoin = legEnd[legSmall] - legBeg[legSmall] + 2;

    // Assign colours to the junction pearl.
    if (event[ iParton[legEnd[legBig]] ].id() < 0) {
      event[iBaryon].col(event[ iParton[legBeg[legBig]] ].acol());
      event[iBaryon].acol(event[ iParton[legBeg[legSmall]] ].col());
    } else {
      event[iBaryon].col(event[ iParton[legBeg[legSmall]] ].acol());
      event[iBaryon].acol(event[ iParton[legBeg[legBig]] ].col());
    }

    // Modify parton list to qqbar system with pearl-on-string.
    iParton = iPartonMax;

    // Done.
    return true;

  }

  // Following procedure for standard junction fragmentation.
  // Find average junction velocity using exponential weighting.
  double weightSoFar = 0.;
  Vec4 vJunction( 0., 0., 0., 0.);
  for (int i = 0; i < int(weightJRF.size()); ++i) {
    vJunction += listJRF[i] * ( exp(- weightSoFar / weightSum)
      - exp(- (weightJRF[i] + weightSoFar) / weightSum ) );
    weightSoFar += weightJRF[i];
  }
  vJunction /= (1 - exp(- weightSoFar / weightSum ));

  // Sum leg four-momenta in original system.
  Vec4 pInLeg[3];
  for (leg = 0; leg < 3; ++leg) {
    pInLeg[leg] = 0.;
    for (int i = legBeg[leg]; i <= legEnd[leg]; ++i)
      pInLeg[leg] += event[ iParton[i] ].p();
  }

  // If failed to construct average JRF, use the CoM frame instead.
  if (isnan(vJunction.e())) {
    loggerPtr->WARNING_MSG("failed to find average JRF, use CoM instead");
    vJunction = pInLeg[0] + pInLeg[1] + pInLeg[2];
  }

  // Construct boost matrix from original system to JRF.
  MtoJRF.reset();
  MtoJRF.bstback(vJunction);

  // Opposite operation: boost from JRF to original system.
  MfromJRF = MtoJRF;
  MfromJRF.invert();

  // Sum leg four-momenta in original system and in JRF.
  Vec4 pInJRF[3];
  for (leg = 0; leg < 3; ++leg) {
    pInJRF[leg] = pInLeg[leg];
    pInJRF[leg].rotbst(MtoJRF);
  }

  // Pick the two legs with lowest momenta in JRF.
  legMin = (pInJRF[0].pAbs() < pInJRF[1].pAbs()) ? 0 : 1;
  int legMax = 1 - legMin;
  if (pInJRF[2].pAbs() < min(pInJRF[0].pAbs(), pInJRF[1].pAbs()) ) legMin = 2;
  else if (pInJRF[2].pAbs() > max(pInJRF[0].pAbs(), pInJRF[1].pAbs()) )
    legMax = 2;
  legMid = 3 - legMin - legMax;

  // Save info on which status codes belong with the three legs.
  int iJunction = (-iParton[0]) / 10 - 1;
  event.statusJunction( iJunction, legMin, 85);
  event.statusJunction( iJunction, legMid, 86);
  event.statusJunction( iJunction, legMax, 83);

  // Temporarily copy the partons on the low-energy legs, into the JRF,
  // in reverse order, so (anti)quark leg end first.
  vector<int> iPartonMin;
  iPartonMinLeg.clear();
  for (int i = legEnd[legMin]; i >= legBeg[legMin]; --i) {
    if (setVertices) iPartonMinLeg.push_back( iParton[i] );
    int iNew = event.append( event[ iParton[i] ] );
    event[iNew].rotbst(MtoJRF);
    iPartonMin.push_back( iNew );
  }
  vector<int> iPartonMid;
  iPartonMidLeg.clear();
  for (int i = legEnd[legMid]; i >= legBeg[legMid]; --i) {
    if (setVertices) iPartonMidLeg.push_back( iParton[i] );
    int iNew = event.append( event[ iParton[i] ] );
    event[iNew].rotbst(MtoJRF);
    iPartonMid.push_back( iNew );
  }

  // Define fictitious opposing partons in JRF and store as string ends.
  Vec4 pOppose = pInJRF[legMin];
  pOppose.flip3();
  bool posId = event[ iPartonMin[0] ].id() > 0 ? true : false;
  int idOppose = - event[ iPartonMin[0] ].id();
  int colOppose = posId ? 0 : event[iPartonMin.back()].acol();
  int acolOppose = posId ? event[iPartonMin.back()].acol() : 0;
  int iOppose = event.append( idOppose, 77, 0, 0, 0, 0, colOppose, acolOppose,
    pOppose, 0.);
  iPartonMin.push_back( iOppose);

  pOppose = pInJRF[legMid];
  pOppose.flip3();
  posId = event[ iPartonMid[0] ].id() > 0 ? true : false;
  idOppose = - event[ iPartonMid[0] ].id();
  colOppose = posId ? 0 : event[iPartonMid.back()].acol();
  acolOppose = posId ? event[iPartonMid.back()].acol() : 0;
  iOppose = event.append( idOppose, 77, 0, 0, 0, 0, colOppose, acolOppose,
    pOppose, 0.);
  iPartonMid.push_back( iOppose);

  // Set up kinematics of string evolution in low-energy temporary systems.
  systemMin.setUp(iPartonMin, event);
  systemMid.setUp(iPartonMid, event);

  // Outer fallback loop, when too little energy left for third leg.
  int idMin = 0;
  int idMid = 0;
  Vec4 pDiquark;
  for ( int iTryOuter = 0; ; ++iTryOuter) {

    // Middle fallback loop, when much unused energy in leg remnants.
    double eLeftMin = 0.;
    double eLeftMid = 0.;
    for ( int iTryMiddle = 0; ; ++iTryMiddle) {

      // Loop over the two lowest-energy legs.
      for (int legLoop = 0; legLoop < 2; ++ legLoop) {
        int legNow = (legLoop == 0) ? legMin : legMid;

        // Find momentum of hadrons created thus far.
        Vec4 pHadSoFar = {0., 0., 0., 0.};
        if (legNow == legMid)
          for (int i = 0; i < hadrons.size(); ++i)
            pHadSoFar += hadrons[i].p();

        // Read in properties specific to this leg.
        StringSystem& systemNow = (legLoop == 0) ? systemMin : systemMid;
        vector<int>& iPartonNow = (legLoop == 0) ? iPartonMin : iPartonMid;
        int idPos = event[ iPartonNow[0] ].id();
        idOppose = event[ iPartonNow.back() ].id();
        double eInJRF = pInJRF[legNow].e();
        int statusHad = (legLoop == 0) ? 85 : 86;

        // Inner fallback loop, when a diquark comes in to junction.
        double eUsed = 0.;
        vector<StringVertex> junctionVertices;
        for ( int iTryInner = 0; ; ++iTryInner) {

          if (iTryInner > 2 * NTRYJNMATCH) {
            if (abs(idPos) > 10)
              loggerPtr->ERROR_MSG("diquark too close to junction");
            else loggerPtr->ERROR_MSG("caught in junction flavour loop");
            event.popBack( iPartonMin.size() + iPartonMid.size() );
            return false;
          }

          bool needBaryon = (abs(idPos) > 10 && iTryInner > NTRYJNMATCH);
          double eExtra   = (iTryInner > NTRYJNMATCH) ? EEXTRAJNMATCH : 0.;

          // Set up two string ends, and begin fragmentation loop.
          setStartEnds(idPos, idOppose, systemNow, legNow);
          eUsed = 0.;
          int nHadrons = 0;
          bool noNegE = true;

          // Keep track of hadron momentum.
          Vec4 hadMom;
          Vec4 pHadLeg(0., 0., 0., 0.);

          // Keep track of strangeness enhancement around the junction.
          bool extraStrange = false;
          Vec4 pHadPrev;

          // Inform the UserHooks about the string to he hadronised.
          if ( userHooksPtr && userHooksPtr->canChangeFragPar() )
            userHooksPtr->setStringEnds(&posEnd, 0, iPartonNow);

          for ( ; ; ++nHadrons) {

            // The FlavourRope treatment changes the fragmentation parameters.
            if (flavRopePtr) {
              if (!flavRopePtr->doChangeFragPar(flavSelPtr, zSelPtr, pTSelPtr,
                hadMom.m2Calc(), (legLoop == 0 ? iPartonMin : iPartonMid ),
                idPos )) loggerPtr->ERROR_MSG(
                  "FlavourRope failed to change fragmentation parameters");
            }

            // Possibility for a user to change the fragmentation parameters.
            if ( (userHooksPtr != 0) && userHooksPtr->canChangeFragPar() ) {
              if ( !userHooksPtr->doChangeFragPar( flavSelPtr, zSelPtr,
                pTSelPtr, idPos, hadMom.m2Calc(), iPartonNow, &posEnd) )
                loggerPtr->ERROR_MSG(
                  "failed to change hadronisation parameters");
            }

            // Check how many nearby string pieces there are for the
            // next hadron.
            kappaRatio = 1.;
            if (closePacking) kappaEffRatio(systemNow, posEnd, true,
              iPartonNow, rapPairs, max(0.,eInJRF - eUsed), event);

            // Construct trial hadron from positive end.
            double strangeFac = (extraStrange) ? strangeParm : 0.;
            posEnd.newHadron(kappaRatio, false, true, strangeFac, probQQmod);
            Vec4 pHad = posEnd.kinematicsHadron(systemNow, newVertex);

            // Possibility for a user to veto the hadron production.
            if ( (userHooksPtr != 0) && userHooksPtr->canChangeFragPar() ) {
              // Provide full particle info for veto decision.
              if ( userHooksPtr->doVetoFragmentation( Particle( posEnd.idHad,
                statusHad, iPos, iNeg, 0, 0, 0, 0,
                pHad, posEnd.mHad), &posEnd ) ) {
                --nHadrons;
                continue;
              }
            }

            // Negative energy signals failure in construction.
            if (pHad.e() < 0. ) { noNegE = false; break; }

            // Break if passed system midpoint ( = junction) in energy.
            // Exceptions: small systems, and/or with diquark end.
            bool delayedBreak = false;
            Vec4 pCheck = pInJRF[legMin] + pInJRF[legMid]
              - (pHadSoFar + hadMom + pHad);

            if ((eUsed + pHad.e() + eExtra > eInJRF) ||
              (pCheck.mCalc() < 0)) {

              // Allow for no breaks on legMin.
              if (nHadrons == 0 && abs(idPos) < 10 && legNow == legMin) break;

              // Allow for strangeness enhancement around the junction.
              else if (nHadrons > 0 && strangeJunc && !extraStrange) {
                junctionVertices.pop_back();
                hadrons.popBack();
                extraStrange = true;
                nHadrons -= 2;
                hadMom -= pHadPrev;
                eUsed -= pHadPrev.e();
                posEnd.updateToPrev();
                continue;
              }
              else if (nHadrons > 0 || !needBaryon) break;
              delayedBreak = true;
            }

            // Else construct kinematics of the new hadron and store it.
            hadrons.append( posEnd.idHad, statusHad, iPos, iNeg,
              0, 0, posEnd.colOld, posEnd.colNew, pHad, posEnd.mHad);

            // Append the new vertex.
            junctionVertices.push_back( newVertex);

            // Update used up energy and momentum.
            hadMom += pHad;
            eUsed += pHad.e();

            // Store information about previous break then update.
            pHadPrev = pHad;
            posEnd.storePrev();
            posEnd.update();

            // Delayed break in small systems, and/or with diquark end.
            if (delayedBreak) {
              ++nHadrons;
              break;
            }
          }

          // Possible to produce zero hadrons if the endpoint is not a diquark.
          if (iTryInner > NTRYJNMATCH && !noNegE && nHadrons == 0 &&
            abs(idPos) < 10) break;

          // End of fragmentation loop. Inner loopback if ends on a diquark.
          if ( noNegE && abs(posEnd.flavOld.id) < 10 ) break;
          hadrons.popBack(nHadrons);
          junctionVertices.clear();
          if (legNow == legMin) legMinVertices.clear();
          else legMidVertices.clear();
        }

        // End of one-leg fragmentation. Store end quark and remnant energy.
        if (legNow == legMin) {
          idMin = posEnd.flavOld.id;
          eLeftMin = eInJRF - eUsed;
        } else {
          idMid = posEnd.flavOld.id;
          eLeftMid = eInJRF - eUsed;
        }

        // Store breakup vertices in two different vectors depending on leg.
        if (setVertices) {
          for (int i = 0; i < int(junctionVertices.size()); ++i) {
            if (legNow == legMin)
                 legMinVertices.push_back( junctionVertices[i]);
            else legMidVertices.push_back( junctionVertices[i]);
          }
        }
      }

      // End of both-leg fragmentation.
      // Middle loopback if too much energy left.
      double eTrial = eBothLeftJunction + rndmPtr->flat() * eMaxLeftJunction;
      if (iTryMiddle > NTRYJNMATCH
        || ( min( eLeftMin, eLeftMid) < eBothLeftJunction
        && max( eLeftMin, eLeftMid) < eTrial ) ) break;
      hadrons.clear();
      legMinVertices.clear();
      legMidVertices.clear();
    }

    // Boost hadrons away from the JRF to the original frame.
    for (int i = 0; i < hadrons.size(); ++i) {
      hadrons[i].rotbst(MfromJRF);

      // Recalculate energy to compensate for numerical precision loss
      // in iterative calculation of MfromJRF.
      hadrons[i].e( hadrons[i].eCalc() );
      pJunctionHadrons += hadrons[i].p();
    }

    // Outer loopback if combined diquark mass too negative
    // or too little energy left in third leg.
    pDiquark = pInLeg[legMin] + pInLeg[legMid] - pJunctionHadrons;
    double m2Left = m2( pInLeg[legMax], pDiquark);
    if (iTryOuter >  NTRYJNMATCH || (pDiquark.mCalc() > MDIQUARKMIN
      && m2Left > eMinLeftJunction * pInLeg[legMax].e()) ) break;
    hadrons.clear();
    legMinVertices.clear();
    legMidVertices.clear();
    pJunctionHadrons = 0.;
  }

  // Now found solution; no more loopback. Remove temporary parton copies.
  event.popBack( iPartonMin.size() + iPartonMid.size() );

  // Construct and store an effective diquark string end from the
  // two remnant quark ends, for temporary usage.
  idDiquark = flavSelPtr->makeDiquark( idMin, idMid);
  int iDiquark  = event.append( idDiquark, 78, 0, 0, 0, 0, 0, 0,
    pDiquark, pDiquark.mCalc());

  // Find the partons on the last leg, again in reverse order.
  iPartonMax.clear();
  for (int i = legEnd[legMax]; i >= legBeg[legMax]; --i)
    iPartonMax.push_back( iParton[i] );
  iPartonMax.push_back( iDiquark );
  nExtraJoin = 1;

  // Recluster gluons nearby to diquark end when taken too much energy.
  int iPsize       = iPartonMax.size();
  double m0Diquark = event[iDiquark].m0();
  while (iPsize > 2) {
    Vec4 pGluNear = event[ iPartonMax[iPsize - 2] ].p();
    if ( (pDiquark + 0.5 * pGluNear).mCalc() > m0Diquark + mJoin ) break;
    pDiquark += pGluNear;
    event[iDiquark].p( pDiquark );
    event[iDiquark].m( pDiquark.mCalc() );
    iPartonMax.pop_back();
    --iPsize;
    iPartonMax[iPsize - 1] = iDiquark;
  }
  if ( idDiquark > 0 )
    event[iDiquark].acol(event[iPartonMax[iPsize - 2]].col());
  else
    event[iDiquark].col(event[iPartonMax[iPsize - 2]].acol());

  // Modify parton list to remaining leg + remnant of the first two.
  iParton = iPartonMax;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Protect against collinear pairs of massless partons in JRF finding.
// Given the energy of the collinear pair, either construct a diquark velocity
// or step to next partons.

bool StringFragmentation::collinearPair(Event& event) {

  // Go to CoM frame.
  Vec4 pCoM[3];
  Vec4 pSumTemp = pLeg[0] + pLeg[1] + pLeg[2];
  for (int i = 0; i < 3; ++i) {
    pCoM[i] = pLeg[i];
    pCoM[i].bstback( pSumTemp );
  }

  // Loop over legs.
  for (int i = 0; i < 3; ++i) {
    int j = (i == 2) ? 0 : i + 1;
    int k = (j == 2) ? 0 : j + 1;

    // Check if partons on legs i and j are massless collinear.
    if ( m2Leg[i] < M2MINJRF && m2Leg[j] < M2MINJRF
      && abs(costheta( pCoM[i], pCoM[j]) - 1) < 1e-4) {

      // For small energy partons, step to the next parton if possible.
      // Stop iterative procedure if we have low energy endpoints.
      for (int l = 0; l < 2; ++l) {
        int leg = (l == 0) ? i : j;
        if (pCoM[leg].e() < EMINJRF) {
          if (!endpoint[leg]) nextParton(event, leg);
          else lastJRF = true;
        }
      }

      // If soft endpoints and no previous JRFs found, use the CoM frame.
      if (lastJRF && int(listJRF.size()) == 0) {
        loggerPtr->WARNING_MSG(
          "soft massless endpoint encountered in junction system.");
        listJRF.push_back(pSumTemp);
        weightJRF.push_back(pNormJunction);
        weightSum = pNormJunction;
      }

      // Exit procedure here given a soft collinear parton.
      if (pCoM[i].e() < EMINJRF || pCoM[j].e() < EMINJRF) return true;

      // Attempt to make a diquark-type velocity with the collinear
      // pair. Store total available energy in CoM to make u/d type
      // diquark of mass mDiq.
      double eSum = pCoM[i].e() + pCoM[j].e();
      double mDiq = particleDataPtr->m0(2101);

      // If sufficient energy to form u/d type diquark mass.
      if (eSum > mDiq) {

        // Find velocity given diquark mass and energy eSum in
        // direction of collinear pair.
        double pDiq = sqrt( eSum*eSum - mDiq*mDiq );
        Vec4 vTemp = (pDiq / eSum) * pCoM[i] / pCoM[i].pAbs();
        vTemp.e(1.);

        // Add two boosts to give final result.
        Vec4 vToJRF( 0., 0., 0., 1.);
        vToJRF.bst(vTemp);
        vToJRF.bst(pSumTemp);
        vToJRF /= vToJRF.e();

        // Store velocity and weight.
        listJRF.push_back(vToJRF);
        updateWeights( pDiq, vToJRF);

        // Update momenta to next partons on legs.
        if (endpoint[i] || endpoint[j] || endpoint[k]) lastJRF = true;
        else for (int l = 0; l < 3; ++l) nextParton(event, l);

        // Finish collinear pairs check.
        return true;

      // If not sufficient energy for a diquark, step to next partons.
      } else {

        // If soft endpoints, stop iterative procedure.
        if (endpoint[i] || endpoint[j]) {
          lastJRF = true;

          // If no previous JRFs found, use the CoM frame.
          if (int(listJRF.size()) == 0) {
            loggerPtr->WARNING_MSG("JRF finding failure, use CoM instead");
            listJRF.push_back(pSumTemp);
            weightJRF.push_back(pNormJunction);
            weightSum = pNormJunction;
          }

        // Else step to next partons.
        } else {
          nextParton(event, i);
          nextParton(event, j);
        }

        // Finish collinear pairs check.
        return true;
      }
    }
  }

  // If no collinear pairs, return false.
  return false;

}

//--------------------------------------------------------------------------

// Check for pearl-on-a-string cases or where multiple partons have the same
// rest frame. For pearl-on-a-string, use the average velocity of the motion.

bool StringFragmentation::perturbedJRF(Event& event) {

  // Loop through legs.
  for (int i = 0; i < 3; ++i) {

    // Only consider endpoint partons with mass.
    if ( m2Leg[i] >= M2MINJRF ) {

      // Store momentum of legs j and k in rest frame of parton i.
      int j = (i == 2) ? 0 : i + 1;
      int k = (j == 2) ? 0 : j + 1;
      Vec4 pRF[2];
      double pSize[2];
      pRF[0] = pLeg[j];
      pRF[0].bstback(pLeg[i]);
      pSize[0] = pRF[0].pAbs();
      pRF[1] = pLeg[k];
      pRF[1].bstback(pLeg[i]);
      pSize[1] = pRF[1].pAbs();

      // If multiple endpoints at rest, stop iterative procedure.
      if (endpoint[i] && ((endpoint[j] && pSize[0] < EMINJRF)
        || (endpoint[k] && pSize[1] <  EMINJRF))) {

        lastJRF = true;
        Vec4 vToRF = pLeg[i];
        vToRF /= vToRF.e();

        // If no previous JRFs found, use the rest frame as the JRF.
        // Often encountered with diquark legs, therefore could not combine
        // diquark and quark endpoint in clustering.
        if (listJRF.size() == 0) {
          loggerPtr->WARNING_MSG("small junction system, two short legs");
          weightJRF.push_back(pNormJunction);
          weightSum = pNormJunction;
          listJRF.push_back(vToRF);

        // Standard update of weights using rest of parton i.
        } else {
          if (pSize[0] < pSize[1]) updateWeights(pSize[0], vToRF);
          else updateWeights(pSize[1], vToRF);
          listJRF.push_back(vToRF);
        }
      }

      // Check for multiple partons with the same RF, i.e. soft massive gluons.
      // Step to next partons on soft gluon legs.
      if (pSize[0] < EMINJRF || pSize[1] < EMINJRF) {
        if (pSize[0] < EMINJRF && !endpoint[j]) nextParton(event, j);
        if (pSize[1] < EMINJRF && !endpoint[k]) nextParton(event, k);
        if (!endpoint[i]) nextParton(event, i);
        return true;
      }

      // Check if the angle between j and k is greater than 120 degrees.
      // i.e. the Mercedes frame does not exist. Instead find the
      // perturbed JRF, calculated using the Ariadne frame.
      if ( costheta(pRF[0], pRF[1]) < -0.5 ) {
        Vec4 pCopy[3];
        for (int l = 0; l < 3; ++l) pCopy[l] = pLeg[l];

        // Boost to CoM frame of j and k partons.
        Vec4 pSumJK = pCopy[j] + pCopy[k];
        pCopy[i].bstback(pSumJK);
        pCopy[j].bstback(pSumJK);
        pCopy[k].bstback(pSumJK);

        // Construct boost to the Ariadne frame from j-k CoM frame by removing
        // parallel component between parton i and j-k.
        Vec4 bstArdne = dot3(pCopy[i], pCopy[j]) * pCopy[j] / pCopy[j].pAbs2();
        bstArdne.e( pCopy[i].e() );
        pCopy[i].bstback(bstArdne);
        pCopy[j].bstback(bstArdne);
        pCopy[k].bstback(bstArdne);

        // Find the shortest leg out of j and k to find determine the
        // timeframe.
        int iMin = ( pCopy[j].pAbs() < pCopy[k].pAbs() ) ? j : k;
        double pMin = pCopy[iMin].pAbs();

        // Find effective time of frame. Store two times: one for within the
        // hadronisation timescale (pNormJunction), and the other after that.
        double pSmall[2];
        pSmall[0] = pMin; pSmall[1] = 0;
        if (pSumJRF + pMin > pNormJunction && pSumJRF < pNormJunction) {
          pSmall[0] = pNormJunction - pSumJRF;
          pSmall[1] = pMin - pSmall[0];
          if (pSumJRF + pMin > pMaxJunction) {
            pSmall[1] = pMaxJunction - pNormJunction;
            lastJRF = true;
          }
        } else if (pSumJRF + pMin > pMaxJunction) {
          pSmall[0] = pMaxJunction - pSumJRF;
          lastJRF = true;
        }

        // Check whether another iteration is necessary.
        if (endpoint[iMin]) lastJRF = true;

        // Approximate pearl motion as dx/dt = v0 * exp(-t * mFac)
        // exact form: dx/dt = 1 / sqrt(1 + m^2 / (p0 - 2 * kappa * x)^2)
        double mFac = 1.75 / sqrt( m2Leg[i] );
        double v0 = pCopy[i].pAbs() / pCopy[i].e();

        // Store unit vector in direction of parton i.
        Vec4 piUnit(0., 0., 0., 1.);
        if (pCopy[i].pAbs() > EMINJRF) piUnit = pCopy[i] / pCopy[i].pAbs();

        // Store change of velocity if a pearl forms from time zero and
        // remains a pearl for the duration of pMaxJunction.
        if (int(weightJRF.size()) == 0 && lastJRF) {
          vPearl = v0 - v0 * exp(- pSmall[0] * mFac);
          legPearl = i;
        }

        // Find the average velocity over time pSmall.
        for (int l = 0; l < 2; ++l) {
          double x0 = 0;
          double xNow = v0 * (1 - exp(- pSmall[l] * mFac)) / mFac;
          if (l == 1) {
            if (pSmall[1] < EMINJRF) break;
            x0 = v0 * (1 - exp(-pSmall[0] * mFac)) / mFac;
            xNow = v0 * (1 - exp(- (pSmall[0] + pSmall[1]) * mFac)) / mFac;
          }
          double vAvg = (xNow - x0) / pSmall[l];

          // Combine vAvg and other boosts for overall boost to perturbed JRF.
          Vec4 bst = piUnit * vAvg;
          bst.e( 1.);
          Vec4 vToJRF( 0., 0., 0., 1.);
          vToJRF.bst(bst);
          vToJRF.bst(bstArdne);
          vToJRF.bst(pSumJK);
          vToJRF /= vToJRF.e();
          listJRF.push_back(vToJRF);

          // Update pSumJRF and weights for JRFs.
          double wVal = pSmall[l] / sqrt( 1 - vToJRF.pAbs2() );
          weightJRF.push_back( wVal );
          if (pSumJRF < pNormJunction) weightSum += wVal;
          pSumJRF += pSmall[l];
        }

        // Do not need to update legs if last iteration.
        if (lastJRF) return true;

        // Update pearl quark to velocity at time pMin.
        double vFinal = v0 * exp(- pMin * mFac);
        double piNew = sqrt( m2Leg[i] ) * vFinal / sqrt(1 - vFinal*vFinal);
        if (pCopy[i].pAbs() > EMINJRF) {
          pCopy[i] *= piNew / pCopy[i].pAbs();
          pCopy[i].e( sqrt(m2Leg[i] + piNew*piNew) );
        }
        pCopy[i].bst(bstArdne);
        pCopy[i].bst(pSumJK);
        pLeg[i] = pCopy[i];

        // Update pMin to next parton on the leg.
        nextParton(event, iMin);

        // Reduce other leg by pMin.
        int l = (j == iMin) ? k : j;
        pCopy[l] = (pCopy[l].pAbs() - pMin) * pCopy[l] / pCopy[l].pAbs();
        pCopy[l].e( sqrt( m2Leg[l] + pCopy[l].pAbs2()) );
        pCopy[l].bst(bstArdne);
        pCopy[l].bst(pSumJK);
        pLeg[l] = pCopy[l];

        // If too little energy, step to next parton or stop iteration.
        if (pLeg[l].e() < EMINJRF) {
          if (!endpoint[l]) nextParton(event, l);
          else lastJRF = true;
        }
        return true;
      }
    }
  }

  // Done checking for perturbed JRFs.
  return false;

}

//--------------------------------------------------------------------------

// Update the weights for iterative procedure. Return the effective pSmall.
// This ensures pSumJRF does not exceed pMaxJunction and weightSum is not
// updated for pSumJRF > pNormJunction.

double StringFragmentation::updateWeights(double pSmall, Vec4 vJunIn) {

  // Multiply weighting of pSmall by gamma factor. Do not exceed pMaxJunction.
  if (pSumJRF + pSmall > pMaxJunction)
    weightJRF.push_back( pMaxJunction - pSumJRF );
  else weightJRF.push_back( pSmall );
  weightJRF.back() = weightJRF.back() / sqrt( 1 - vJunIn.pAbs2() );

  // Sum gamma-factor scaled weightings up to pSumJRF = pNormJunction for
  // use in exponential averaging.
  if (pSumJRF < pNormJunction && pSumJRF + pSmall > pNormJunction) {
    weightSum += ( pNormJunction - pSumJRF ) / sqrt( 1 - vJunIn.pAbs2() );
  } else if (pSumJRF < pNormJunction) weightSum += weightJRF.back();

  // Do not exceed pMaxJunction.
  double pSmallEff = pSmall;
  if (pSumJRF + pSmall > pMaxJunction) {
    lastJRF = true;
    pSmallEff = pMaxJunction - pSumJRF;
  }

  // Update pSumJRF.
  pSumJRF += pSmallEff;

  return pSmallEff;

}

//--------------------------------------------------------------------------

// Update momenta of each junction leg for iterative procedure.

int StringFragmentation::updateLegs(Event& event, Vec4 vJunIn, bool juncCoM) {

  // Boost each momenta and store the size of the momenta in the iterative JRF.
  // For endpoints, use factor of 2 to allow for oscillations on the junction.
  double pSize[3];
  for (int i = 0; i < 3; ++i) {
    pLeg[i].bstback(vJunIn);
    pSize[i] = (endpoint[i]) ? 2 * pLeg[i].pAbs() : pLeg[i].pAbs();
  }

  // Find smallest momenta.
  int iMin = (pSize[0] < pSize[1]) ? 0 : 1;
  if (pSize[2] < pSize[iMin] ) iMin = 2;

  // If shortest leg is a massless soft endpoint, stop iterative procedure.
  if (m2Leg[iMin] < M2MINJRF && pSize[iMin] < EMINJRF && endpoint[iMin])
    lastJRF = true;

  // Update the frame weights and find the shortest effective momentum.
  double pSmall = updateWeights(pSize[iMin], vJunIn);

  // Do not update leg momenta if ending iterative procedure..
  if (lastJRF) return 0;

  // If using failsafe of centre-of-mass frame, update legs to the next parton.
  if (juncCoM) {
    for (int i = 0; i < 3; ++i) nextParton(event, i);
    return 0;
  }

  // Update momenta of each junction leg if continuing iterative procedure.
  for (int i = 0; i < 3; ++i) {

    // Special treatment of shortest leg.
    if (i == iMin) {

      // If an endpoint, allow for oscillations.
      if (endpoint[iMin]) {
        pLeg[i].flip3();
        pLeg[i].bst(vJunIn);

      // Step to next parton if possible.
      } else nextParton(event, i);
    }

    // Reduce momenta of legs by pSmall. If endpoint oscillation starts,
    // momentum direction will flip naturally as pLeg[i].pAbs() - pSmall < 0.
    else {
      pLeg[i] = (pLeg[i].pAbs() - pSmall) * pLeg[i] / pLeg[i].pAbs();
      pLeg[i].e( sqrtpos( pLeg[i].pAbs2() + m2Leg[i] ) );

      // Protect against small energy partons.
      if ( pLeg[i].e() > EMINJRF ) pLeg[i].bst(vJunIn);
      else {

        // If not at endpoint, step to next parton on leg.
        if ( !endpoint[i] ) nextParton(event, i);

        // Else stop frame finding here.
        else { lastJRF = true; break; }
      }
    }
  }

  return 0;

}

//--------------------------------------------------------------------------

// Step to the next parton on the given junction leg.

void StringFragmentation::nextParton(Event& event, int leg) {

  iLeg[leg] = iLeg[leg] + 1;
  int i = iParton[ iLeg[leg] ];
  pLeg[leg] = event[i].p();
  m2Leg[leg] = max( event[i].m2Calc(), 0.);
  idLeg[leg] = event[i].id();
  endpoint[leg] = (iLeg[leg] == legEnd[leg]);

}

//--------------------------------------------------------------------------

// Decide whether to use the gluon-approximation for pearl-on-a-string cases.

bool StringFragmentation::pearlOnAString(Event& event, int iMin) {

  // Only allow gluon-approx if the endpoint has at least its constituent mass.
  double mCon = particleDataPtr->m0( event[ iParton[legEnd[iMin]] ].id() );
  if (event[ iParton[legEnd[iMin]] ].mCalc() < mCon ) return false;

  // Boost the pearl quark to the Ariadne frame.
  Vec4 pCopy[3];
  for (int i = 0; i < 3; ++i) pCopy[i] = event[ iParton[legEnd[i]] ].p();
  int i = iMin;
  int j = (i == 2) ? 0 : i + 1;
  int k = (j == 2) ? 0 : j + 1;

  // Boost to CoM frame of j and k partons.
  Vec4 pSumJK = pCopy[j] + pCopy[k];
  pCopy[i].bstback(pSumJK);
  pCopy[j].bstback(pSumJK);

  // Construct boost to the Ariadne frame from j-k CoM frame by removing
  // parallel component between parton i and j-k.
  Vec4 bstArdne = dot3(pCopy[i], pCopy[j]) * pCopy[j] / pCopy[j].pAbs2();
  bstArdne.e( pCopy[i].e() );
  pCopy[i].bstback(bstArdne);

  // Probabilistically choose whether to use gluon-approximation.
  // Change in pearl velocity, with v=1/2 giving maximal probability.
  double probPearl = 1 - 1 / (1 + pow(4 * vPearl, pearlFac));

  // Find the momentum to give to pearl gluon based on Ariadne frame momentum.
  if ( rndmPtr->flat() < probPearl ) {

    // Store massless 4-vector for gPearl, defined by Ariadne frame momentum.
    // Store excess energy in pPearl.
    double eArdne = pCopy[i].e();
    double pArdne = pCopy[i].pAbs();
    pCopy[i].e( pArdne );
    gPearl = pCopy[i];
    pPearl = {0., 0., 0., eArdne - pArdne };

    // Boost back to original frame.
    gPearl.bst( bstArdne );
    gPearl.bst( pSumJK );
    pPearl.bst( bstArdne );
    pPearl.bst( pSumJK );

    return true;
  }

  // Done.
  return false;

}

//--------------------------------------------------------------------------

// Find the boost matrix to the mercedes rest frame of a junction.

Vec4 StringFragmentation::junctionRestFrame(Vec4& p0, Vec4& p1,
  Vec4& p2, bool angleCheck) {

  // Calculate masses and other invariants.
  Vec4 pSumJun  = p0 + p1 + p2;
  Vec4 pList[3] = {p0, p1, p2};
  double sHat   = pSumJun.m2Calc();
  double pp[3][3];
  pp[0][0]      = p0.m2Calc();
  pp[1][1]      = p1.m2Calc();
  pp[2][2]      = p2.m2Calc();
  pp[0][1] = pp[1][0] = p0 * p1;
  pp[0][2] = pp[2][0] = p0 * p2;
  pp[1][2] = pp[2][1] = p1 * p2;

  // Check angle between partons in rest frame of each parton.
  // If more than 120 degrees, use the rest frame as the initial JRF.
  // Note: only call for CR methods, not for iterative JRF finding procedure.
  if (angleCheck) {
    for (int i = 0; i < 3; ++i) {

      // Only consider partons with mass.
      if ( pp[i][i] > M2MINJRF ) {

        // Store momentum of legs j and k in rest frame of parton i.
        int j = (i == 2) ? 0 : i + 1;
        int k = (j == 2) ? 0 : j + 1;
        Vec4 pRF[2];
        pRF[0] = pList[j];
        pRF[0].bstback(pList[i]);
        pRF[1] = pList[k];
        pRF[1].bstback(pList[i]);

        // Check if the angle between j and k is greater than 120 degrees.
        // If so, use the rest frame of i as the initial time JRF.
        if ( costheta(pRF[0], pRF[1]) < -0.5 ) {
          Vec4 vToJRF = pList[i] / pList[i].e();
          return vToJRF;
        }
      }
    }
  }

  // Requirement (eiMax)_j = pi*pj/mj < (eiMax)_k = pi*pk/mk, used below,
  // here rewritten as pi*pj * mk < pi*pk * mj and squared.
  double eMax01 = pow2(pp[0][1]) * pp[2][2];
  double eMax02 = pow2(pp[0][2]) * pp[1][1];
  double eMax12 = pow2(pp[1][2]) * pp[0][0];

  // Initially pick i to be the most massive parton. but allow other tries.
  int i = (pp[1][1] > pp[0][0]) ? 1 : 0;
  if (pp[2][2] > max(pp[0][0], pp[1][1])) i = 2;
  int j, k;
  double ei     = 0.;
  double ej     = 0.;
  double ek     = 0.;
  for (int iTry = 0; iTry < 3; ++iTry) {

    // Pick j to give minimal eiMax, and k the third vector.
    if (i == 0) j = (eMax02 < eMax01) ? 2 : 1;
    else if (i == 1) j = (eMax12 < eMax01) ? 2 : 0;
    else j = (eMax12 < eMax02) ? 1 : 0;
    k = 3 - i - j;

    // Alternative names according to i, j, k conventions.
    double m2i  = pp[i][i];
    double m2j  = pp[j][j];
    double m2k  = pp[k][k];
    double pipj = pp[i][j];
    double pipk = pp[i][k];
    double pjpk = pp[j][k];

    // Trivial to find new parton energies if all three partons are massless.
    if (m2i < M2MINJRF) {
      ei        = sqrt( 2. * pipk * pipj / (3. * pjpk) );
      ej        = sqrt( 2. * pjpk * pipj / (3. * pipk) );
      ek        = sqrt( 2. * pipk * pjpk / (3. * pipj) );

    // Else find three-momentum range for parton i and values at extremes.
    } else {
      // Minimum when i is at rest.
      double piMin = 0.;
      double eiMin = sqrt(m2i);
      double ejMin = pipj / eiMin;
      double ekMin = pipk / eiMin;
      double pjMin = sqrtpos( ejMin*ejMin - m2j );
      double pkMin = sqrtpos( ekMin*ekMin - m2k );
      double fMin  = ejMin * ekMin + 0.5 * pjMin * pkMin - pjpk;
      // Maximum estimated when j + k is at rest, alternatively j at rest.
      double eiMax = (pipj + pipk)
                   / sqrt( max( M2MINJRF, m2j + m2k + 2. * pjpk) );
      if (m2j > M2MINJRF) eiMax = min( eiMax, pipj / sqrt(m2j) );
      double piMax = sqrtpos( eiMax*eiMax - m2i );
      double temp  = eiMax*eiMax - 0.25 *piMax*piMax;
      double pjMax = (eiMax * sqrtpos( pipj*pipj - m2j * temp )
        - 0.5 * piMax * pipj) / temp;
      double pkMax = (eiMax * sqrtpos( pipk*pipk - m2k * temp )
        - 0.5 * piMax * pipk) / temp;
      double ejMax = sqrt(pjMax*pjMax + m2j);
      double ekMax = sqrt(pkMax*pkMax + m2k);
      double fMax  = ejMax * ekMax + 0.5 * pjMax * pkMax - pjpk;

      // If unexpected values at upper endpoint then pick another parton.
      if (fMax > 0.) {
        int iPrel = (i + 1)%3;
        if (pp[iPrel][iPrel] > M2MINJRF) {i = iPrel; continue;}
        ++iTry;
        iPrel = (i + 2)%3;
        if (iTry < 3 && pp[iPrel][iPrel] > M2MINJRF) {i = iPrel; continue;}
      }

      // Start binary + linear search to find solution inside range.
      int iterMin = 0;
      int iterMax = 0;
      double pi   = 0.5 * (piMin + piMax);
      for (int iter = 0; iter < NTRYJRFEQ; ++iter) {

        // Derive momentum of other two partons and distance to root.
        ei = sqrt(pi*pi + m2i);
        temp = ei*ei - 0.25 * pi*pi;
        double pj = (ei * sqrtpos( pipj*pipj - m2j * temp )
          - 0.5 * pi * pipj) / temp;
        double pk = (ei * sqrtpos( pipk*pipk - m2k * temp )
          - 0.5 * pi * pipk) / temp;
        ej = sqrt(pj*pj + m2j);
        ek = sqrt(pk*pk + m2k);
        double fNow = ej * ek + 0.5 * pj * pk - pjpk;

        // Replace lower or upper bound by new value.
        if (fNow > 0.) { ++iterMin; piMin = pi; fMin = fNow;}
        else {++iterMax; piMax = pi; fMax = fNow;}

        // Pick next i momentum to explore, hopefully closer to root.
        if (2 * iter < NTRYJRFEQ
          && (iterMin < 2 || iterMax < 2 || 4 * iter < NTRYJRFEQ))
          { pi = 0.5 * (piMin + piMax); continue;}
        if (fMin < 0. || fMax > 0. || abs(fNow) < CONVJRFEQ * sHat) break;
        pi = piMin + (piMax - piMin) * fMin / (fMin - fMax);
      }

    // If arrived here then either succeeded or exhausted possibilities.
    } break;
  }

  // Now we know the energies in the junction rest frame.
  double eNew[3] = { 0., 0., 0.};
  eNew[i] = ei;
  eNew[j] = ej;
  eNew[k] = ek;

  // Boost (copy of) partons to their rest frame.
  Vec4 p0cm = p0;
  Vec4 p1cm = p1;
  Vec4 p2cm = p2;
  p0cm.bstback(pSumJun);
  p1cm.bstback(pSumJun);
  p2cm.bstback(pSumJun);

  // Construct difference vectors and the boost to junction rest frame.
  Vec4 pDir01      = p0cm / p0cm.e() - p1cm / p1cm.e();
  Vec4 pDir02      = p0cm / p0cm.e() - p2cm / p2cm.e();
  double pDiff01   = pDir01.pAbs2();
  double pDiff02   = pDir02.pAbs2();
  double pDiff0102 = dot3(pDir01, pDir02);
  double eDiff01   = eNew[0] / p0cm.e() - eNew[1] / p1cm.e();
  double eDiff02   = eNew[0] / p0cm.e() - eNew[2] / p2cm.e();
  double denom     = pDiff01 * pDiff02 - pDiff0102*pDiff0102;
  double coef01    = (eDiff01 * pDiff02 - eDiff02 * pDiff0102) / denom;
  double coef02    = (eDiff02 * pDiff01 - eDiff01 * pDiff0102) / denom;
  Vec4 vTemp   = coef01 * pDir01 + coef02 * pDir02;
  vTemp.e( sqrt(1. + vTemp.pAbs2()) );

  // Add two boosts to give final result.
  Vec4 vToJRF( 0., 0., 0., 1.);
  vToJRF.bstback(vTemp);
  vToJRF.bst(pSumJun);
  vToJRF /= vToJRF.e();
  return vToJRF;

}

//--------------------------------------------------------------------------

// When string fragmentation has failed several times,
// try to join some more nearby partons.

int StringFragmentation::extraJoin(double facExtra, Event& event) {

  // Keep on looping while pairs found below joining threshold.
  int nJoin  = 0;
  int iPsize = iParton.size();
  while (iPsize > 2) {

    // Look for the pair of neighbour partons (along string) with
    // the smallest invariant mass (subtracting quark masses).
    int iJoinMin    = -1;
    double mJoinMin = 2. * facExtra * mJoin;
    for (int i = 0; i < iPsize - 1; ++i) {
      Particle& parton1 = event[ iParton[i] ];
      Particle& parton2 = event[ iParton[i + 1] ];
      Vec4 pSumNow;
      pSumNow += (parton2.isGluon()) ? 0.5 * parton1.p() : parton1.p();
      pSumNow += (parton2.isGluon()) ? 0.5 * parton2.p() : parton2.p();
      double mJoinNow = pSumNow.mCalc();
      if (!parton1.isGluon()) mJoinNow -= parton1.m0();
      if (!parton2.isGluon()) mJoinNow -= parton2.m0();
      if (mJoinNow < mJoinMin) { iJoinMin = i; mJoinMin = mJoinNow; }
    }

    // Decide whether to join, if not finished.
    if (iJoinMin < 0 || mJoinMin > facExtra * mJoin) return nJoin;
    ++nJoin;

    // Create new joined parton.
    int iJoin1  = iParton[iJoinMin];
    int iJoin2  = iParton[iJoinMin + 1];
    int idNew   = (event[iJoin1].isGluon()) ? event[iJoin2].id()
                                            : event[iJoin1].id();
    int colNew  = event[iJoin1].col();
    int acolNew = event[iJoin2].acol();
    if (colNew == acolNew) {
      colNew    = event[iJoin2].col();
      acolNew   = event[iJoin1].acol();
    }
    Vec4 pNew   = event[iJoin1].p() + event[iJoin2].p();

    // Append joined parton to event record and reduce parton list.
    int iNew = event.append( idNew, 73, min(iJoin1, iJoin2),
      max(iJoin1, iJoin2), 0, 0, colNew, acolNew, pNew, pNew.mCalc() );
    iParton[iJoinMin] = iNew;
    for (int i = iJoinMin + 1; i < iPsize - 1; ++i)
      iParton[i] = iParton[i + 1];
    iParton.pop_back();
    --iPsize;

  // Done.
  }
  return nJoin;
}

//--------------------------------------------------------------------------

// Implemented in close packing to calculate the change in the string tension
// due to surrounding strings. Returns the ratio of kappaEff to kappa0.
// Note: optional exponential scaling is not taken into account here and
// instead applied afterwards.

void StringFragmentation::kappaEffRatio(StringSystem& systemNow,
  StringEnd end, bool fromPos, vector<int> partonList,
  vector< vector< pair<double,double> > >& rapPairs,
  double mRem, Event& event) {

  // Get temporary hadron momentum.
  // Elaborate check for thermal model.
  Vec4   pHad     = Vec4(0., 0., 0., -1.);
  double phi      = 2.0 * M_PI * rndmPtr->flat();
  double mult     = -1.0;
  int    nTryMax  = 100;
  double multStep = 5.0 / ((double)nTryMax/2);
  double multNow  = 1.0 + multStep;
  for (int i = 1; i <= nTryMax; i++) {
    pHad = end.kinematicsHadronTmp(systemNow, mRem, phi, mult);
    // If valid momentum found, done.
    if (pHad.e() > 0.0) break;
    // Set mult as multiplicative factor. Alternate between adding and
    // subtracting multStep.
    mult = 1.0;
    if (i%2 == 0) {
      mult    *= multNow;
      multNow += multStep;
    } else mult /= multNow;
  }

  // In case of failure, use remnant momentum.
  if (pHad.e() < 0.0) pHad = pRem;

  // Extract pT2 of hadron.
  double pT2Had     = pHad.pT2();
  // Now loop through the list of rapidity pairs and count strings
  // sitting at the hadron rapidity.
  Particle hadron = Particle();
  hadron.p(pHad); hadron.m(pHad.mCalc());
  double yHad = hadron.y();

  // p and q have same and opposite flux as fragmenting string respectively.
  // Start p at -1 to subtract off current string.
  int p = -1;
  int q = 0;

  // Find rapidity values defining current region.
  int iP = (fromPos) ? posEnd.iPosOld : posEnd.iNegOld;
  int iN = (fromPos) ? posEnd.iPosOld + 1 : posEnd.iNegOld + 1;
  double yPos = yMax(event[ partonList[iP] ], 0.1);
  double yNeg = yMax(event[ partonList[iN] ], 0.1);
  if (abs(yHad) > 1000) yHad = (fromPos) ? yPos : yNeg;

  // Find flux direction of current string piece.
  bool fluxDirPos = (event[posEnd.iEnd].colType() == 1) ? true : false;
  if (event[posEnd.iEnd].colType() == 2)
    if (event[iP].col() == event[iN].acol() && event[iP].col() != 0)
      fluxDirPos = true;

  // Redefine flux direction based on rapidity ordering.
  if (yNeg < yPos) fluxDirPos = !fluxDirPos;

  // Loop through systems of rapidity pairs.
  for (int iSub = 0; iSub < int(rapPairs.size()); iSub++) {
    vector< pair<double,double> > pairNow = rapPairs[iSub];

    // Loop over rapidity pairs and find overlaps.
    for (int iPair = 0; iPair < int(pairNow.size()); iPair++) {

      // Order y1 and y2 in same flux direction as current string.
      // pairNow[iPair].first/second has yCol/yAcol.
      double y1 = pairNow[iPair].first;
      double y2 = pairNow[iPair].second;
      if (!fluxDirPos) swap( y1, y2);

      // Allow for partial overlaps near endpoints. For continuous string,
      // each piece linearly takes over from the last, adding to 1.
      if ( (y1 < y2) && (y1 - 0.5 < yHad) && (yHad < y2 + 0.5) ) {
        // Linear tapering off near ends. At endpoint, overlap = 0.5.
        double ov1 = min( yHad - y1 + 0.5, 0.5);
        double ov2 = min( y2 - yHad + 0.5, 0.5);
        p += ov1 + ov2;
      } else if ( (y2 < y1) && (y2 - 0.5 < yHad) && (yHad < y1 + 0.5) ) {
        double ov1 = min( y1 - yHad + 0.5, 0.5);
        double ov2 = min( yHad - y2 + 0.5, 0.5);
        q += ov1 + ov2;
      }
    }
  }

  if (p == -1) p = 0;

  // Number of strings effective contribution taking pT into account.
  double kappaRatioEff = closePackingTension * (p + closePackingTensionRatio *
    q) / (1.0 + pT2Had / closePackingPT20);

  // Add back current string.
  kappaRatio = kappaRatioEff + 1.0;

  // pT suppress the probability of fluctuation reconnections.
  double qqFacPmod = qqFacP / (1.0 + pT2Had / closePackingPT20);
  double qqFacQmod = qqFacQ / (1.0 + pT2Had / closePackingPT20);

  // Survival probability for a diquark created via colour fluctuations.
  probQQmod = pow( 1 - qqFacPmod, p/9) * pow( 1 - qqFacQmod, q/9);

}

//==========================================================================

} // end namespace Pythia8
