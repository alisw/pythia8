// BeamSetup.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the BeamSetup class.

#include "Pythia8/BeamSetup.h"
#include "Pythia8/Plugins.h"

namespace Pythia8 {

//==========================================================================

// The BeamSetup class.

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

bool BeamSetup::setPDFPtr( PDFPtr pdfAPtrIn, PDFPtr pdfBPtrIn,
  PDFPtr pdfHardAPtrIn, PDFPtr pdfHardBPtrIn, PDFPtr pdfPomAPtrIn,
  PDFPtr pdfPomBPtrIn, PDFPtr pdfGamAPtrIn, PDFPtr pdfGamBPtrIn,
  PDFPtr pdfHardGamAPtrIn, PDFPtr pdfHardGamBPtrIn, PDFPtr pdfUnresAPtrIn,
  PDFPtr pdfUnresBPtrIn, PDFPtr pdfUnresGamAPtrIn, PDFPtr pdfUnresGamBPtrIn,
  PDFPtr pdfVMDAPtrIn, PDFPtr pdfVMDBPtrIn) {

  // Default is no pointer to different PDF kinds.
  pdfAPtr = pdfBPtr = pdfHardAPtr = pdfHardBPtr = pdfPomAPtr = pdfPomBPtr
    = pdfGamAPtr = pdfGamBPtr = pdfHardGamAPtr = pdfHardGamBPtr = pdfUnresAPtr
    = pdfUnresBPtr = pdfUnresGamAPtr = pdfUnresGamBPtr = pdfVMDAPtr
    = pdfVMDBPtr = nullptr;

  // Switch off external PDF's by zero as input.
  if ( !pdfAPtrIn && !pdfBPtrIn) return true;

  // The two PDF objects cannot be one and the same.
  if (pdfAPtrIn == pdfBPtrIn) return false;

  // Save pointers.
  pdfAPtr       = pdfAPtrIn;
  pdfBPtr       = pdfBPtrIn;

  // By default same pointers for hard-process PDF's.
  pdfHardAPtr   = pdfAPtrIn;
  pdfHardBPtr   = pdfBPtrIn;

  // Optionally allow separate pointers for hard process.
  if (pdfHardAPtrIn && pdfHardBPtrIn) {
    if (pdfHardAPtrIn == pdfHardBPtrIn) return false;
    pdfHardAPtr = pdfHardAPtrIn;
    pdfHardBPtr = pdfHardBPtrIn;
  }

  // Optionally allow pointers for Pomerons in the proton.
  if (pdfPomAPtrIn && pdfPomBPtrIn) {
    if (pdfPomAPtrIn == pdfPomBPtrIn) return false;
    pdfPomAPtr  = pdfPomAPtrIn;
    pdfPomBPtr  = pdfPomBPtrIn;
  }

  // Optionally allow pointers for Gammas in the leptons.
  if (pdfGamAPtrIn && pdfGamBPtrIn) {
    if (pdfGamAPtrIn == pdfGamBPtrIn) return false;
    pdfGamAPtr  = pdfGamAPtrIn;
    pdfGamBPtr  = pdfGamBPtrIn;
  }

  // Optionally allow pointers for Hard PDFs for photons in the leptons.
  if (pdfHardGamAPtrIn && pdfHardGamBPtrIn) {
    if (pdfHardGamAPtrIn == pdfHardGamBPtrIn) return false;
    pdfHardGamAPtr  = pdfHardGamAPtrIn;
    pdfHardGamBPtr  = pdfHardGamBPtrIn;
  }

  // Optionally allow pointers for unresolved PDFs.
  if (pdfUnresAPtrIn && pdfUnresBPtrIn) {
    if (pdfUnresAPtrIn == pdfUnresBPtrIn) return false;
    pdfUnresAPtr = pdfUnresAPtrIn;
    pdfUnresBPtr = pdfUnresBPtrIn;
  }

  // Optionally allow pointers for unresolved PDFs for photons from leptons.
  if (pdfUnresGamAPtrIn && pdfUnresGamBPtrIn) {
    if (pdfUnresGamAPtrIn == pdfUnresGamBPtrIn) return false;
    pdfUnresGamAPtr = pdfUnresGamAPtrIn;
    pdfUnresGamBPtr = pdfUnresGamBPtrIn;
  }

  // Optionally allow pointers for VMD in the gamma.
  if (pdfVMDAPtrIn && pdfVMDBPtrIn) {
    if (pdfVMDAPtrIn == pdfVMDBPtrIn) return false;
    pdfVMDAPtr  = pdfVMDAPtrIn;
    pdfVMDBPtr  = pdfVMDBPtrIn;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

bool BeamSetup::setPDFAPtr( PDFPtr pdfAPtrIn ) {

  // Reset pointers to be empty.
  pdfAPtr = pdfBPtr = pdfHardAPtr = pdfHardBPtr = pdfPomAPtr = pdfPomBPtr
    = pdfGamAPtr = pdfGamBPtr = pdfHardGamAPtr = pdfHardGamBPtr = pdfUnresAPtr
    = pdfUnresBPtr = pdfUnresGamAPtr = pdfUnresGamBPtr = pdfVMDAPtr
    = pdfVMDBPtr = nullptr;

  // Switch off external PDF's by zero as input.
  if (!pdfAPtrIn) return true;

  // Save pointers.
  pdfAPtr       = pdfAPtrIn;
  // By default same pointers for hard-process PDF's.
  pdfHardAPtr   = pdfAPtrIn;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Routine to pass in pointers to PDF's. Usage optional.

bool BeamSetup::setPDFBPtr( PDFPtr pdfBPtrIn ) {

  // Reset pointers to be empty.
  pdfAPtr = pdfBPtr = pdfHardAPtr = pdfHardBPtr = pdfPomAPtr = pdfPomBPtr
    = pdfGamAPtr = pdfGamBPtr = pdfHardGamAPtr = pdfHardGamBPtr = pdfUnresAPtr
    = pdfUnresBPtr = pdfUnresGamAPtr = pdfUnresGamBPtr = pdfVMDAPtr
    = pdfVMDBPtr = nullptr;

  // Switch off external PDF's by zero as input.
  if (!pdfBPtrIn) return true;

  // Save pointers.
  pdfBPtr       = pdfBPtrIn;
  // By default same pointers for hard-process PDF's.
  pdfHardBPtr   = pdfBPtrIn;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Switch to new beam particle identities; for similar hadrons only.

bool BeamSetup::setBeamIDs( int idAIn, int idBIn) {

  // Do nothing if nothing changed.
  bool switchA = (idAIn != 0) && (idAIn != idA);
  bool switchB = (idBIn != 0) && (idBIn != idB);
  hasSwitchedIDs = switchA || switchB;
  if (!hasSwitchedIDs) return true;

  // Optionally perform checks to see if new are close relatives to old.
  // Empty for now. Could be based on below, to check that same PDF is used.

  // For allowIDAswitch on one may need to set a new PDF for A.
  // Note that cases are (have to be!) synchronized with the idAList order.
  int iPDFAnew = -1;
  if (switchA && allowIDAswitch) {
    switch ((abs(idAIn) / 10) % 1000) {
      case 11: case 21: iPDFAnew = 1; break;
      case 31: case 32: case 13: iPDFAnew = 2; break;
      case 22: iPDFAnew = (abs(idAIn) == 221) ? 3 : 1; break;
      case 33: iPDFAnew = (abs(idAIn) == 331) ? 4 : 5; break;
      case 41: case 42: iPDFAnew = 6; break;
      case 43: iPDFAnew = 7; break;
      case 44: iPDFAnew = 8; break;
      case 51: case 52: iPDFAnew = 9; break;
      case 53: iPDFAnew = 10; break;
      case 54: iPDFAnew = 11; break;
      case 55: iPDFAnew = 12; break;
      case 222: case 221: case 211: case 111: iPDFAnew = 0; break;
      case 311: case 312: case 321: case 322: case 213: iPDFAnew = 13; break;
      case 331: case 332: iPDFAnew = 14; break;
      case 333: iPDFAnew = 15; break;
      case 411: case 412: case 421: case 422: iPDFAnew = 16; break;
      case 431: case 413: case 432: case 423: iPDFAnew = 17; break;
      case 433: iPDFAnew = 18; break;
      case 511: case 512: case 521: case 522: iPDFAnew = 19; break;
      case 531: case 513: case 532: case 523: iPDFAnew = 20; break;
      case 533: iPDFAnew = 21; break;
    }

    // It should have worked, but error if not.
    if (iPDFAnew == -1 || iPDFAnew >= (int)pdfASavePtrs.size()) {
      loggerPtr->ERROR_MSG(
        "did not find PDF", "for idA = " + to_string(idAIn));
      switchA = false;
      if (!switchB) return false;
    }
  }

  // Store the new identities, also in Info.
  if (switchA) idA = idAIn;
  if (switchB) idB = idBIn;
  infoPtr->setBeamIDs( idA, idB);

  // Modify beam particles. Possibly also PDF for idA.
  if (switchA) {
    if (allowIDAswitch && iPDFAnew != iPDFAsave) {
      beamA.initPDFPtr( pdfASavePtrs[iPDFAnew], pdfASavePtrs[iPDFAnew]);
      iPDFAsave = iPDFAnew;
    }
    beamA.setBeamID( idA);
    beamA.initBeamKind();
  }
  if (switchB) {
    beamB.setBeamID( idB);
    beamB.initBeamKind();
  }

  return true;
}

//--------------------------------------------------------------------------

// Set beam CM energy.

bool BeamSetup::setKinematics(double eCMIn) {

  // Check that the frameType matches the input provided.
  if (frameType != 1) {
    loggerPtr->ABORT_MSG("input parameters do not match frame type");
    return false;
  }

  // Save input value.
  eCM = eCMIn;
  return true;

}

//--------------------------------------------------------------------------

// Set beam energies.

bool BeamSetup::setKinematics(double eAIn, double eBIn) {

  // Check that the frameType matches the input provided.
  if (frameType != 2) {
    loggerPtr->ABORT_MSG("input parameters do not match frame type");
    return false;
  }

  // Save input values.
  eA = eAIn;
  eB = eBIn;
  return true;

}

//--------------------------------------------------------------------------

// Set beam momenta.

bool BeamSetup::setKinematics(double pxAIn, double pyAIn, double pzAIn,
  double pxBIn, double pyBIn, double pzBIn) {

  // Check that the frameType matches the input provided.
  if (frameType != 3) {
    loggerPtr->ABORT_MSG("input parameters do not match frame type");
    return false;
  }

  // Save input values.
  pxA = pxAIn;
  pyA = pyAIn;
  pzA = pzAIn;
  pxB = pxBIn;
  pyB = pyBIn;
  pzB = pzBIn;
  return true;

}

//--------------------------------------------------------------------------

// Set beam momenta.

bool BeamSetup::setKinematics(Vec4 pAIn, Vec4 pBIn) {

  // Check that the frameType matches the input provided.
  if (frameType != 3) {
    loggerPtr->ABORT_MSG("input parameters do not match frame type");
    return false;
  }

  // Save input values.
  pxA = pAIn.px();
  pyA = pAIn.py();
  pzA = pAIn.pz();
  pxB = pBIn.px();
  pyB = pBIn.py();
  pzB = pBIn.pz();
  return true;

}

//--------------------------------------------------------------------------

// Set up frame of beams, Les Houches input, and switches for beam handling.

bool BeamSetup::initFrame() {

  // Find which frame type to use.
  frameType = mode("Beams:frameType");

  // Initialization with internal processes: read in and set values.
  doVarEcm       = false;
  allowIDAswitch = false;
  iPDFAsave      = 0;
  if (frameType < 4 ) {
    doLHA     = false;
    boostType = frameType;
    idA       = mode("Beams:idA");
    idB       = mode("Beams:idB");
    eCM       = parm("Beams:eCM");
    eA        = parm("Beams:eA");
    eB        = parm("Beams:eB");
    pxA       = parm("Beams:pxA");
    pyA       = parm("Beams:pyA");
    pzA       = parm("Beams:pzA");
    pxB       = parm("Beams:pxB");
    pyB       = parm("Beams:pyB");
    pzB       = parm("Beams:pzB");

    // Special option with variable incoming projectile.
    doVarEcm       = flag("Beams:allowVariableEnergy");
    allowIDAswitch = flag("Beams:allowIDAswitch");
    if (allowIDAswitch && !doVarEcm) {
      loggerPtr->ABORT_MSG(
        "allowed idA switch also requires Beams:allowVariableEnergy = on");
      return false;
    }

  // Initialization with a Les Houches Event File or an LHAup object.
  } else {
    doLHA     = true;
    boostType = 2;
    string lhef        = word("Beams:LHEF");
    string lhefHeader  = word("Beams:LHEFheader");
    bool   readHeaders = flag("Beams:readLHEFheaders");
    bool   setScales   = flag("Beams:setProductionScalesFromLHEF")
      || flag("Beams:setDipoleShowerStartingScalesFromLHEF");
    skipInit           = flag("Beams:newLHEFsameInit");
    int    nSkipAtInit = mode("Beams:nSkipLHEFatInit");

    // For file input: renew file stream or (re)new Les Houches object.
    if (frameType == 4) {
      const char* cstring1 = lhef.c_str();
      bool useExternal = (lhaUpPtr && !useNewLHA && lhaUpPtr->useExternal());
      if (!useExternal && useNewLHA && skipInit)
        lhaUpPtr->newEventFile(cstring1);
      else if (!useExternal) {
        // Header is optional, so use NULL pointer to indicate no value.
        const char* cstring2 = (lhefHeader == "void")
          ? nullptr : lhefHeader.c_str();
        lhaUpPtr = make_shared<LHAupLHEF>(infoPtr, cstring1, cstring2,
          readHeaders, setScales);
      }

      // Check that file was properly opened.
      if (!lhaUpPtr->fileFound()) {
        loggerPtr->ABORT_MSG("Les Houches Event File not found");
        return false;
      }

    // For object input: at least check that not null pointer.
    } else {
      if (lhaUpPtr == 0) {
        loggerPtr->ABORT_MSG("LHAup object not found");
        return false;
      }

      // LHAup object generic abort using fileFound() routine.
      if (!lhaUpPtr->fileFound()) {
        loggerPtr->ABORT_MSG("LHAup initialisation error");
        return false;
      }
    }

    // Send in pointer to info. Store or replace LHA pointer in other classes.
    lhaUpPtr->setPtr( infoPtr);

    // If second time around, only with new file, then simplify.
    // Optionally skip ahead a number of events at beginning of file.
    if (skipInit) {
      if (nSkipAtInit > 0) lhaUpPtr->skipEvent(nSkipAtInit);
      return true;
    }

    // Set LHAinit information (in some external program).
    if ( !lhaUpPtr->setInit()) {
      loggerPtr->ABORT_MSG("Les Houches initialization failed");
      return false;
    }

    // Extract beams from values set in an LHAinit object.
    idA = lhaUpPtr->idBeamA();
    idB = lhaUpPtr->idBeamB();
    int idRenameBeams = mode("LesHouches:idRenameBeams");
    if (abs(idA) == idRenameBeams) idA = 16;
    if (abs(idB) == idRenameBeams) idB = -16;
    if (idA == 0 || idB == 0) doProcessLevel = false;
    eA  = lhaUpPtr->eBeamA();
    eB  = lhaUpPtr->eBeamB();

    // Optionally skip ahead a number of events at beginning of file.
    if (nSkipAtInit > 0) lhaUpPtr->skipEvent(nSkipAtInit);
  }

  // Find out if beams are or have a resolved photon beam.
  // The PDF:lepton2gamma is kept for backwards compatibility, now
  // beamA2gamma and beamB2gamma are the master switches.
  bool lepton2gamma = flag("PDF:lepton2gamma");
  if (lepton2gamma && ( abs(idA) == 11 || abs(idA) == 13 || abs(idA) == 15 ))
      settingsPtr->flag("PDF:beamA2gamma", true);
  if (lepton2gamma && ( abs(idB) == 11 || abs(idB) == 13 || abs(idB) == 15 ))
      settingsPtr->flag("PDF:beamB2gamma", true);
  beamA2gamma = flag("PDF:beamA2gamma");
  beamB2gamma = flag("PDF:beamB2gamma");
  gammaMode   = mode("Photon:ProcessType");

  // Check if resolved photons are needed.
  beamAResGamma = (beamA2gamma || idA == 22)
    && ( (gammaMode == 1) || (gammaMode == 2) || (gammaMode == 0) );
  beamBResGamma = (beamB2gamma || idB == 22)
    && ( (gammaMode == 1) || (gammaMode == 3) || (gammaMode == 0) );

  // Check if unresolved photons are needed.
  beamAUnresGamma = (beamA2gamma || idA == 22)
    && ( (gammaMode == 4) || (gammaMode == 3) || (gammaMode == 0) );
  beamBUnresGamma = (beamB2gamma || idB == 22)
    && ( (gammaMode == 4) || (gammaMode == 2) || (gammaMode == 0) );

  // Check if VMD sampling is required for beam A and/or B.
  doDiffraction = flag("SoftQCD:all") || flag("SoftQCD:inelastic")
    || flag("SoftQCD:centralDiffractive") || flag("SoftQCD:singleDiffractive")
    || flag("SoftQCD:singleDiffractiveXB")
    || flag("SoftQCD:singleDiffractiveAX")
    || flag("SoftQCD:doubleDiffractive");
  doSoftQCD  = doDiffraction || flag("SoftQCD:elastic")
    || flag("SoftQCD:nonDiffractive");
  doHardDiff = flag("Diffraction:doHard");
  doVMDsideA = doSoftQCD && beamAResGamma;
  doVMDsideB = doSoftQCD && beamBResGamma;

  // Some other necessary setup.
  doProcessLevel   = flag("ProcessLevel:all");
  doMomentumSpread = flag("Beams:allowMomentumSpread");
  if (doVarEcm) doMomentumSpread = false;
  doVertexSpread   = flag("Beams:allowVertexSpread");
  doPartonVertex   = flag("PartonVertex:setVertex");
  doVertexPlane    = flag("PartonVertex:randomPlane");

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Initialize kinematics and PDFs of beams.

bool BeamSetup::initBeams(bool doNonPertIn, StringFlav* flavSelPtr) {

  // Info on process kinds.
  doNonPert = doNonPertIn;

  // Set up values related to beam shape.
  if (!beamShapePtr) beamShapePtr = make_shared<BeamShape>();
  beamShapePtr->init( *settingsPtr, rndmPtr);

  // Check that beams and beam combination can be handled.
  if (!checkBeams()) {
    loggerPtr->ABORT_MSG("checkBeams initialization failed");
    return false;
  }

  // Simplified beam setup when no process level.
  if (doNonPert && !doSoftQCD) {
    beamA.initID( idA);
    beamB.initID( idB);
    if (!initKinematics()) {
      loggerPtr->ABORT_MSG("kinematics initialization failed");
      return false;
    }
  }
  else if (!doProcessLevel) boostType = 1;

  // Full beam setup: first beam kinematics.
  else {

    if (!initKinematics()) {
      loggerPtr->ABORT_MSG("kinematics initialization failed");
      return false;
    }

    // Set up pointers to PDFs.
    if (!initPDFs()) {
      loggerPtr->ABORT_MSG("PDF initialization failed");
      return false;
    }

    // Set up the two beams and the common remnant system.
    beamA.init( idA, pzAcm, eA, mA, pdfAPtr, pdfHardAPtr,
      isUnresolvedA, flavSelPtr);
    beamB.init( idB, pzBcm, eB, mB, pdfBPtr, pdfHardBPtr,
      isUnresolvedB, flavSelPtr);

    // Special setup to allow switching between beam PDFs.
    if (allowIDAswitch) beamA.initSwitchID( pdfASavePtrs);

    // Pass information whether the beam will contain a photon beam.
    if (beamA2gamma) beamA.initGammaInBeam();
    if (beamB2gamma) beamB.initGammaInBeam();

    // Init also unresolved PDF pointers for photon beams when needed.
    if (beamAUnresGamma) beamA.initUnres( pdfUnresAPtr);
    if (beamBUnresGamma) beamB.initUnres( pdfUnresBPtr);

    // Optionally set up new alternative beams for these Pomerons.
    if ( doDiffraction || doHardDiff ) {
      beamPomA.init( 990,  0.5 * eCM, 0.5 * eCM, 0.,
        pdfPomAPtr, pdfPomAPtr, false, flavSelPtr);
      beamPomB.init( 990, -0.5 * eCM, 0.5 * eCM, 0.,
        pdfPomBPtr, pdfPomBPtr, false, flavSelPtr);
    }

    // Initialise VMD beams from gammas (in leptons). Use pion PDF for VMDs.
    if (doVMDsideA) beamVMDA.init( 111,  0.5 * eCM, 0.5 * eCM, 0.,
      pdfVMDAPtr, pdfVMDAPtr, false, flavSelPtr);
    if (doVMDsideB) beamVMDB.init( 111,  0.5 * eCM, 0.5 * eCM, 0.,
      pdfVMDBPtr, pdfVMDBPtr, false, flavSelPtr);

    // Optionally set up photon beams from lepton beams if resolved photons.
    if ( !(beamA.isGamma()) && beamA2gamma) {
      if ( gammaMode < 4 ) {
        beamGamA.init( 22,  0.5 * eCM, 0.5 * eCM, 0.,
          pdfGamAPtr, pdfHardGamAPtr, false, flavSelPtr);
      }
      if ( beamAUnresGamma ) beamGamA.initUnres( pdfUnresGamAPtr);
    }
    if ( !(beamB.isGamma()) && beamB2gamma) {
      if ( gammaMode < 4 ) {
        beamGamB.init( 22, -0.5 * eCM, 0.5 * eCM, 0.,
          pdfGamBPtr, pdfHardGamBPtr, false, flavSelPtr);
      }
      if ( beamBUnresGamma ) beamGamB.initUnres( pdfUnresGamBPtr);
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Clear all beams.

void BeamSetup::clear() {
  beamA.clear();
  beamB.clear();
  beamPomA.clear();
  beamPomB.clear();
  beamGamA.clear();
  beamGamB.clear();
  beamVMDA.clear();
  beamVMDB.clear();
}

//--------------------------------------------------------------------------

// Pick new beam valence flavours (for pi0, eta, K0S, Pomeron, etc.).

void BeamSetup::newValenceContent() {
  beamA.newValenceContent();
  beamB.newValenceContent();
  if ( doDiffraction || doHardDiff) {
    beamPomA.newValenceContent();
    beamPomB.newValenceContent();
  }
  if (doVMDsideA) beamVMDA.newValenceContent();
  if (doVMDsideB) beamVMDB.newValenceContent();
}

//--------------------------------------------------------------------------

// Recalculate kinematics for each event when beam momentum has a spread.

void BeamSetup::nextKinematics() {

  // Pick beam momentum spread and beam vertex. May be all that is needed.
  if (doMomentumSpread || doVertexSpread) beamShapePtr->pick();
  if (!doMomentumSpread && !doVarEcm) return;

  // Read out masses, since the particle id's may have changed.
  mA      = particleDataPtr->m0(idA);
  mB      = particleDataPtr->m0(idB);

  // Momentum spread: read out momentum shift to give current beam momenta.
  if (doMomentumSpread) {
    pAnow = pAinit + beamShapePtr->deltaPA();
    pAnow.e( sqrt(pAnow.pAbs2() + mA*mA) );
    pBnow = pBinit + beamShapePtr->deltaPB();
    pBnow.e( sqrt(pBnow.pAbs2() + mB*mB) );
    eCM   = (pAnow + pBnow).mCalc();

  // For variable energy in rest frame only need new eCM value, already set.
  } else if (frameType == 1) {

  // Variable energy but collinear beams: give current beam momenta.
  } else if (frameType == 2) {
    pAnow = Vec4( 0., 0.,  sqrtpos( eA*eA - mA*mA), eA);
    pBnow = Vec4( 0., 0., -sqrtpos( eB*eB - mB*mB), eB);
    eCM   = (pAnow + pBnow).mCalc();
    betaZ  = (pAnow.pz() + pBnow.pz()) / (eA + eB);
    gammaZ = (eA + eB) / eCM;

  // Variable three-momenta stored and energy calculated.
  } else if (frameType == 3) {
    pAnow = Vec4( pxA, pyA, pzA, sqrt(pxA*pxA + pyA*pyA + pzA*pzA + mA*mA) );
    pBnow = Vec4( pxB, pyB, pzB, sqrt(pxB*pxB + pyB*pyB + pzB*pzB + mB*mB) );
    eCM   = (pAnow + pBnow).mCalc();

  // Other possibilites not supported.
  } else {
    loggerPtr->ERROR_MSG("unsupported frameType");
    return;
  }

  // Construct CM frame kinematics.
  pzAcm = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
        * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  pzBcm = -pzAcm;
  eA    = sqrt(mA*mA + pzAcm*pzAcm);
  eB    = sqrt(mB*mB + pzBcm*pzBcm);

  // Set relevant info for other classes to use.
  infoPtr->setBeamA( idA, pzAcm, eA, mA);
  infoPtr->setBeamB( idB, pzBcm, eB, mB);
  infoPtr->setECM( eCM);
  beamA.newPzE( pzAcm, eA);
  beamB.newPzE( pzBcm, eB);

  // Set boost/rotation matrices from/to CM frame.
  if (frameType != 1) {
    MfromCM.reset();
    MfromCM.fromCMframe( pAnow, pBnow);
    MtoCM = MfromCM;
    MtoCM.invert();
  }

}

//--------------------------------------------------------------------------

// Boost from CM frame to lab frame, or inverse. Set production vertex.

void BeamSetup::boostAndVertex( Event& process, Event& event, bool toLab,
  bool setVertex) {

  // Optionally rotate event around its axis to randomize parton vertices.
  if (toLab && doPartonVertex && event.size() > 2) {
    if (process.size() > 2) {
      process[1].vProd( event[1].vProd() );
      process[2].vProd( event[2].vProd() );
    }
    if (doVertexPlane) {
      double phiVert = 2. * M_PI * rndmPtr->flat();
      process.rot( 0., phiVert);
      event.rot( 0., phiVert);
    }
  }

  // Boost process from CM frame to lab frame.
  if (toLab) {
    if      (boostType == 2) process.bst(0., 0., betaZ, gammaZ);
    else if (boostType == 3) process.rotbst(MfromCM);

    // Boost nonempty event from CM frame to lab frame.
    if (event.size() > 0) {
      if      (boostType == 2) event.bst(0., 0., betaZ, gammaZ);
      else if (boostType == 3) event.rotbst(MfromCM);
    }

  // Boost process from lab frame to CM frame.
  } else {
    if      (boostType == 2) process.bst(0., 0., -betaZ, gammaZ);
    else if (boostType == 3) process.rotbst(MtoCM);

    // Boost nonempty event from lab frame to CM frame.
    if (event.size() > 0) {
      if      (boostType == 2) event.bst(0., 0., -betaZ, gammaZ);
      else if (boostType == 3) event.rotbst(MtoCM);
    }
  }

  // Fix energy from mass and three-momentum, to patch up large boosts.
  for (int i = 1; i < event.size(); ++i)
    event[i].e( sqrtpos(event[i].m2() + event[i].pAbs2()) );

  // Set production vertex; assumes particles are in lab frame and at origin.
  if (setVertex && doVertexSpread) {
    Vec4 vertex = beamShapePtr->vertex();
    for (int i = 0; i < process.size(); ++i) process[i].vProdAdd( vertex);
    for (int i = 0; i < event.size(); ++i) event[i].vProdAdd( vertex);
  }

}

//--------------------------------------------------------------------------

// Check that beams and beam combination can be handled. Set up unresolved.

bool BeamSetup::checkBeams() {

  // Absolute flavours. If not to do process level then no check needed.
  int idAabs = abs(idA);
  int idBabs = abs(idB);
  if (!doProcessLevel) return true;

  // Special case for low-energy nonperturbative processes.
  if (doNonPert) {
    if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB)) {
      loggerPtr->ERROR_MSG("non-perturbative processes defined "
        "only for hadron-hadron collisions.");
      return false;
    }
    if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM) {
      loggerPtr->ERROR_MSG("beam particles have higher mass than eCM");
      return false;
    }
    return true;
  }

  // Neutrino beams always unresolved, charged lepton ones conditionally.
  bool isLeptonA    = (idAabs > 10 && idAabs < 17);
  bool isLeptonB    = (idBabs > 10 && idBabs < 17);
  bool isUnresLep   = !flag("PDF:lepton");
  bool isGammaA     = idAabs == 22;
  bool isGammaB     = idBabs == 22;
  isUnresolvedA     = (isLeptonA && isUnresLep);
  isUnresolvedB     = (isLeptonB && isUnresLep);

  // Also photons may be unresolved.
  if ( idAabs == 22 && !beamAResGamma ) isUnresolvedA = true;
  if ( idBabs == 22 && !beamBResGamma ) isUnresolvedB = true;

  // But not if resolved photons present.
  if ( beamAResGamma ) isUnresolvedA = false;
  if ( beamBResGamma ) isUnresolvedB = false;

  // Equate Dark Matter "beams" with incoming neutrinos.
  if (idAabs > 50 && idAabs < 61) isLeptonA = isUnresolvedA = true;
  if (idBabs > 50 && idBabs < 61) isLeptonB = isUnresolvedB = true;

  // Photon-initiated processes.
  if ( beamA2gamma || beamB2gamma || isGammaA || isGammaB ) {

    // No photon inside photon beams.
    if ( (beamA2gamma && isGammaA) || (beamB2gamma && isGammaB) ) {
      loggerPtr->ERROR_MSG("not possible to have a photon sub-beam"
         " within a photon beam");
      return false;
    }

    // Only gm+gm in lepton+lepton collisions.
    if ( isLeptonA && isLeptonB && ( !beamA2gamma || !beamB2gamma ) ) {
      loggerPtr->ERROR_MSG(
        "DIS with resolved photons currently not supported");
      return false;
    }

    // Photon beam and photon sub-beam not simultaneously allowed.
    if ( ( beamA2gamma && isGammaB ) || ( beamB2gamma && isGammaA ) ) {
      loggerPtr->ERROR_MSG("photoproduction together with pure photon "
        "beam currently not supported");
      return false;
    }

    // Allow soft QCD processes only when no direct photons present.
    bool isSoft = flag("SoftQCD:all") || flag("SoftQCD:nonDiffractive")
      || flag("SoftQCD:elastic") || flag("SoftQCD:singleDiffractive")
      || flag("SoftQCD:singleDiffractiveXB")
      || flag("SoftQCD:singleDiffractiveAX")
      || flag("SoftQCD:DoubleDiffractive")
      || flag("SoftQCD:CentralDiffractive") || flag("SoftQCD:inelastic");
    if (isSoft) {
      if ( ( (beamA2gamma || isGammaA) && !beamAResGamma )
        || ( (beamB2gamma || isGammaB) && !beamBResGamma ) ) {
        loggerPtr->ERROR_MSG("soft QCD only with resolved photons");
        return false;

      // Soft processes OK with resolved photons and hadrons.
      } else {
        return true;
      }

    // Otherwise OK.
    } else {
      return true;
    }
  }

  // Lepton-lepton collisions.
  if (isLeptonA && isLeptonB ) {

    // Lepton-lepton collisions OK (including neutrinos) if both (un)resolved
    if (isUnresolvedA == isUnresolvedB) return true;
  }

  // MBR model only implemented for pp/ppbar/pbarp collisions.
  int PomFlux     = mode("SigmaDiffractive:PomFlux");
  if (PomFlux == 5) {
    bool ispp       = (idAabs == 2212 && idBabs == 2212);
    bool ispbarpbar = (idA == -2212 && idB == -2212);
    if (ispp && !ispbarpbar) return true;
    loggerPtr->ERROR_MSG("cannot handle this beam combination with "
      "PomFlux == 5");
    return false;
  }

  // Hadron-hadron collisions OK, with Pomeron counted as hadron.
  bool isHadronA = particleDataPtr->isHadron(idA) || idA == 990;
  bool isHadronB = particleDataPtr->isHadron(idB) || idB == 990;
  int modeUnresolvedHadron = mode("BeamRemnants:unresolvedHadron");
  if (isHadronA && modeUnresolvedHadron%2 == 1) isUnresolvedA = true;
  if (isHadronB && modeUnresolvedHadron > 1)    isUnresolvedB = true;
  if (isHadronA && isHadronB) return true;

  // Lepton-hadron collisions OK for DIS processes or LHEF input,
  // although still primitive.
  if ( (isLeptonA && isHadronB) || (isHadronA && isLeptonB) ) {
    bool doDIS = flag("WeakBosonExchange:all")
      || flag("WeakBosonExchange:ff2ff(t:gmZ)")
      || flag("WeakBosonExchange:ff2ff(t:W)")
      || flag("Check:beams") || (frameType == 4);
    if (doDIS) return true;
  }

  // Allow to explicitly omit beam check for LHEF input.
  if ( mode("Beams:frameType") == 4 && !flag("Check:beams")) return true;

  // If no case above then failed.
  loggerPtr->ERROR_MSG("cannot handle this beam combination");
  return false;

}

//--------------------------------------------------------------------------

// Calculate kinematics at initialization. Store beam four-momenta.

bool BeamSetup::initKinematics() {

  // Find masses. Initial guess that we are in CM frame.
  mA       = particleDataPtr->m0(idA);
  mB       = particleDataPtr->m0(idB);
  betaZ    = 0.;
  gammaZ   = 1.;

  // Collinear beams not in CM frame: find CM energy.
  if (boostType == 2) {
    eA     = max(eA, mA);
    eB     = max(eB, mB);
    pzA    = sqrt(eA*eA - mA*mA);
    pzB    = -sqrt(eB*eB - mB*mB);
    pAinit = Vec4( 0., 0., pzA, eA);
    pBinit = Vec4( 0., 0., pzB, eB);
    eCM    = sqrt( pow2(eA + eB) - pow2(pzA + pzB) );

    // Find boost to rest frame.
    betaZ  = (pzA + pzB) / (eA + eB);
    gammaZ = (eA + eB) / eCM;
  }

  // Completely general beam directions: find CM energy.
  else if (boostType == 3) {
    eA     = sqrt( pxA*pxA + pyA*pyA + pzA*pzA + mA*mA);
    eB     = sqrt( pxB*pxB + pyB*pyB + pzB*pzB + mB*mB);
    pAinit = Vec4( pxA, pyA, pzA, eA);
    pBinit = Vec4( pxB, pyB, pzB, eB);
    eCM = (pAinit + pBinit).mCalc();

    // Find boost+rotation needed to move from/to CM frame.
    MfromCM.reset();
    MfromCM.fromCMframe( pAinit, pBinit);
    MtoCM = MfromCM;
    MtoCM.invert();
  }

  // Fail if CM energy below beam masses.
  if (eCM < mA + mB) {
    loggerPtr->ERROR_MSG("too low energy");
    return false;
  }

  // Set up CM-frame kinematics with beams along +-z axis.
  pzAcm    = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
  pzBcm    = -pzAcm;
  eA       = sqrt(mA*mA + pzAcm*pzAcm);
  eB       = sqrt(mB*mB + pzBcm*pzBcm);

  // If in CM frame then store beam four-vectors (else already done above).
  if (boostType != 2 && boostType != 3) {
    pAinit = Vec4( 0., 0., pzAcm, eA);
    pBinit = Vec4( 0., 0., pzBcm, eB);
  }

  // Store main info for access in process generation.
  infoPtr->setBeamA( idA, pzAcm, eA, mA);
  infoPtr->setBeamB( idB, pzBcm, eB, mB);
  infoPtr->setECM( eCM);

  // Must allow for generic boost+rotation when beam momentum spread.
  if (doMomentumSpread) boostType = 3;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Set up pointers to PDFs.

bool BeamSetup::initPDFs() {

  // Optionally set up photon PDF's for lepton -> gamma collisions. Done before
  // the main PDFs so that the gamma pointer can be used for the main PDF
  // (lepton). Both set also in case that only one of the photons is resolved.
  bool setupGammaBeams = ( (beamA2gamma || beamB2gamma) && (gammaMode < 4) );
  if (setupGammaBeams) {
    if ( beamA2gamma && pdfGamAPtr == 0 ) {
      pdfGamAPtr = getPDFPtr(22, 1, "A");
      if (!pdfGamAPtr->isSetup()) return false;

      // Set also unresolved photon beam when also unresolved photons.
      if (gammaMode != 1) {
        pdfUnresGamAPtr = getPDFPtr(22, 1, "A", false);
        if (!pdfUnresGamAPtr->isSetup()) return false;
      }

      // Set up optional hard photon PDF pointers.
      if (flag("PDF:useHard")) {
        pdfHardGamAPtr = getPDFPtr(22, 2);
        if (!pdfHardGamAPtr->isSetup()) return false;
      } else pdfHardGamAPtr = pdfGamAPtr;
    }
    if ( beamB2gamma && pdfGamBPtr == 0 ) {
      pdfGamBPtr = getPDFPtr(22, 1, "B");
      if (!pdfGamBPtr->isSetup()) return false;

      // Set also unresolved photon beam when also unresolved photons.
      if (gammaMode != 1) {
        pdfUnresGamBPtr = getPDFPtr(22, 1, "B", false);
        if (!pdfUnresGamBPtr->isSetup()) return false;
      }

      // Set up optional hard photon PDF pointers.
      if (flag("PDF:useHard")) {
        pdfHardGamBPtr = getPDFPtr(22, 2, "B");
        if (!pdfHardGamBPtr->isSetup()) return false;
      } else pdfHardGamBPtr = pdfGamBPtr;
    }
  }

  // Special setup for variable incoming idA hadron beam.
  if (allowIDAswitch) {
    pdfASavePtrs = vector<PDFPtr>(idAList.size());
    for (size_t iPA = 0; iPA < idAList.size(); ++iPA)
      pdfASavePtrs[iPA] = getPDFPtr( idAList[iPA], 1, "A" );
    pdfAPtr     = pdfASavePtrs[0];
    pdfBPtr     = getPDFPtr(idB, 1, "B");
    pdfHardAPtr = pdfAPtr;
    pdfHardBPtr = pdfBPtr;
    pdfPomAPtr  = getPDFPtr(990);
    pdfPomBPtr  = getPDFPtr(990);
    return true;
  }

  // Set up the PDF's, if not already done.
  if (pdfAPtr == 0) {
    pdfAPtr     = getPDFPtr(idA);
    if (pdfAPtr == 0 || !pdfAPtr->isSetup()) {
      loggerPtr->ERROR_MSG("could not set up PDF for beam A");
      return false;
    }
    pdfHardAPtr = pdfAPtr;
  }
  if (pdfBPtr == 0) {
    pdfBPtr     = getPDFPtr(idB, 1, "B");
    if (pdfBPtr == 0 || !pdfBPtr->isSetup()) {
      loggerPtr->ERROR_MSG("could not set up PDF for beam B");
      return false;
    }
    pdfHardBPtr = pdfBPtr;
  }

  // Optionally set up separate PDF's for hard process.
  if (flag("PDF:useHard")) {
    pdfHardAPtr = getPDFPtr(idA, 2);
    if (!pdfHardAPtr->isSetup()) return false;
    pdfHardBPtr = getPDFPtr(idB, 2, "B");
    if (!pdfHardBPtr->isSetup()) return false;
  }

  // Optionally use nuclear modifications for hard process PDFs.
  if (flag("PDF:useHardNPDFA")) {
    int idANucleus = mode("PDF:nPDFBeamA");
    pdfHardAPtr = getPDFPtr(idANucleus, 2, "A");
    if (!pdfHardAPtr->isSetup()) {
      loggerPtr->ERROR_MSG("could not set up nuclear PDF for beam A");
      return false;
    }
  }
  if (flag("PDF:useHardNPDFB")) {
    int idBNucleus = mode("PDF:nPDFBeamB");
    pdfHardBPtr = getPDFPtr(idBNucleus, 2, "B");
    if (!pdfHardBPtr->isSetup()) {
      loggerPtr->ERROR_MSG("could not set up nuclear PDF for beam B");
      return false;
    }
  }

  // Set up additional unresolved PDFs for photon beams when relevant.
  if ( (idA == 22 || beamA2gamma) && (gammaMode != 1 && gammaMode != 2) ) {
    if ( pdfUnresAPtr == 0 ) {
      pdfUnresAPtr = getPDFPtr(idA, 1, "A", false);
      if (!pdfUnresAPtr->isSetup()) return false;
    }
  }
  if ( (idB == 22 || beamB2gamma) && (gammaMode != 1 && gammaMode != 3) ) {
    if ( pdfUnresBPtr == 0 ) {
      pdfUnresBPtr = getPDFPtr(idB, 1, "B", false);
      if (!pdfUnresBPtr->isSetup()) return false;
    }
  }

  // Optionally set up Pomeron PDF's for diffractive physics.
  if ( doDiffraction || doHardDiff) {
    if (pdfPomAPtr == 0) {
      pdfPomAPtr    = getPDFPtr(990);
    }
    if (pdfPomBPtr == 0) {
      pdfPomBPtr    = getPDFPtr(990);
    }
  }

  // Optionally set up VMD PDF's for photon physics.
  if ( doSoftQCD && (doVMDsideA || doVMDsideB)) {
    if (pdfVMDAPtr == 0) {
      pdfVMDAPtr    = getPDFPtr(111);
    }
    if (pdfVMDBPtr == 0) {
      pdfVMDBPtr    = getPDFPtr(111);
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Create an LHAPDF plugin PDF.

PDFPtr BeamSetup::initLHAPDF(int idIn, string cfg) {

  // Determine the LHAPDF plugin library name.
  if (cfg.size() < 8) {
    loggerPtr->ERROR_MSG("invalid pSet " + cfg); return nullptr;}
  string className = cfg.substr(0, 7);
  if (className != "LHAPDF5" && className != "LHAPDF6") {
    loggerPtr->ERROR_MSG("invalid pSet " + cfg); return nullptr;}
  string libName = "libpythia8lhapdf" + cfg.substr(6, 1) + ".so";

  // Determine the PDF set name and member.
  string   setName = cfg.substr(8);
  int      setMem = 0;
  size_t   pos = setName.find_last_of("/");
  if (pos != string::npos) setMem = stoi(setName.substr(pos + 1));
  setName = setName.substr(0, pos);

  // Load and initialize the PDF.
  PDFPtr pdf = make_plugin<PDF>(
    libName, className, nullptr, settingsPtr, loggerPtr);
  if (pdf != nullptr && !pdf->init(idIn, setName, setMem, loggerPtr))
    return nullptr;
  else return pdf;

}

//--------------------------------------------------------------------------

// Routine to set up a PDF pointer.

PDFPtr BeamSetup::getPDFPtr(int idIn, int sequence, string beam,
  bool resolved) {

  // Temporary pointer to be returned.
  PDFPtr tempPDFPtr = 0;

  // Data file directory.
  string xmlPath = word( "xmlPath");
  string pdfdataPath = xmlPath.substr(0, xmlPath.length() - 7) + "pdfdata/";

  // One option is to treat a Pomeron like a pi0.
  if (idIn == 990 && word("PDF:PomSet") == "2") idIn = 111;

  // Check if photon beam inside proton.
  bool proton2gamma = (abs(idIn) == 2212) && ( ( beamA2gamma && (beam == "A") )
                    || ( beamB2gamma && (beam == "B") ) );

  // Nucleon-like beam, normal or hard choice.
  int idAbs = abs(idIn);
  int qContent = (idAbs / 10) % 1000;
  if ( ( qContent == 222 || qContent == 111
      || qContent == 221 || qContent == 211 ) && !proton2gamma ) {
    string pWord = word("PDF:p"
      + string(sequence == 1 ? "" : "Hard") + "Set"
      + string(beam == "A" ? "" : "B") ) ;
    if (pWord == "void" && sequence != 1 && beam == "B")
      pWord = word("PDF:pHardSet");
    if (pWord == "void") pWord = word("PDF:pSet");
    istringstream pStream(pWord);
    int pSet = 0;
    pStream >> pSet;

    // Use internal LHAgrid1 implementation for LHAPDF6 files.
    if (pSet == 0 && pWord.length() > 9
      && toLower(pWord).substr(0,9) == "lhagrid1:")
      tempPDFPtr = make_shared<LHAGrid1>
        (idIn, pWord, pdfdataPath, loggerPtr);

    // Use sets from LHAPDF.
    else if (pSet == 0) tempPDFPtr = initLHAPDF(idIn, pWord);

    // Use internal sets.
    else if (pSet == 1) tempPDFPtr = make_shared<GRV94L>(idIn);
    else if (pSet == 2) tempPDFPtr = make_shared<CTEQ5L>(idIn);
    else if (pSet <= 6)
      tempPDFPtr = make_shared<MSTWpdf>
        (idIn, pSet - 2, pdfdataPath, loggerPtr);
    else if (pSet <= 12)
      tempPDFPtr = make_shared<CTEQ6pdf>(idIn, pSet - 6, 1.,
        pdfdataPath, loggerPtr);
    else if (pSet <= 24)
      tempPDFPtr = make_shared<LHAGrid1>
        (idIn, pWord, pdfdataPath, loggerPtr);
    else tempPDFPtr = 0;
  }

  // Quasi-real photons inside a (anti-)proton beam.
  else if (proton2gamma) {

    // Find the resolved photon PDF to combine with the flux.
    PDFPtr tempGammaPDFPtr = nullptr;

    // Set up the combination of flux and PDF for resolved photons.
    if (resolved) {

      // Find the pre-set photon PDF, hard or normal.
      if (beam == "A") {
        tempGammaPDFPtr = (sequence == 1) ? pdfGamAPtr : pdfHardGamAPtr;
      } else {
        tempGammaPDFPtr = (sequence == 1) ? pdfGamBPtr : pdfHardGamBPtr;
      }
    }

    // Set the photon flux pointer and construct approximation.
    // Use the existing machinery for external fluxes.
    PDFPtr tempGammaFluxPtr = nullptr;
    double m2beam = pow2(particleDataPtr->m0(idIn));

    // Externally provided flux.
    if (mode("PDF:proton2gammaSet") == 0) {

      // Find the correct flux for given beam set with setPhotonFluxPtr().
      tempGammaFluxPtr = (beam == "A") ? pdfGamFluxAPtr : pdfGamFluxBPtr;

      // Check that external flux exist and complain if not.
      if (tempGammaFluxPtr == 0) {
        tempPDFPtr = 0;
        loggerPtr->ERROR_MSG(
          "no external photon flux provided with PDF:proton2gammaSet == 0",
          "for beam " + beam );
      }

    // Classic EPA proton by Budnev, Ginzburg, Meledin and Serbo.
    } else if (mode("PDF:proton2gammaSet") == 1) {

      // Check if Q^2 sampling on and turn off if necessary.
      if (flag("Photon:sampleQ2") == true ) {
        settingsPtr->flag("Photon:sampleQ2", false);
        loggerPtr->WARNING_MSG("photon virtuality sampling turned off as "
          "chosen flux is Q2 independent");
      }
      tempGammaFluxPtr = make_shared<ProtonPoint>(idIn, loggerPtr);

    // EPA approximation by Drees and Zeppenfeld.
    } else if (mode("PDF:proton2gammaSet") == 2) {
      tempGammaFluxPtr = make_shared<Proton2gammaDZ>(idIn);
    } else {
      loggerPtr->ERROR_MSG("invalid option for photon flux from proton");
    }

    // Construct flux object when pointer succesfully created.
    if ( tempGammaFluxPtr != 0) {
      tempPDFPtr = make_shared<EPAexternal>(idIn, m2beam, tempGammaFluxPtr,
        tempGammaPDFPtr, infoPtr, loggerPtr);
    } else {
      tempPDFPtr = 0;
      loggerPtr->ERROR_MSG("failed to set photon flux from protons");
    }
  }

  // Pion-like beam (or, in one option, Pomeron beam).
  else if ( qContent == 21 || qContent == 11
        || (qContent == 22 && idAbs != 221)) {
    string piWord = word("PDF:piSet"
                  + string(beam == "A" ? "" : "B") ) ;
    if (piWord == "void" && beam == "B") piWord = word("PDF:piSet");
    istringstream piStream(piWord);
    int piSet = 0;
    piStream >> piSet;

    // If VMD process then scale PDF accordingly:
    // f_a^VMD = alphaEM * (1/f_rho^2 + 1/f_omega^2 + 1/f_phi^2 + 1/f_J/psi)
    //         * f_a^pi0.
    // COR: New value here includes J/psi
    double rescale = (doVMDsideA || doVMDsideB) ? 0.0046549 : 1.;

    // Use internal LHAgrid1 implementation for LHAPDF6 files.
    if (piSet == 0 && piWord.length() > 9
      && toLower(piWord).substr(0,9) == "lhagrid1:")
      tempPDFPtr = make_shared<LHAGrid1>
        (idIn, piWord, pdfdataPath, loggerPtr);

    // Use sets from LHAPDF.
    else if (piSet == 0) tempPDFPtr = initLHAPDF(idIn, piWord);

    // Use internal set.
    else if (piSet == 1) tempPDFPtr = make_shared<GRVpiL>(idIn, rescale);
    else if (piSet == 2) tempPDFPtr = make_shared<GRSpiL>(idIn, rescale);
    else if (piSet == 3)
      tempPDFPtr = make_shared<LHAGrid1>(idIn, "lhagrid1:SU21piplus.dat",
        pdfdataPath, loggerPtr);
    else tempPDFPtr = nullptr;
  }

  // Pomeron beam, if not treated like a pi0 beam.
  else if (idIn == 990) {
    string pomWord = word("PDF:PomSet");
    double rescale = parm("PDF:PomRescale");
    istringstream pomStream(pomWord);
    int pomSet = 0;
    pomStream >> pomSet;

    // Use internal LHAgrid1 implementation for LHAPDF6 files.
    if (pomSet == 0 && pomWord.length() > 9
      && toLower(pomWord).substr(0,9) == "lhagrid1:")
      tempPDFPtr = make_shared<LHAGrid1>
        (idIn, pomWord, pdfdataPath, loggerPtr);

    // Use sets from LHAPDF.
    else if (pomSet == 0) tempPDFPtr = initLHAPDF(idIn, pomWord);

    // A generic Q2-independent parametrization.
    else if (pomSet == 1) {
      double gluonA      = parm("PDF:PomGluonA");
      double gluonB      = parm("PDF:PomGluonB");
      double quarkA      = parm("PDF:PomQuarkA");
      double quarkB      = parm("PDF:PomQuarkB");
      double quarkFrac   = parm("PDF:PomQuarkFrac");
      double strangeSupp = parm("PDF:PomStrangeSupp");
      tempPDFPtr = make_shared<PomFix>( 990, gluonA, gluonB, quarkA, quarkB,
        quarkFrac, strangeSupp);
    }

    // The H1 Q2-dependent parametrizations. Initialization requires files.
    else if (pomSet == 3 || pomSet == 4) tempPDFPtr =
      make_shared<PomH1FitAB>( 990, pomSet - 2, rescale, pdfdataPath,
        loggerPtr);
    else if (pomSet == 5) tempPDFPtr =
      make_shared<PomH1Jets>( 990, 1, rescale, pdfdataPath, loggerPtr);
    else if (pomSet == 6) tempPDFPtr =
      make_shared<PomH1FitAB>( 990, 3, rescale, pdfdataPath, loggerPtr);
    // The parametrizations of Alvero, Collins, Terron and Whitmore.
    else if (pomSet > 6 && pomSet < 11)  {
      tempPDFPtr = make_shared<CTEQ6pdf>( 990, pomSet + 4, rescale,
        pdfdataPath, loggerPtr);
      loggerPtr->WARNING_MSG("pomeron flux parameters forced for ACTW PDFs");
      settingsPtr->mode("SigmaDiffractive:PomFlux", 4);
      double pomFluxEps = (pomSet == 10) ? 0.19 : 0.14;
      settingsPtr->parm("SigmaDiffractive:PomFluxEpsilon", pomFluxEps);
      settingsPtr->parm("SigmaDiffractive:PomFluxAlphaPrime", 0.25);
    }
    else if (pomSet == 11 ) tempPDFPtr =
      make_shared<PomHISASD>(990, getPDFPtr(2212), *settingsPtr, loggerPtr);
    else if (pomSet >= 12 && pomSet <= 15) tempPDFPtr =
      make_shared<LHAGrid1>(idIn, "1" + pomWord, pdfdataPath, loggerPtr);
    else tempPDFPtr = 0;
  }

  // Set up nuclear PDFs.
  else if (idIn > 100000000) {

    // Which nPDF set to use.
    int nPDFSet = (beam == "B") ? mode("PDF:nPDFSetB")
                                : mode("PDF:nPDFSetA");

    // Temporary pointer for storing proton PDF pointer.
    PDFPtr tempProtonPDFPtr = (beam == "B") ? pdfHardBPtr : pdfHardAPtr;
    if (nPDFSet == 0)
      tempPDFPtr = make_shared<Isospin>(idIn, tempProtonPDFPtr);
    else if (nPDFSet == 1 || nPDFSet == 2)
      tempPDFPtr = make_shared<EPS09>(idIn, nPDFSet, 1, pdfdataPath,
        tempProtonPDFPtr, loggerPtr);
    else if (nPDFSet == 3)
      tempPDFPtr = make_shared<EPPS16>(idIn, 1, pdfdataPath,
        tempProtonPDFPtr, loggerPtr);
    else tempPDFPtr = 0;
  }

  // Photon beam, either point-like (unresolved) or resolved.
  else if (abs(idIn) == 22) {

    // For unresolved beam use the point-like PDF.
    if (!resolved) {
      tempPDFPtr = make_shared<GammaPoint>(idIn);
    } else {
      int gammaSet = mode("PDF:GammaSet");

      // Point-like beam if unresolved photons.
      bool beamIsPoint
        = ( !beamAResGamma && beamAUnresGamma && !(beam == "B") )
        || ( !beamBResGamma && beamBUnresGamma && (beam == "B") );

      // Use different PDFs for hard process.
      if ( sequence == 2) {

        // Find the name or number of the hard PDF set.
        string gmWord = word("PDF:GammaHardSet");
        int gmSet     = 0;
        if (gmWord == "void") gmSet = mode("PDF:GammaSet");
        else {
          istringstream gmStream(gmWord);
          gmStream >> gmSet;
        }

        // Use sets from LHAPDF. Only available for hard processes.
        if (gmSet == 0 && !beamIsPoint) {
          tempPDFPtr = initLHAPDF(idIn, gmWord);
          return tempPDFPtr;
        }

        // Or set up an internal set (though currently only one).
        gammaSet = gmSet;
      }

      // Set up the PDF.
      if      (beamIsPoint)   tempPDFPtr = make_shared<GammaPoint>(idIn);
      else if (gammaSet == 1) tempPDFPtr = make_shared<CJKL>(idIn, rndmPtr);
      else                    tempPDFPtr = 0;
    }
  }

  // Lepton beam: neutrino, resolved charged lepton or unresolved ditto.
  // Also photon inside lepton PDFs.
  else if (abs(idIn) > 10 && abs(idIn) < 17) {

    // For neutrinos only point-like PDF.
    if (abs(idIn)%2 == 0) {
      tempPDFPtr = make_shared<NeutrinoPoint>(idIn);

    // Set up resolved photon inside lepton for beam A.
    } else if ( beamAResGamma && (beam == "A") && resolved ) {

      // Find the pre-set photon PDF, hard or normal.
      PDFPtr tempGammaPDFPtr = 0;
      if ( sequence == 2) tempGammaPDFPtr = pdfHardGamAPtr;
      else                tempGammaPDFPtr = pdfGamAPtr;

      // Get the mass of lepton and maximum virtuality of the photon.
      double m2beam     = pow2(particleDataPtr->m0(idIn));
      double Q2maxGamma = parm("Photon:Q2max");

      // Initialize the gamma-inside-lepton PDFs with internal photon flux.
      if (mode("PDF:lepton2gammaSet") == 1) {
        tempPDFPtr = make_shared<Lepton2gamma>(idIn, m2beam, Q2maxGamma,
          tempGammaPDFPtr, infoPtr);

      // Initialize the gamma-inside-lepton PDFs with external photon flux.
      // Requires that the pointer to the flux set.
      } else if ( mode("PDF:lepton2gammaSet") == 0 ) {
        PDFPtr tempGammaFluxPtr = pdfGamFluxAPtr;
        if ( tempGammaFluxPtr != 0)
          tempPDFPtr = make_shared<EPAexternal>(idIn, m2beam,
            tempGammaFluxPtr, tempGammaPDFPtr, infoPtr, loggerPtr);
        else {
          tempPDFPtr = 0;
          loggerPtr->ERROR_MSG(
            "no external photon flux provided with PDF:lepton2gammaSet == 0");
        }
      } else tempPDFPtr = 0;

    // Set up resolved photon inside lepton for beam B.
    } else if ( beamBResGamma && (beam == "B") && resolved ) {

      // Find the pre-set photon PDF, hard or normal.
      PDFPtr tempGammaPDFPtr = 0;
      if ( sequence == 2) tempGammaPDFPtr = pdfHardGamBPtr;
      else                tempGammaPDFPtr = pdfGamBPtr;

      // Get the mass of lepton and maximum virtuality of the photon.
      double m2beam     = pow2(particleDataPtr->m0(idIn));
      double Q2maxGamma = parm("Photon:Q2max");

      // Initialize the gamma-inside-lepton PDFs with internal photon flux.
      if (mode("PDF:lepton2gammaSet") == 1) {
        tempPDFPtr = make_shared<Lepton2gamma>(idIn, m2beam, Q2maxGamma,
          tempGammaPDFPtr, infoPtr);

      // Initialize the gamma-inside-lepton PDFs with external photon flux.
      } else if ( mode("PDF:lepton2gammaSet") == 0 ) {
        PDFPtr tempGammaFluxPtr = pdfGamFluxBPtr;
        if ( tempGammaFluxPtr != 0)
          tempPDFPtr = make_shared<EPAexternal>(idIn, m2beam,
            tempGammaFluxPtr, tempGammaPDFPtr, infoPtr, loggerPtr);
        else {
          tempPDFPtr = 0;
          loggerPtr->ERROR_MSG(
            "no external photon flux provided with PDF:lepton2gammaSet == 0");
        }
      } else tempPDFPtr = 0;

    // Usual lepton PDFs.
    } else if (flag("PDF:lepton")) {
      double m2beam = pow2(particleDataPtr->m0(idIn));
      double Q2maxGamma = parm("Photon:Q2max");
      if (mode("PDF:lepton2gammaSet") == 1 ) {
        tempPDFPtr = make_shared<Lepton>(idIn, Q2maxGamma, infoPtr);

      // External photon flux for direct-photon processes.
      } else if (mode("PDF:lepton2gammaSet") == 0 ) {
        PDFPtr tempGammaPDFPtr;
        PDFPtr tempGammaFluxPtr = (beam == "B") ?
          pdfGamFluxBPtr : pdfGamFluxAPtr;
        if ( tempGammaFluxPtr != 0) tempPDFPtr =
          make_shared<EPAexternal>(idIn, m2beam,
            tempGammaFluxPtr, tempGammaPDFPtr, infoPtr, loggerPtr);
        else {
          tempPDFPtr = 0;
          loggerPtr->ERROR_MSG(
            "no external photon flux provided with PDF:lepton2gammaSet == 0");
        }
      } else tempPDFPtr = 0;
    }
    else tempPDFPtr = make_shared<LeptonPoint>(idIn);

  // Dark matter beam set up as pointlike lepton.
  } else if (abs(idIn) > 50 && abs(idIn) < 60) {
    tempPDFPtr = make_shared<LeptonPoint>(idIn);

  // Further hadronic beams.
  } else if (particleDataPtr->isHadron(idIn)) {
    string baseParticle;
    switch ((idAbs / 10) % 1000) {
      case 11: case 21: baseParticle = "piplus"; break;
      case 22: baseParticle = "eta"; break;
      case 33: baseParticle = (idAbs == 331) ? "eta" : "phi"; break;
      case 31: case 32: case 13: baseParticle = "Kplus"; break;
      case 41: case 42: baseParticle = "Dzero"; break;
      case 43: baseParticle = "Dsplus"; break;
      case 44: baseParticle = "Jpsi"; break;
      case 51: case 52: baseParticle = "Bplus"; break;
      case 53: baseParticle = "Bszero"; break;
      case 54: baseParticle = "Bcplus"; break;
      case 55: baseParticle = "Upsilon"; break;
      case 222: case 221: case 211: case 111: baseParticle = "proton"; break;
      case 311: case 312: case 321: case 322: case 213: baseParticle
        = "Sigmaplus"; break;
      case 331: case 332: baseParticle = "Xizero"; break;
      case 333: baseParticle = "Omega"; break;
      case 411: case 412: case 421: case 422: baseParticle = "Sigmacplusplus";
        break;
      case 431: case 413: case 432: case 423: baseParticle = "Xicplus"; break;
      case 433: baseParticle = "Omegac"; break;
      case 511: case 512: case 521: case 522: baseParticle = "Sigmabplus";
        break;
      case 531: case 513: case 532: case 523: baseParticle = "Xibzero"; break;
      case 533: baseParticle = "Omegab"; break;
    }
    tempPDFPtr = make_shared<LHAGrid1>(idIn,
      "lhagrid1:SU21"+baseParticle+".dat", pdfdataPath, loggerPtr);
  }

  // Optionally allow extrapolation beyond x and Q2 limits.
  if (tempPDFPtr)
    tempPDFPtr->setExtrapolate( flag("PDF:extrapolate") );

  // Done.
  return tempPDFPtr;
}

//--------------------------------------------------------------------------

// Return a map of the PDF pointers.

map<string, PDFPtr> BeamSetup::getPDFPtr() {
  return {
    {"A", pdfAPtr}, {"B", pdfBPtr},
    {"HardA", pdfHardAPtr}, {"HardB", pdfHardBPtr},
    {"PomA", pdfPomAPtr}, {"PomB", pdfPomBPtr},
    {"GamA", pdfGamAPtr}, {"GamB", pdfGamBPtr},
    {"HardGamA", pdfHardGamAPtr}, {"HardGamB", pdfHardGamBPtr},
    {"UnresA", pdfUnresAPtr}, {"UnresB", pdfUnresBPtr},
    {"UnresGamA", pdfUnresGamAPtr}, {"UnresGamB", pdfUnresGamBPtr},
    {"VMDA", pdfVMDAPtr}, {"VMDB", pdfVMDBPtr}};
}

//==========================================================================

} // end namespace Pythia8
