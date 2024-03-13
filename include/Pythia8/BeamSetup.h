// BeamSetup.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains a helper class for the setup of beam flavour,
// kinematics and PDFs.

#ifndef Pythia8_BeamSetup_H
#define Pythia8_BeamSetup_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/BeamShape.h"
#include "Pythia8/HadronLevel.h"
#include "Pythia8/Info.h"
#include "Pythia8/LesHouches.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PartonLevel.h"
#include "Pythia8/PhysicsBase.h"
#include "Pythia8/ProcessLevel.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The BeamSetup class contains a number of routines auxiliary to Pythia
// to set up beam flavour, kinematics and PDFs.

class BeamSetup : public PhysicsBase {

public:

  // Constructor.
  BeamSetup() = default;

  // Possibility to pass in pointers to PDF's.
  bool setPDFPtr( PDFPtr pdfAPtrIn, PDFPtr pdfBPtrIn,
    PDFPtr pdfHardAPtrIn = nullptr, PDFPtr pdfHardBPtrIn = nullptr,
    PDFPtr pdfPomAPtrIn = nullptr, PDFPtr pdfPomBPtrIn = nullptr,
    PDFPtr pdfGamAPtrIn = nullptr, PDFPtr pdfGamBPtrIn = nullptr,
    PDFPtr pdfHardGamAPtrIn = nullptr, PDFPtr pdfHardGamBPtrIn = nullptr,
    PDFPtr pdfUnresAPtrIn = nullptr, PDFPtr pdfUnresBPtrIn = nullptr,
    PDFPtr pdfUnresGamAPtrIn = nullptr, PDFPtr pdfUnresGamBPtrIn = nullptr,
    PDFPtr pdfVMDAPtrIn = nullptr, PDFPtr pdfVMDBPtrIn = nullptr);
  bool setPDFAPtr( PDFPtr pdfAPtrIn );
  bool setPDFBPtr( PDFPtr pdfBPtrIn );

  // Set photon fluxes externally. Used with option "PDF:lepton2gammaSet = 2".
  bool setPhotonFluxPtr( PDFPtr photonFluxAIn, PDFPtr photonFluxBIn) {
    if ( photonFluxAIn ) pdfGamFluxAPtr = photonFluxAIn;
    if ( photonFluxBIn ) pdfGamFluxBPtr = photonFluxBIn;
    return true;}

  // Possibility to pass in pointer to external LHA-interfaced generator.
  bool setLHAupPtr( LHAupPtr lhaUpPtrIn) {lhaUpPtr = lhaUpPtrIn;
    useNewLHA = false; return true;}

  // Switch to new beam particle identities; for similar hadrons only.
  bool setBeamIDs( int idAin, int idBin = 0);

  // Switch beam kinematics.
  bool setKinematics(double eCMIn);
  bool setKinematics(double eAIn, double eBIn);
  bool setKinematics(double pxAIn, double pyAIn, double pzAIn,
                     double pxBIn, double pyBIn, double pzBIn);
  bool setKinematics(Vec4 pAIn, Vec4 pBIn);

  // Possibility to pass in pointer for beam shape.
  bool setBeamShapePtr( BeamShapePtr beamShapePtrIn) {
    beamShapePtr = beamShapePtrIn; return true;}

  // Possibility to access the pointer to the BeamShape object.
  BeamShapePtr getBeamShapePtr() { return beamShapePtr; }

  // Return a parton density set among list of possibilities.
  PDFPtr getPDFPtr(int idIn, int sequence = 1, string beam = "A",
    bool resolved = true);

  // Return a map of the PDF pointers.
  map<string, PDFPtr> getPDFPtr();

  // Set up frame of beams, Les Houches input, and switches for beam handling.
  bool initFrame();

  // Initialize kinematics and PDFs of beams.
  bool initBeams(bool doNonPertIn, StringFlav* flavSelPtr);

  // Return whether VMD states sampled.
  bool getVMDsideA() { return doVMDsideA; }
  bool getVMDsideB() { return doVMDsideB; }

  // Clear all beams.
  void clear();

  // Pick new beam valence flavours (for pi0, eta, K0S, Pomeron, etc.).
  void newValenceContent();

  // Recalculate kinematics for each event when beam momentum has a spread.
  void nextKinematics();

  // Boost from CM frame to lab frame, or inverse. Set production vertex.
  void boostAndVertex( Event& process, Event& event, bool toLab,
    bool setVertex);

  // Print parton lists for the main beams. For debug mainly.
  void list() const { beamA.list(); beamB.list(); }

  // Some data values are kept public so that the Pythia class can access them.
  bool   doLHA = false, useNewLHA = false, skipInit = false,
         doMomentumSpread = {}, doVertexSpread = {}, doVarEcm = {},
         allowIDAswitch = {}, hasSwitchedIDs = {}, beamA2gamma = {},
         beamB2gamma = {};
  int    idA = {}, idB = {}, frameType = {}, boostType = {}, iPDFAsave = {},
         gammaMode = {};
  double mA = {}, mB = {}, pxA = {}, pxB = {}, pyA = {}, pyB = {}, pzA = {},
         pzB = {}, eA = {}, eB = {}, pzAcm = {}, pzBcm = {}, eCM = {},
         betaZ = {}, gammaZ = {};
  Vec4   pAinit = {}, pBinit = {}, pAnow = {}, pBnow = {};
  RotBstMatrix MfromCM = {}, MtoCM = {};
  LHAupPtr lhaUpPtr = {};

  // The two incoming beams.
  BeamParticle   beamA = {};
  BeamParticle   beamB = {};

  // Alternative Pomeron beam-inside-beam.
  BeamParticle beamPomA = {};
  BeamParticle beamPomB = {};

  // Alternative photon beam-inside-beam.
  BeamParticle beamGamA = {};
  BeamParticle beamGamB = {};

  // Alternative VMD beam-inside-beam.
  BeamParticle beamVMDA = {};
  BeamParticle beamVMDB = {};

  // Hadron types for rapid switching.
  vector<int> idAList = { 2212, 211, 311, 221,
         331, 333, 411, 431, 443, 511, 531, 541, 553, 3212, 3312, 3334,
         4112, 4312, 4332, 5112, 5312, 5332};

protected:

  void onInitInfoPtr() override {
    registerSubObject(beamA);
    registerSubObject(beamB);
    registerSubObject(beamPomA);
    registerSubObject(beamPomB);
    registerSubObject(beamGamA);
    registerSubObject(beamGamB);
    registerSubObject(beamVMDA);
    registerSubObject(beamVMDB);
  }

private:

  // Initialization data, plus some event-specific.
  bool   doNonPert = {}, doDiffraction = {}, doSoftQCD = {},
         doHardDiff = {}, doProcessLevel = {}, doPartonVertex = {},
         doVertexPlane = {}, isUnresolvedA = {}, isUnresolvedB = {},
         doVMDsideA = {}, doVMDsideB = {}, beamAResGamma = {},
         beamBResGamma = {}, beamAUnresGamma = {}, beamBUnresGamma = {};

  // Pointers to the PDFs of beams, with several alternatives.
  PDFPtr pdfAPtr = {}, pdfBPtr = {}, pdfHardAPtr = {}, pdfHardBPtr = {},
         pdfPomAPtr = {}, pdfPomBPtr = {}, pdfGamAPtr = {}, pdfGamBPtr = {},
         pdfHardGamAPtr = {}, pdfHardGamBPtr = {}, pdfUnresAPtr = {},
         pdfUnresBPtr = {}, pdfUnresGamAPtr = {}, pdfUnresGamBPtr = {},
         pdfGamFluxAPtr = {}, pdfGamFluxBPtr = {}, pdfVMDAPtr = {},
         pdfVMDBPtr = {};

  // Array of PDFs to be used when idA can be changed between events.
  vector<PDFPtr> pdfASavePtrs = {};

  // Pointer to BeamShape object for beam momentum and interaction vertex.
  BeamShapePtr  beamShapePtr = {};

  // Check that beams and beam combination can be handled.
  bool checkBeams();

  // Calculate kinematics at initialization.
  bool initKinematics();

  // Set up pointers to PDFs.
  bool initPDFs();

  // Create an LHAPDF plugin PDF.
  PDFPtr initLHAPDF(int idIn, string cfg);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_BeamSetup_H
