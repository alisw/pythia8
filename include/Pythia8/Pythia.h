// Pythia.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class for event generation.
// Pythia: provide the main user interface to everything else.

#ifndef Pythia8_Pythia_H
#define Pythia8_Pythia_H

// Version number defined for use in macros and for consistency checks.
#define PYTHIA_VERSION 8.315
#define PYTHIA_VERSION_INTEGER 8315

// Header files for the Pythia class and for what else the user may need.
#include "Pythia8/Analysis.h"
#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/BeamSetup.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/HadronLevel.h"
#include "Pythia8/HadronWidths.h"
#include "Pythia8/Info.h"
#include "Pythia8/JunctionSplitting.h"
#include "Pythia8/LesHouches.h"
#include "Pythia8/Logger.h"
#include "Pythia8/SigmaLowEnergy.h"
#include "Pythia8/Merging.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/PartonLevel.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonDistributions.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/PhysicsBase.h"
#include "Pythia8/ProcessLevel.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/RHadrons.h"
#include "Pythia8/Ropewalk.h"
#include "Pythia8/Settings.h"
#include "Pythia8/ShowerModel.h"
#include "Pythia8/SigmaTotal.h"
#include "Pythia8/SimpleSpaceShower.h"
#include "Pythia8/SimpleTimeShower.h"
#include "Pythia8/SpaceShower.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/StringInteractions.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/SLHAinterface.h"
#include "Pythia8/TimeShower.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/VinciaCommon.h"
#include "Pythia8/Weights.h"

namespace Pythia8 {

//==========================================================================

// Forward declaration of the HeavyIons and HIUserHooks classes.
class HeavyIons;
class HIUserHooks;

// Forward declaration of PythiaParallel class, to be friended.
class PythiaParallel;

// The Pythia class contains the top-level routines to generate an event.

class Pythia {

public:

  // Constructor. (See Pythia.cc file.)
  Pythia(string xmlDir = "../share/Pythia8/xmldoc", bool printBanner = true);

  // Constructor to copy settings and particle database from another Pythia
  // object instead of XML files (to speed up multiple initialisations).
  Pythia(Settings& settingsIn, ParticleData& particleDataIn,
    bool printBanner = true);

  // Constructor taking input from streams instead of files.
  Pythia( istream& settingsStrings, istream& particleDataStrings,
    bool printBanner = true);

  // Destructor.
  ~Pythia() {}

  // Copy and = constructors cannot be used.
  Pythia(const Pythia&) = delete;
  Pythia& operator=(const Pythia&) = delete;

  // Check consistency of version numbers (called by constructors).
  bool checkVersion();

  // Read in one update for a setting or particle data from a single line.
  inline bool readString(string line, bool warn = true,
    int subrun = SUBRUNDEFAULT) {
    return isConstructed ? settings.readString(line, warn, subrun) : false;}

  // Read in updates for settings or particle data from user-defined file.
  inline bool readFile(string fileName, bool warn = true,
    int subrun = SUBRUNDEFAULT) {
    return isConstructed ? settings.readFile(fileName, warn, subrun) : false;}
  inline bool readFile(string fileName, int subrun) {
    return readFile(fileName, true, subrun);}
  inline bool readFile(istream& is = cin, bool warn = true,
    int subrun = SUBRUNDEFAULT) {
    return isConstructed ? settings.readFile(is, warn, subrun) : false;}
  inline bool readFile(istream& is, int subrun) {
    return readFile(is, true, subrun);}

  // Possibility to pass in pointers to PDF's.
  inline bool setPDFPtr(PDFPtr pdfAPtrIn, PDFPtr pdfBPtrIn,
    PDFPtr pdfHardAPtrIn = nullptr, PDFPtr pdfHardBPtrIn = nullptr,
    PDFPtr pdfPomAPtrIn = nullptr, PDFPtr pdfPomBPtrIn = nullptr,
    PDFPtr pdfGamAPtrIn = nullptr, PDFPtr pdfGamBPtrIn = nullptr,
    PDFPtr pdfHardGamAPtrIn = nullptr, PDFPtr pdfHardGamBPtrIn = nullptr,
    PDFPtr pdfUnresAPtrIn = nullptr, PDFPtr pdfUnresBPtrIn = nullptr,
    PDFPtr pdfUnresGamAPtrIn = nullptr, PDFPtr pdfUnresGamBPtrIn = nullptr,
    PDFPtr pdfVMDAPtrIn = nullptr, PDFPtr pdfVMDBPtrIn = nullptr) {
      return beamSetup.setPDFPtr( pdfAPtrIn, pdfBPtrIn, pdfHardAPtrIn,
      pdfHardBPtrIn, pdfPomAPtrIn, pdfPomBPtrIn, pdfGamAPtrIn, pdfGamBPtrIn,
      pdfHardGamAPtrIn, pdfHardGamBPtrIn, pdfUnresAPtrIn, pdfUnresBPtrIn,
      pdfUnresGamAPtrIn, pdfUnresGamBPtrIn, pdfVMDAPtrIn, pdfVMDBPtrIn); }
  inline bool setPDFAPtr( PDFPtr pdfAPtrIn ) {
    return beamSetup.setPDFAPtr( pdfAPtrIn); }
  inline bool setPDFBPtr( PDFPtr pdfBPtrIn ) {
    return beamSetup.setPDFBPtr( pdfBPtrIn); }

  // Set photon fluxes externally. Used with option "PDF:lepton2gammaSet = 2".
  inline bool setPhotonFluxPtr(PDFPtr photonFluxAIn, PDFPtr photonFluxBIn) {
    return beamSetup.setPhotonFluxPtr( photonFluxAIn, photonFluxBIn); }

  // Possibility to pass in pointer to external LHA-interfaced generator.
  inline bool setLHAupPtr(LHAupPtr lhaUpPtrIn) {
    lhaUpPtr = lhaUpPtrIn;
    useNewLHA = false;
    return beamSetup.setLHAupPtr(lhaUpPtrIn);}

  // Possibility to pass in pointer for external handling of some decays.
  inline bool setDecayPtr(DecayHandlerPtr decayHandlePtrIn,
    vector<int> handledParticlesIn = {}) {
    decayHandlePtr = decayHandlePtrIn;
    handledParticles = handledParticlesIn.size() == 0 ?
      decayHandlePtrIn->handledParticles() : handledParticlesIn; return true;}

  // Possibility to pass in pointer for external random number generation.
  inline bool setRndmEnginePtr(RndmEnginePtr rndmEnginePtrIn) {
    return rndm.rndmEnginePtr(rndmEnginePtrIn);}

  // Possibility to pass in pointer for user hooks.
  inline bool setUserHooksPtr(UserHooksPtr userHooksPtrIn) {
    userHooksPtr = userHooksPtrIn; return true;}

  // Possibility to add further pointers to allow multiple user hooks.
  inline bool addUserHooksPtr(UserHooksPtr userHooksPtrIn) {
    if ( !userHooksPtrIn ) return false;
    if ( !userHooksPtr ) return setUserHooksPtr(userHooksPtrIn);
    shared_ptr<UserHooksVector> uhv =
      dynamic_pointer_cast<UserHooksVector>(userHooksPtr);
    if ( !uhv ) { uhv = make_shared<UserHooksVector>();
      uhv->hooks.push_back(userHooksPtr); userHooksPtr = uhv; }
    uhv->hooks.push_back(userHooksPtrIn); return true;}

  // Possibility to insert a user hook.
  inline bool insertUserHooksPtr(int idx, UserHooksPtr userHooksPtrIn) {
    if ( !userHooksPtrIn || !userHooksPtr ) return false;
    shared_ptr<UserHooksVector> uhv =
      dynamic_pointer_cast<UserHooksVector>(userHooksPtr);
    if ( !uhv || idx < 0 || idx > (int)uhv->hooks.size() ) return false;
    uhv->hooks.insert(uhv->hooks.begin() + idx, userHooksPtrIn); return true;}

  // Possibility to pass in pointer for full merging class.
  inline bool setMergingPtr(MergingPtr mergingPtrIn) {
    mergingPtr = mergingPtrIn; return true;}

  // Possibility to pass in pointer for merging hooks.
  inline bool setMergingHooksPtr(MergingHooksPtr mergingHooksPtrIn) {
    mergingHooksPtr = mergingHooksPtrIn; return true;}

  // Possibility to pass in pointer for beam shape.
  inline bool setBeamShapePtr(BeamShapePtr beamShapePtrIn) {
    return beamSetup.setBeamShapePtr(beamShapePtrIn);}

  // Possibility to pass in pointer for external cross section,
  // with option to include external phase-space generator.
  inline bool setSigmaPtr(SigmaProcessPtr sigmaPtrIn,
    PhaseSpacePtr phaseSpacePtrIn = nullptr) {
    sigmaPtrs.resize(0), phaseSpacePtrs.resize(0);
    sigmaPtrs.push_back(sigmaPtrIn);
    phaseSpacePtrs.push_back(phaseSpacePtrIn); return true;}

  // Possibility to add further pointers to allow for multiple cross sections.
  inline bool addSigmaPtr(SigmaProcessPtr sigmaPtrIn,
    PhaseSpacePtr phaseSpacePtrIn = nullptr) {
    sigmaPtrs.push_back(sigmaPtrIn);
    phaseSpacePtrs.push_back(phaseSpacePtrIn); return true;}

  // Possibility to insert further pointers to allow for multiple
  // cross sections.
  inline bool insertSigmaPtr(int idx, SigmaProcessPtr sigmaPtrIn,
    PhaseSpacePtr phaseSpacePtrIn = nullptr) {
    if (idx < 0 || idx > (int)sigmaPtrs.size()) return false;
    sigmaPtrs.insert(sigmaPtrs.begin() + idx, sigmaPtrIn);
    phaseSpacePtrs.insert(phaseSpacePtrs.begin() + idx, phaseSpacePtrIn);
    return true;}

  // Possibility to pass in pointer for external resonance.
  inline bool setResonancePtr(ResonanceWidthsPtr resonancePtrIn) {
    resonancePtrs.resize(0);
    resonancePtrs.push_back( resonancePtrIn); return true;}

  // Possibility to add further pointers to allow for multiple resonances.
  inline bool addResonancePtr(ResonanceWidthsPtr resonancePtrIn) {
    resonancePtrs.push_back( resonancePtrIn); return true;}

  // Possibility to insert further pointers to allow for multiple resonances.
  inline bool insertResonancePtr(int idx, ResonanceWidthsPtr resonancePtrIn) {
    if (idx < 0 || idx > (int)resonancePtrs.size()) return false;
    resonancePtrs.insert( resonancePtrs.begin() + idx, resonancePtrIn);
    return true;}

  // Possibility to pass in pointer for external showers.
  inline bool setShowerModelPtr(ShowerModelPtr showerModelPtrIn) {
    showerModelPtr = showerModelPtrIn; return true;}

  // Possibility to pass in pointer for external fragmentation model.
  inline bool setFragmentationPtr(FragmentationModelPtr fragmentationPtrIn) {
    fragPtrs.resize(0);
    fragPtrs.push_back(fragmentationPtrIn); return true;}

  // Possibility to allow for multiple external fragmentation models.
  inline bool addFragmentationPtr(FragmentationModelPtr fragmentationPtrIn) {
    fragPtrs.push_back(fragmentationPtrIn); return true;}

  // Possibility to insert external fragmentation model, in specific position.
  inline bool insertFragmentationPtr(int idx,
    FragmentationModelPtr fragmentationPtrIn) {
    if (idx < 0 || idx > (int)fragPtrs.size()) return false;
    fragPtrs.insert( fragPtrs.begin() + idx, fragmentationPtrIn);
    return true;}

  // Possibility to pass in pointer for modelling of heavy ion collisions.
  inline bool setHeavyIonsPtr(HeavyIonsPtr heavyIonsPtrIn) {
    heavyIonsPtr = heavyIonsPtrIn; return true;}

  // Possibility to pass a HIUserHooks pointer for modifying the
  // behavior of the heavy ion modelling.
  inline bool setHIHooks(HIUserHooksPtr hiHooksPtrIn) {
    hiHooksPtr = hiHooksPtrIn; return true; }

  // Possibility to get the pointer to a object modelling heavy ion
  // collisions.
  inline HeavyIonsPtr getHeavyIonsPtr() { return heavyIonsPtr;}

  // Possibility to access the pointer to the BeamShape object.
  inline BeamShapePtr getBeamShapePtr() { return beamSetup.getBeamShapePtr();}

  // Possibility to get the pointer to the parton-shower model.
  inline ShowerModelPtr getShowerModelPtr() { return showerModelPtr;}

  // Possibility to get the pointer to the LHA accessor.
  inline LHAupPtr getLHAupPtr() { return lhaUpPtr;}

  // Possibility to pass in pointer for setting of parton space-time vertices.
  inline bool setPartonVertexPtr( PartonVertexPtr partonVertexPtrIn) {
    partonVertexPtr = partonVertexPtrIn; return true;}

  // Initialize.
  bool init();

  // Generate the next event.
  inline bool next() { return next(0); }
  bool next(int procTypeIn);

  // Switch to new beam particle identities; for similar hadrons only.
  bool setBeamIDs( int idAin, int idBin = 0);

  // Switch beam kinematics.
  bool setKinematics(double eCMIn);
  bool setKinematics(double eAIn, double eBIn);
  bool setKinematics(double pxAIn, double pyAIn, double pzAIn,
                     double pxBIn, double pyBIn, double pzBIn);
  bool setKinematics(Vec4 pAIn, Vec4 pBIn);

  // Generate only a single timelike shower as in a decay.
  inline int forceTimeShower( int iBeg, int iEnd, double pTmax,
    int nBranchMax = 0) {
    if (!isInit) {
      logger.ERROR_MSG("Pythia is not properly initialized");
      return 0;
    }
    partonSystems.clear();
    infoPrivate.setScalup( 0, pTmax);
    return timesDecPtr->shower(iBeg, iEnd, event, pTmax, nBranchMax);}

  // Generate only the hadronization/decay stage.
  bool forceHadronLevel( bool findJunctions = true);

  // Special routine to allow more decays if on/off switches changed.
  bool moreDecays() {return hadronLevel.moreDecays(event);}
  bool moreDecays(int index) {return hadronLevel.decay(index, event);}

  // Special routine to force R-hadron decay when not done before.
  bool forceRHadronDecays() {return doRHadronDecays();}

  // Do a low-energy collision between two hadrons in the event record.
  inline bool doLowEnergyProcess(int i1, int i2, int procTypeIn) {
    if (!isInit) {
      logger.ERROR_MSG("Pythia is not properly initialized"); return false; }
    return hadronLevel.doLowEnergyProcess( i1, i2, procTypeIn, event); }

  // Get total cross section for two hadrons in the event record or standalone.
  inline double getSigmaTotal() {
    return getSigmaTotal(beamSetup.idA, beamSetup.idB, beamSetup.eCM, 0);}
  inline double getSigmaTotal(double eCM12, int mixLoHi = 0) {
    return getSigmaTotal(beamSetup.idA, beamSetup.idB, eCM12, mixLoHi);}
  inline double getSigmaTotal(int id1, int id2, double eCM12,
    int mixLoHi = 0) {
    return getSigmaTotal(id1, id2, eCM12, particleData.m0(id1),
      particleData.m0(id2), mixLoHi); }
  inline double getSigmaTotal(int id1, int id2, double eCM12, double m1,
    double m2, int mixLoHi = 0) {
    if (!isInit) {
      logger.ERROR_MSG("Pythia is not properly initialized"); return 0.; }
    return sigmaCmb.sigmaTotal(id1, id2, eCM12, m1, m2, mixLoHi); }

  // Get partial (elastic, diffractive, non-diffractive, ...) cross sections
  // for two hadrons in the event record or standalone.
  inline double getSigmaPartial(int procTypeIn) {
    return getSigmaPartial(beamSetup.idA, beamSetup.idB, beamSetup.eCM,
      procTypeIn, 0); }
  inline double getSigmaPartial(double eCM12, int procTypeIn,
    int mixLoHi = 0) {return getSigmaPartial(
      beamSetup.idA, beamSetup.idB, eCM12, procTypeIn, mixLoHi); }
  inline double getSigmaPartial(int id1, int id2, double eCM12, int procTypeIn,
    int mixLoHi = 0) {return getSigmaPartial(id1, id2, eCM12,
      particleData.m0(id1), particleData.m0(id2), procTypeIn, mixLoHi); }
  inline double getSigmaPartial(int id1, int id2, double eCM12, double m1,
    double m2, int procTypeIn, int mixLoHi = 0) {
    if (!isInit) {
      logger.ERROR_MSG("Pythia is not properly initialized"); return 0.;}
    return sigmaCmb.sigmaPartial(
      id1, id2, eCM12, m1, m2, procTypeIn, mixLoHi);}

  // Return a parton density set among list of possibilities.
  inline PDFPtr getPDFPtr(int idIn, int sequence = 1, string beam = "A",
    bool resolved = true) {
    return beamSetup.getPDFPtr( idIn, sequence, beam, resolved); }

  // List the current Les Houches event.
  inline void LHAeventList() { if (lhaUpPtr != 0) lhaUpPtr->listEvent();}

  // Skip a number of Les Houches events at input.
  inline bool LHAeventSkip(int nSkip) {
    if (lhaUpPtr != 0) return lhaUpPtr->skipEvent(nSkip);
    return false;}

  // Main routine to provide final statistics on generation.
  void stat();

  // Read in settings values: shorthand, not new functionality.
  bool   flag(string key) {return settings.flag(key);}
  int    mode(string key) {return settings.mode(key);}
  double parm(string key) {return settings.parm(key);}
  string word(string key) {return settings.word(key);}

  // The event record for the parton-level central process.
  Event           process = {};

  // The event record for the complete event history.
  Event           event = {};

  // Public information and statistic on the generation.
  const Info&     info = infoPrivate;

  // Logger: for diagnostic messages, errors, statistics, etc.
  Logger          logger = {};

  // Settings: databases of flags/modes/parms/words to control run.
  Settings        settings = {};

  // ParticleData: the particle data table/database.
  ParticleData    particleData = {};

  // Random number generator.
  Rndm            rndm = {};

  // Standard Model couplings, including alphaS and alphaEM.
  CoupSM          coupSM = {};

  // SUSY couplings.
  CoupSUSY        coupSUSY = {};

  // SLHA Interface
  SLHAinterface   slhaInterface = {};

  // The partonic content of each subcollision system (auxiliary to event).
  PartonSystems   partonSystems = {};

  // Merging object as wrapper for matrix element merging routines.
  MergingPtr      mergingPtr = {};

  // Pointer to MergingHooks object for user interaction with the merging.
  // MergingHooks also more generally steers the matrix element merging.
  MergingHooksPtr  mergingHooksPtr;

  // Pointer to a HeavyIons object for generating heavy ion collisions
  HeavyIonsPtr   heavyIonsPtr = {};

  // Pointer to a HIUserHooks object to modify heavy ion modelling.
  HIUserHooksPtr hiHooksPtr = {};

  // HadronWidths: the hadron widths data table/database.
  HadronWidths    hadronWidths = {};

  // The two incoming beams.
  const BeamParticle& beamA = beamSetup.beamA;
  const BeamParticle& beamB = beamSetup.beamB;

private:

  // Friend PythiaParallel to give full access to underlying info.
  friend class PythiaParallel;

  // The collector of all event generation weights that should eventually
  // be transferred to the final output.
  WeightContainer weightContainer = {};

  // The main keeper/collector of information, accessible from all
  // PhysicsBase objects. The information is available from the
  // outside through the public info object.
  Info infoPrivate = {};

  // Initialise new Pythia object (called by constructors).
  void initPtrs();

  // Initialise user provided plugins.
  void initPlugins();

  // Functions to be called at the beginning and end of a next() call.
  void beginEvent();
  void endEvent(PhysicsBase::Status status);

  // Register a PhysicsBase object and give it a pointer to the info object.
  inline void registerPhysicsBase(PhysicsBase &pb) {
    if (find(physicsPtrs.begin(), physicsPtrs.end(), &pb) != physicsPtrs.end())
      return;
    pb.initInfoPtr(infoPrivate);
    physicsPtrs.push_back(&pb);
  }

  // If new pointers are set in Info propagate this to all
  // PhysicsBase objects.
  void pushInfo();

  // Constants: could only be changed in the code itself.
  static const double VERSIONNUMBERHEAD, VERSIONNUMBERCODE;
  // Maximum number of tries to produce parton level from given input.
  // Negative integer to denote that no subrun has been set.
  static const int    NTRY = 10;

  // Initialization data, extracted from database.
  string xmlPath = {};
  bool   doProcessLevel = {}, doPartonLevel = {}, doHadronLevel = {},
         doLowEnergy = {}, doSoftQCDall = {}, doSoftQCDinel = {},
         doCentralDiff = {}, doDiffraction = {}, doSoftQCD = {},
         doHardDiff = {}, doResDec = {}, doFSRinRes = {}, decayRHadrons = {},
         doPartonVertex = {}, abortIfVeto = {}, checkEvent = {},
         checkHistory = {}, doNonPert = {};
  int    nErrList = {};
  double epTolErr = {}, epTolWarn = {}, mTolErr = {}, mTolWarn = {};

  // Initialization data, from init(...) call, plus some event-specific.
  bool   isConstructed = {}, isInit = {}, showSaV = {}, showMaD = {},
         doReconnect = {}, forceHadronLevelCR = {};
  int    nCount = {}, nShowLHA = {}, nShowInfo = {}, nShowProc = {},
         nShowEvt = {}, reconnectMode = {};

  // information for error checkout.
  int    nErrEvent = {};
  vector<int> iErrId = {}, iErrCol = {}, iErrEpm = {}, iErrNan = {},
         iErrNanVtx = {};

  // Setup of beams: flavours, kinematics and PDFs.
  BeamSetup beamSetup = {};

  // LHAup object for generating external events.
  bool     doLHA = false;
  bool     useNewLHA = false;
  LHAupPtr lhaUpPtr = {};

  // Pointer to external decay handler and list of particles it handles.
  DecayHandlerPtr decayHandlePtr = {};
  vector<int>     handledParticles = {};

  // Pointer to UserHooks object for user interaction with program.
  UserHooksPtr userHooksPtr = {};
  bool doVetoProcess = {}, doVetoPartons = {},
       retryPartonLevel = {}, canVetoHadronization = {};

  // Pointer to BeamShape object for beam momentum and interaction vertex.
  BeamShapePtr beamShapePtr = {};
  bool         doVertexSpread = {};
  double       eMinPert = {}, eWidthPert = {};

  // Pointers to external processes derived from the Pythia base classes.
  vector<SigmaProcessPtr> sigmaPtrs = {};

  // Pointers to external phase-space generators derived from Pythia
  // base classes.
  vector<PhaseSpacePtr> phaseSpacePtrs = {};

  // Pointers to external calculation of resonance widths.
  vector<ResonanceWidthsPtr> resonancePtrs = {};

  // Pointers to timelike and spacelike showers, including Vincia and
  // Dire. Note, the showerModelPtr must be declared before the
  // individual shower pointers, partonLevel, and hadronLevel. This
  // ensures shared pointers that belong to showerModelPtr are
  // unloaded correctly.
  ShowerModelPtr showerModelPtr = {};
  TimeShowerPtr  timesDecPtr = {};
  TimeShowerPtr  timesPtr = {};
  SpaceShowerPtr spacePtr = {};

  // Pointers to fragmentation models.
  vector<FragmentationModelPtr> fragPtrs = {};

  // Pointer to assign space-time vertices during parton evolution.
  PartonVertexPtr partonVertexPtr;

  // The main generator class to define the core process of the event.
  ProcessLevel processLevel = {};

  // The main generator class to produce the parton level of the event.
  PartonLevel partonLevel = {};

  // The main generator class to perform trial showers of the event.
  PartonLevel trialPartonLevel = {};

  // Flags for defining the merging scheme.
  bool        doMerging = {};

  // The StringInteractions class.
  StringIntPtr stringInteractionsPtr;

  // The Colour reconnection class.
  ColRecPtr colourReconnectionPtr = {};

  // The junction splitting class.
  JunctionSplitting junctionSplitting = {};

  // The main generator class to produce the hadron level of the event.
  HadronLevel hadronLevel = {};

  // The total cross section classes are used both on process and parton level.
  SigmaTotal         sigmaTot = {};
  SigmaLowEnergy     sigmaLowEnergy;
  NucleonExcitations nucleonExcitations = {};
  SigmaCombined      sigmaCmb = {};

  // The fragmentation pointer is used in low energy processes and HadronLevel.
  LundFragmentationPtr fragPtr{};

  // The RHadrons class is used both at PartonLevel and HadronLevel.
  RHadronsPtr rHadronsPtr{};

  // Flags for handling generation of heavy ion collisons.
  bool        hasHeavyIons = {}, doHeavyIons = {};

  // Write the Pythia banner, with symbol and version information.
  void banner();

  // Check that combinations of settings are allowed; change if not.
  void checkSettings();

  // Simplified treatment for low-energy nonperturbative collisions.
  bool nextNonPert(int procTypeIn = 0);

  // Perform R-hadron decays.
  bool doRHadronDecays();

  // Check that the final event makes sense.
  bool check();

  // Initialization of SLHA data.
  bool initSLHA();
  stringstream particleDataBuffer;

  // Keep track of and initialize all pointers to PhysicsBase-derived objects.
  vector<PhysicsBase*> physicsPtrs = {};

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Pythia_H
