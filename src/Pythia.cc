// Pythia.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Pythia class.

#include "Pythia8/Pythia.h"
#include "Pythia8/ColourReconnection.h"
#include "Pythia8/Dire.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8/History.h"
#include "Pythia8/StringInteractions.h"
#include "Pythia8/Vincia.h"
#include "Pythia8/Plugins.h"

namespace Pythia8 {

//==========================================================================

// The Pythia class.

//--------------------------------------------------------------------------

// The current Pythia (sub)version number, to agree with XML version.
const double Pythia::VERSIONNUMBERHEAD = PYTHIA_VERSION;
const double Pythia::VERSIONNUMBERCODE = 8.311;

//--------------------------------------------------------------------------

// Constructor.

Pythia::Pythia(string xmlDir, bool printBanner) {

  // Initialise / reset pointers and global variables.
  initPtrs();

  // Find path to data files, i.e. xmldoc directory location.
  // Environment variable takes precedence, then constructor input,
  // and finally the pre-processor constant XMLDIR.
  const char* envPath = getenv("PYTHIA8DATA");
  xmlPath = envPath ? envPath : "";
  if (xmlPath == "") {
    if (xmlDir.length() && xmlDir[xmlDir.length() - 1] != '/') xmlDir += "/";
    xmlPath = xmlDir;
    ifstream xmlFile((xmlPath + "Index.xml").c_str());
    if (!xmlFile.good()) xmlPath = XMLDIR;
    xmlFile.close();
  }
  if (xmlPath.empty() || (xmlPath.length() && xmlPath[xmlPath.length() - 1]
      != '/')) xmlPath += "/";

  // Read in files with all flags, modes, parms and words.
  settings.initPtrs(&logger);
  string initFile = xmlPath + "Index.xml";
  isConstructed = settings.init( initFile);
  if (!isConstructed) {
    logger.ABORT_MSG("settings unavailable");
    return;
  }

  // Save XML path in settings.
  settings.addWord("xmlPath", xmlPath);

  // Allow include files to be saved to settings.
  settings.addWord("include", "");

  // Check that XML and header version numbers match code version number.
  if (!checkVersion()) return;

  // Read in files with all particle data.
  particleData.initPtrs( &infoPrivate);
  string dataFile = xmlPath + "ParticleData.xml";
  isConstructed = particleData.init( dataFile);
  if (!isConstructed) {
    logger.ABORT_MSG("particle data unavailable");
    return;
  }

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  infoPrivate.addCounter(0);

  // Special settings needed for heavy ion setup.
  HeavyIons::addSpecialSettings(settings);

}

//--------------------------------------------------------------------------

// Constructor from pre-initialised ParticleData and Settings objects.

Pythia::Pythia(Settings& settingsIn, ParticleData& particleDataIn,
  bool printBanner) {

  // Initialise / reset pointers and global variables.
  initPtrs();

  // Copy XML path from existing Settings database.
  xmlPath = settingsIn.word("xmlPath");

  // Copy settings database and redirect pointers.
  settings = settingsIn;
  settings.initPtrs(&logger);
  isConstructed = settings.getIsInit();
  if (!isConstructed) {
    logger.ABORT_MSG("settings unavailable");
    return;
  }

  // Check XML and header version numbers match code version number.
  if (!checkVersion()) return;

  // Copy particleData database and redirect pointers.
  particleData = particleDataIn;
  particleData.initPtrs( &infoPrivate);
  isConstructed = particleData.getIsInit();
  if (!isConstructed) {
    logger.ABORT_MSG("particle data unavailable");
    return;
  }

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  infoPrivate.addCounter(0);

}

//--------------------------------------------------------------------------

// Constructor from string streams.

Pythia::Pythia( istream& settingsStrings, istream& particleDataStrings,
  bool printBanner) {

  // Initialise / reset pointers and global variables.
  initPtrs();

  // Copy settings database
  settings.init( settingsStrings );

  // Reset pointers to pertain to this PYTHIA object.
  settings.initPtrs( &logger);
  isConstructed = settings.getIsInit();
  if (!isConstructed) {
    logger.ABORT_MSG("settings unavailable");
    return;
  }

  // Check XML and header version numbers match code version number.
  if (!checkVersion()) return;

  // Read in files with all particle data.
  particleData.initPtrs( &infoPrivate);
  isConstructed = particleData.init( particleDataStrings );
  if (!isConstructed) {
    logger.ABORT_MSG("particle data unavailable");
    return;
  }

  // Write the Pythia banner to output.
  if (printBanner) banner();

  // Not initialized until at the end of the init() call.
  isInit = false;
  infoPrivate.addCounter(0);

}

//--------------------------------------------------------------------------

// Initialise new Pythia object (common code called by constructors).

void Pythia::initPtrs() {

  // Setup of Info.
  infoPrivate.setPtrs( &settings, &particleData, &logger, &rndm, &beamSetup,
    &coupSM, &coupSUSY, &partonSystems, &sigmaTot, &sigmaCmb,
    &hadronWidths, &weightContainer);
  registerPhysicsBase(processLevel);
  registerPhysicsBase(partonLevel);
  registerPhysicsBase(trialPartonLevel);
  registerPhysicsBase(hadronLevel);
  registerPhysicsBase(sigmaTot);
  registerPhysicsBase(nucleonExcitations);
  registerPhysicsBase(sigmaLowEnergy);
  registerPhysicsBase(sigmaCmb);
  registerPhysicsBase(hadronWidths);
  registerPhysicsBase(junctionSplitting);
  registerPhysicsBase(rHadrons);
  registerPhysicsBase(beamSetup);

}

//--------------------------------------------------------------------------

// Initialise user provided plugins.

void Pythia::initPlugins() {

  // Create map of pointers to pass to Pythia.
  map<string, PDFPtr> pdfs = beamSetup.getPDFPtr();

  // Store the previous plugin type.
  string objLast;

  // Loop over the plugins.
  for (string plugin : settings.wvec("Init:plugins")) {

    // Split by the plugin entry by the "::" delimiter.
    vector<string> vals;
    size_t pos(0);
    string cmnd(plugin);
    while (pos != string::npos) {
      pos = plugin.find("::");
      vals.push_back(plugin.substr(0, pos));
      plugin = plugin.substr(pos + 2);
    }
    if (vals.size() < 2) {
      logger.ERROR_MSG("plugin setting " + plugin +
        " must at least specify library and class");
      continue;
    }
    string libName(vals[0]), className(vals[1]), key("default"), fileName("");
    int sub(SUBRUNDEFAULT);
    if (vals.size() > 2) key = vals[2];
    if (vals.size() > 3) fileName = vals[3];
    if (vals.size() > 4) sub = stoi(vals[4]);

    // Determine the plugin type.
    shared_ptr<void> libPtr = dlopen_plugin(libName, &logger);
    if (libPtr == nullptr) continue;
    string objType = type_plugin(libName, className, &logger);
    if (objType == "") continue;

    // Handle the case of a PDF plugin.
    if (objType == typeid(PDF).name()) {
      // If default beam type, set both A and B.
      if (key == "default") {
        for (string keyNow : {"A", "B"}) {
          if (pdfs[keyNow] != nullptr)
            logger.WARNING_MSG("replacing PDF pointer " + keyNow);
          pdfs[keyNow] = make_plugin<PDF>(
            libName, className, this, fileName, sub);
        }
        continue;
      }
      // Otherwise, set the specified PDF beam.
      if (pdfs.find(key) == pdfs.end()) {
        logger.ERROR_MSG("cannot set unknown PDF pointer " + key);
        continue;
      } else if (pdfs[key] != nullptr)
        logger.WARNING_MSG("replacing PDF pointer " + key);
      pdfs[key] = make_plugin<PDF>(
        libName, className, this, fileName, sub);
      continue;
    }

    // Check if the plugin should be set or added.
    if (key == "default") key = "set";
    if (key != "set" && key != "add") {
      logger.ERROR_MSG("the third argument must be set or add");
      continue;
    }

    // Load the other plugin types and set appropriately.
    if (vals.size() > 3) fileName = vals[3];
    if (vals.size() > 4) sub = stoi(vals[4]);
    if (objType == typeid(LHAup).name()) {setLHAupPtr(
        make_plugin<LHAup>(libName, className, this, fileName, sub));
    } else if (objType == typeid(DecayHandler).name()) {setDecayPtr(
        make_plugin<DecayHandler>(libName, className, this, fileName, sub));
    } else if (objType == typeid(RndmEngine).name()) {setRndmEnginePtr(
        make_plugin<RndmEngine>(libName, className, this, fileName, sub));
    } else if (objType == typeid(UserHooks).name()) {
      if (key == "add") addUserHooksPtr(
        make_plugin<UserHooks>(libName, className, this, fileName, sub));
      else setUserHooksPtr(
        make_plugin<UserHooks>(libName, className, this, fileName, sub));
    } else if (objType == typeid(Merging).name()) {setMergingPtr(
        make_plugin<Merging>(libName, className, this, fileName, sub));
    } else if (objType == typeid(MergingHooks).name()) {setMergingHooksPtr(
        make_plugin<MergingHooks>(libName, className, this, fileName, sub));
    } else if (objType == typeid(BeamShape).name()) {setBeamShapePtr(
        make_plugin<BeamShape>(libName, className, this, fileName, sub));
    } else if (objType == typeid(SigmaProcess).name()) {
      if (key == "add") addSigmaPtr(
        make_plugin<SigmaProcess>(libName, className, this, fileName, sub));
      else setSigmaPtr(
        make_plugin<SigmaProcess>(libName, className, this, fileName, sub));
    } else if (objType == typeid(PhaseSpace).name() &&
      objLast == typeid(SigmaProcess).name()) {
      SigmaProcessPtr sigmaPtr = sigmaPtrs.back();
      if (key == "add") {sigmaPtrs.pop_back(); phaseSpacePtrs.pop_back();}
      else {sigmaPtrs.resize(0); phaseSpacePtrs.resize(0);}
      addSigmaPtr(sigmaPtr,
        make_plugin<PhaseSpace>(libName, className, this, fileName, sub));
    } else if (objType == typeid(ResonanceWidths).name()) {
      if (key == "add") addResonancePtr(
       make_plugin<ResonanceWidths>(libName, className, this, fileName, sub));
      else setResonancePtr(
       make_plugin<ResonanceWidths>(libName, className, this, fileName, sub));
    } else if (objType == typeid(ShowerModel).name()) {setShowerModelPtr(
        make_plugin<ShowerModel>(libName, className, this, fileName, sub));
    } else if (objType == typeid(HeavyIons).name()) {setHeavyIonsPtr(
        make_plugin<HeavyIons>(libName, className, this, fileName, sub));
    } else if (objType == typeid(HIUserHooks).name()) {setHIHooks(
        make_plugin<HIUserHooks>(libName, className, this, fileName, sub));
    } else {logger.ERROR_MSG("the class " + demangle(objType) + " cannot be "
        "passed to a Pythia object");
    }
    objLast = objType;
  }

  // Set the PDF pointers.
  setPDFPtr(pdfs["A"], pdfs["B"], pdfs["HardA"], pdfs["HardB"],
    pdfs["PomA"], pdfs["PomB"], pdfs["GamA"], pdfs["GamB"],
    pdfs["HardGamA"], pdfs["HardGamB"], pdfs["UnresA"], pdfs["UnresB"],
    pdfs["UnresGamA"], pdfs["UnresGamB"], pdfs["VMDA"], pdfs["VMDB"]);

}

//--------------------------------------------------------------------------

// Check for consistency of version numbers (called by constructors).

bool Pythia::checkVersion() {

  // Check that XML version number matches code version number.
  double versionNumberXML = parm("Pythia:versionNumber");
  isConstructed = (abs(versionNumberXML - VERSIONNUMBERCODE) < 0.0005);
  if (!isConstructed) {
    ostringstream errCode;
    errCode << fixed << setprecision(3) << ": in code " << VERSIONNUMBERCODE
            << " but in XML " << versionNumberXML;
    logger.ABORT_MSG("unmatched version numbers", errCode.str());
    return false;
  }

  // Check that header version number matches code version number.
  isConstructed = (abs(VERSIONNUMBERHEAD - VERSIONNUMBERCODE) < 0.0005);
  if (!isConstructed) {
    ostringstream errCode;
    errCode << fixed << setprecision(3) << ": in code " << VERSIONNUMBERCODE
            << " but in header " << VERSIONNUMBERHEAD;
    logger.ABORT_MSG("unmatched version numbers",  errCode.str());
    return false;
  }

  // All is well that ends well.
  return true;

}

//--------------------------------------------------------------------------

// Read in one update for a setting or particle data from a single line.

bool Pythia::readString(string line, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // If empty line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // If Settings input stretches over several lines then continue with it.
  if (settings.unfinishedInput()) return settings.readString(line, warn);

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalnum(line[firstChar])) return true;

  // Send on particle data to the ParticleData database.
  if (isdigit(line[firstChar])) {
    bool passed = particleData.readString(line, warn);
    if (passed) particleDataBuffer << line << endl;
    return passed;
  }

  // Include statements.
  if (line.find("include") == 0 && settings.readString(line, warn) &&
    word("include") != "") {

    // Try normal path first.
    string fileName = word("include");
    settings.word("include", "");
    ifstream isUser(fileName.c_str());
    if (!isUser.good()) {

      // Split the paths from PYTHIA8CMND.
      vector<string> paths;
      size_t pos(0);
      const char* envChar = getenv("PYTHIA8CMND");
      string envPath = envChar ? envChar : "";
      while (envPath != "" && pos != string::npos) {
        pos = envPath.find(":");
        paths.push_back(envPath.substr(0, pos));
        envPath = envPath.substr(pos + 1);
      }

      // Add the Pythia settings directory.
      paths.push_back(word("xmlPath"). substr(0, xmlPath.length() - 7)
        + "settings");

      // Try the different paths.
      for (string path : paths) {
        ifstream isPath((path + "/" + fileName).c_str());
        if (isPath.good()) return readFile(isPath, warn, subrun);
      }
      logger.ERROR_MSG("did not find file", fileName);
      return false;
    } else return readFile(isUser, warn, subrun);
  }

  // Everything else sent on to Settings.
  return settings.readString(line, warn);

}

//--------------------------------------------------------------------------

// Read in updates for settings or particle data from user-defined file.

bool Pythia::readFile(string fileName, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Open file for reading.
  ifstream is(fileName.c_str());
  if (!is.good()) {
    logger.ERROR_MSG("did not find file", fileName);
    return false;
  }

  // Hand over real work to next method.
  return readFile( is, warn, subrun);

}

//--------------------------------------------------------------------------

// Read in updates for settings or particle data
// from user-defined stream (or file).

bool Pythia::readFile(istream& is, bool warn, int subrun) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Read in one line at a time.
  string line;
  bool isCommented = false;
  bool accepted = true;
  int subrunNow = SUBRUNDEFAULT;
  while ( getline(is, line) ) {

    // Check whether entering, leaving or inside commented-commands section.
    int commentLine = readCommented( line);
    if      (commentLine == +1)  isCommented = true;
    else if (commentLine == -1)  isCommented = false;
    else if (isCommented) ;

    else {
      // Check whether entered new subrun.
      int subrunLine = readSubrun( line, warn);
      if (subrunLine >= 0) subrunNow = subrunLine;

      // Process the line if in correct subrun.
      if ( (subrunNow == subrun || subrunNow == SUBRUNDEFAULT)
         && !readString( line, warn) ) accepted = false;
    }

  // Reached end of input file.
  };
  return accepted;

}

//--------------------------------------------------------------------------

// Routine to initialize with the variable values of the Beams kind.

bool Pythia::init() {

  // Check that constructor worked.
  if (!isConstructed) {
    logger.ABORT_MSG("constructor initialization failed");
    isInit = false;
    return false;
  }

  // Check if this is the first call to init.
  if (isInit)
    logger.WARNING_MSG("be aware that successive "
      "calls to init() do not clear previous settings");
  // Only create plugins the first time.
  else initPlugins();
  isInit = false;

  // Early catching of heavy ion mode.
  doHeavyIons = HeavyIons::isHeavyIon(settings) || mode("HeavyIon:mode") == 2;
  if ( doHeavyIons ) {
    if ( !heavyIonsPtr ) heavyIonsPtr = make_shared<Angantyr>(*this);
    registerPhysicsBase(*heavyIonsPtr);
    if ( hiHooksPtr ) heavyIonsPtr->setHIUserHooksPtr(hiHooksPtr);
    if ( !heavyIonsPtr->init() ) {
      doHeavyIons = false;
      isInit = false;
      logger.ABORT_MSG("the heavy ion model "
      "failed to initialize. Double check settings");
      return false;
    }
  }

  // Early readout, if return false or changed when no beams.
  doProcessLevel = flag("ProcessLevel:all");

  // Check that changes in Settings and ParticleData have worked.
  if (settings.unfinishedInput()) {
    logger.ABORT_MSG("opening { not matched by closing }");
    return false;
  }
  if (settings.readingFailed()) {
    logger.ABORT_MSG("some user settings did not make sense");
    return false;
  }
  if (particleData.readingFailed()) {
    logger.ABORT_MSG("some user particle data did not make sense");
    return false;
  }

  // Initialize error printing settings.
  logger.init(settings);

  // Initialize the random number generator.
  if ( flag("Random:setSeed") ) rndm.init( mode("Random:seed") );
  else                          rndm.init(Rndm::DEFAULTSEED);

  // Count up number of initializations.
  infoPrivate.addCounter(1);

  // Master choice of shower model.
  int showerModel = mode("PartonShowers:model");

  // Set up values related to CKKW-L merging.
  bool doUserMerging     = flag("Merging:doUserMerging");
  bool doMGMerging       = flag("Merging:doMGMerging");
  bool doKTMerging       = flag("Merging:doKTMerging");
  bool doPTLundMerging   = flag("Merging:doPTLundMerging");
  bool doCutBasedMerging = flag("Merging:doCutBasedMerging");
  // Set up values related to unitarised CKKW merging
  bool doUMEPSTree       = flag("Merging:doUMEPSTree");
  bool doUMEPSSubt       = flag("Merging:doUMEPSSubt");
  // Set up values related to NL3 NLO merging
  bool doNL3Tree         = flag("Merging:doNL3Tree");
  bool doNL3Loop         = flag("Merging:doNL3Loop");
  bool doNL3Subt         = flag("Merging:doNL3Subt");
  // Set up values related to unitarised NLO merging
  bool doUNLOPSTree      = flag("Merging:doUNLOPSTree");
  bool doUNLOPSLoop      = flag("Merging:doUNLOPSLoop");
  bool doUNLOPSSubt      = flag("Merging:doUNLOPSSubt");
  bool doUNLOPSSubtNLO   = flag("Merging:doUNLOPSSubtNLO");
  bool doXSectionEst     = flag("Merging:doXSectionEstimate");
  doMerging = doUserMerging || doMGMerging || doKTMerging
    || doPTLundMerging || doCutBasedMerging || doUMEPSTree || doUMEPSSubt
    || doNL3Tree || doNL3Loop || doNL3Subt || doUNLOPSTree
    || doUNLOPSLoop || doUNLOPSSubt || doUNLOPSSubtNLO || doXSectionEst;
  doMerging = doMerging || flag("Merging:doMerging");

  // If not using Vincia, Merging:Process must not be enclosed in {}.
  if (doMerging && showerModel != 2) {
    string mergingProc = word("Merging:Process");
    if (mergingProc.front()=='{' && mergingProc.back() == '}') {
      stringstream ss;
      ss<<"invalid Merging:Process for PartonShower:Model = "<<showerModel;
      logger.ABORT_MSG(ss.str());
      return false;
    }
  }

  // Set/Reset the weights
  weightContainer.initPtrs(&infoPrivate);
  weightContainer.init(doMerging);

  // Set up MergingHooks object for simple shower model.
  if (doMerging && showerModel == 1) {
    if (!mergingHooksPtr) mergingHooksPtr = make_shared<MergingHooks>();
    registerPhysicsBase(*mergingHooksPtr);
    mergingHooksPtr->setLHEInputFile("");
  }

  // Set up Merging object for simple shower model.
  if ( doMerging  && showerModel == 1 && !mergingPtr ) {
    mergingPtr = make_shared<Merging>();
    registerPhysicsBase(*mergingPtr);
  }

  // Set up beams and kinematics.
  if (!beamSetup.initFrame()) return false;

  // Spread LHA information generated in initFrame.
  doLHA     = beamSetup.doLHA;
  useNewLHA = beamSetup.useNewLHA;
  lhaUpPtr  = beamSetup.lhaUpPtr;
  if (doLHA) processLevel.setLHAPtr( lhaUpPtr);

  // Done if only new LHA file.
  if (beamSetup.skipInit) {
    isInit = true;
    infoPrivate.addCounter(2);
    return true;
  }

  // Store the name of the input LHEF for merging.
  if ( (beamSetup.frameType == 4 || beamSetup.frameType == 5)
    && doMerging && showerModel == 1) {
    string lhef = (beamSetup.frameType == 4) ? word("Beams:LHEF") : "";
    mergingHooksPtr->setLHEInputFile( lhef);
  }

  // Set up values related to user hooks.
  doVetoProcess        = false;
  doVetoPartons        = false;
  retryPartonLevel     = false;
  canVetoHadronization = false;

  if ( userHooksPtr ) {
    infoPrivate.userHooksPtr = userHooksPtr;
    registerPhysicsBase(*userHooksPtr);
    pushInfo();
    if (!userHooksPtr->initAfterBeams()) {
      logger.ABORT_MSG("could not initialise UserHooks");
      return false;
    }
    doVetoProcess       = userHooksPtr->canVetoProcessLevel();
    doVetoPartons       = userHooksPtr->canVetoPartonLevel();
    retryPartonLevel    = userHooksPtr->retryPartonLevel();
    canVetoHadronization = userHooksPtr->canVetoAfterHadronization();
  }

  // Setup objects for string interactions (swing, colour
  // reconnection, shoving and rope hadronisation).
  if ( !stringInteractionsPtr ) {
    if ( flag("Ropewalk:RopeHadronization") )
      stringInteractionsPtr = make_shared<Ropewalk>();
    else
      stringInteractionsPtr = make_shared<StringInteractions>();
    registerPhysicsBase(*stringInteractionsPtr);

  }
  stringInteractionsPtr->init();

  // Back to common initialization. Reset error counters.
  nErrEvent = 0;
  infoPrivate.setTooLowPTmin(false);
  infoPrivate.sigmaReset();

  // Initialize data members extracted from database.
  doPartonLevel    = flag("PartonLevel:all");
  doHadronLevel    = flag("HadronLevel:all");
  doLowEnergy      = flag("LowEnergyQCD:all")
                  || flag("LowEnergyQCD:nonDiffractive")
                  || flag("LowEnergyQCD:elastic")
                  || flag("LowEnergyQCD:singleDiffractiveXB")
                  || flag("LowEnergyQCD:singleDiffractiveAX")
                  || flag("LowEnergyQCD:doubleDiffractive")
                  || flag("LowEnergyQCD:excitation")
                  || flag("LowEnergyQCD:annihilation")
                  || flag("LowEnergyQCD:resonant");
  doSoftQCDall     = flag("SoftQCD:all");
  doSoftQCDinel    = flag("SoftQCD:inelastic");
  doCentralDiff    = flag("SoftQCD:centralDiffractive");
  doDiffraction    = flag("SoftQCD:singleDiffractive")
                  || flag("SoftQCD:singleDiffractiveXB")
                  || flag("SoftQCD:singleDiffractiveAX")
                  || flag("SoftQCD:doubleDiffractive")
                  || doSoftQCDall || doSoftQCDinel || doCentralDiff;
  doSoftQCD        = doDiffraction ||
                     flag("SoftQCD:nonDiffractive") ||
                     flag("SoftQCD:elastic");
  doHardDiff       = flag("Diffraction:doHard");
  doResDec         = flag("ProcessLevel:resonanceDecays");
  doFSRinRes       = doPartonLevel && flag("PartonLevel:FSR")
                  && flag("PartonLevel:FSRinResonances");
  decayRHadrons    = flag("RHadrons:allowDecay");
  doPartonVertex   = flag("PartonVertex:setVertex");
  eMinPert         = parm("Beams:eMinPert");
  eWidthPert       = parm("Beams:eWidthPert");
  abortIfVeto      = flag("Check:abortIfVeto");
  checkEvent       = flag("Check:event");
  checkHistory     = flag("Check:history");
  nErrList         = mode("Check:nErrList");
  epTolErr         = parm("Check:epTolErr");
  epTolWarn        = parm("Check:epTolWarn");
  mTolErr          = parm("Check:mTolErr");
  mTolWarn         = parm("Check:mTolWarn");

  // Warn/abort for unallowed process and beam combinations.
  bool doHardProc  = settings.hasHardProc() || doLHA;
  if (doSoftQCD && doHardProc) {
    logger.WARNING_MSG("should not combine softQCD processes with hard ones");
  }
  if (beamSetup.doVarEcm && doHardProc) {
    logger.ABORT_MSG("variable energy only works for softQCD processes");
    return false;
  }
  if (doLowEnergy && doHardProc) {
    logger.ABORT_MSG("lowEnergy and hard processes cannot be used together");
    return false;
  }

  // Check that combinations of settings are allowed; change if not.
  checkSettings();

  // Initialize the SM couplings (needed to initialize resonances).
  coupSM.init( settings, &rndm );

  // Initialize SLHA interface (including SLHA/BSM couplings).
  bool useSLHAcouplings = false;
  slhaInterface = SLHAinterface();
  slhaInterface.setPtr( &infoPrivate);
  slhaInterface.init( useSLHAcouplings, particleDataBuffer );

  // Reset couplingsPtr to the correct memory address.
  particleData.initPtrs( &infoPrivate);

  // Set headers to distinguish the two event listing kinds.
  int startColTag = mode("Event:startColTag");
  process.init("(hard process)", &particleData, startColTag);
  event.init("(complete event)", &particleData, startColTag);

  // Final setup stage of particle data, notably resonance widths.
  particleData.initWidths( resonancePtrs);

  // Read in files with particle widths.
  string dataFile = xmlPath + "HadronWidths.dat";
  if (!hadronWidths.init(dataFile)) {
    logger.ABORT_MSG("hadron widths unavailable");
    return false;
  }
  if (!hadronWidths.check()) {
    logger.ABORT_MSG("hadron widths are invalid");
    return false;
  }

  // Set up R-hadrons particle data, where relevant.
  rHadrons.init();

  // Set up and initialize setting of parton production vertices.
  if (doPartonVertex) {
    if (!partonVertexPtr) partonVertexPtr = make_shared<PartonVertex>();
    registerPhysicsBase(*partonVertexPtr);
    partonVertexPtr->init();
  }

  // Prepare for variable-beam and -energy cross sections.
  string dataFileNucl = xmlPath + "NucleonExcitations.dat";
  if (!nucleonExcitations.init(dataFileNucl)) {
    logger.ABORT_MSG("nucleon excitation data unavailable");
    return false;
  }
  sigmaLowEnergy.init( &nucleonExcitations);
  sigmaCmb.init( &sigmaLowEnergy);

  // Prepare for low-energy QCD processes.
  doNonPert = hadronLevel.initLowEnergyProcesses();
  if (doNonPert && !doSoftQCD && !doHardProc) doProcessLevel = false;

  // Initialise merging hooks for simple shower model.
  if ( doMerging && mergingHooksPtr && showerModel == 1 ) {
    mergingHooksPtr->init();
  }

  // Set up and initialize the ShowerModel (if not provided by user).
  if ( !showerModelPtr ) {
    if ( showerModel == 2 ) showerModelPtr = make_shared<Vincia>();
    else if (showerModel == 3 ) showerModelPtr = make_shared<Dire>();
    else showerModelPtr = make_shared<SimpleShowerModel>();
  }

  // Register shower model as physicsBase object (also sets pointers)
  registerPhysicsBase(*showerModelPtr);

  // Initialise shower model
  if ( !showerModelPtr->init(mergingPtr, mergingHooksPtr,
    partonVertexPtr, &weightContainer) ) {
    logger.ABORT_MSG("shower model failed to initialise");
    return false;
  }

  // Vincia adds a userhook -> need to push this to all physics objects.
  if ( showerModel == 2 ) pushInfo();

  // Get relevant pointers from shower models.
  if (doMerging && showerModel != 1) {
    mergingPtr  = showerModelPtr->getMerging();
    mergingHooksPtr = showerModelPtr->getMergingHooks();
  }
  timesPtr    = showerModelPtr->getTimeShower();
  timesDecPtr = showerModelPtr->getTimeDecShower();
  spacePtr    = showerModelPtr->getSpaceShower();

  // Initialize pointers in showers.
  if (showerModel == 1 || showerModel == 3) {
    if ( timesPtr )
      timesPtr->initPtrs( mergingHooksPtr, partonVertexPtr,
        &weightContainer);
    if ( timesDecPtr && timesDecPtr != timesPtr )
      timesDecPtr->initPtrs( mergingHooksPtr, partonVertexPtr,
        &weightContainer);
    if ( spacePtr )
      spacePtr->initPtrs( mergingHooksPtr, partonVertexPtr,
        &weightContainer);
  }

  // At this point, the mergingPtr should be set. If that is not the
  // case, then the initialization should be aborted.
  if (doMerging && !mergingPtr) {
    logger.ABORT_MSG(
      "merging requested, but merging class not correctly created");
    return false;
  }

  // Set up the beams.
  StringFlav* flavSelPtr = hadronLevel.getStringFlavPtr();
  if (!beamSetup.initBeams(doNonPert, flavSelPtr)) return false;

  // Spread information on beam switching from beamSetup.
  if ( beamSetup.allowIDAswitch)
    partonLevel.initSwitchID( beamSetup.idAList);

  // Turn off central diffraction for VMD processes.
  if (beamSetup.getVMDsideA() || beamSetup.getVMDsideB()) {
    if (doCentralDiff) {
      logger.WARNING_MSG(
        "central diffractive events not implemented for gamma + p/gamma");
      return false;
    }
    if (doSoftQCDall) {
      logger.WARNING_MSG(
        "central diffractive events not implemented for gamma + p/gamma");
      settings.flag("SoftQCD:all", false);
      settings.flag("SoftQCD:elastic", true);
      settings.flag("SoftQCD:nonDiffractive", true);
      settings.flag("SoftQCD:singleDiffractive", true);
      settings.flag("SoftQCD:doubleDiffractive", true);
    }
    if (doSoftQCDinel) {
      logger.WARNING_MSG(
        "central diffractive events not implemented for gamma + p/gamma");
      settings.flag("SoftQCD:inelastic", false);
      settings.flag("SoftQCD:nonDiffractive", true);
      settings.flag("SoftQCD:singleDiffractive", true);
      settings.flag("SoftQCD:doubleDiffractive", true);
    }
  }

  // Send info/pointers to process level for initialization.
  if ( doProcessLevel ) {
    sigmaTot.init();
    if (!processLevel.init(doLHA, &slhaInterface, sigmaPtrs, phaseSpacePtrs)) {
      logger.ABORT_MSG("processLevel initialization failed");
      return false;
    }
  }

  // Initialize shower handlers using initialized beams.
  if (!showerModelPtr->initAfterBeams()) {
    string message = "Fail to initialize ";
    if      (showerModel==1) message += "default";
    else if (showerModel==2) message += "Vincia";
    else if (showerModel==3) message += "Dire";
    message += " shower.";
    logger.ABORT_MSG(message);
    return false;
  }

  // Initialize timelike showers already here, since needed in decays.
  // The pointers to the beams are needed by some external plugin showers.
  timesDecPtr->init( &beamSetup.beamA, &beamSetup.beamB);

  // Alternatively only initialize resonance decays.
  if ( !doProcessLevel) processLevel.initDecays(lhaUpPtr);

  // Send info/pointers to parton level for initialization.
  if ( doPartonLevel && doProcessLevel && !partonLevel.init(timesDecPtr,
    timesPtr, spacePtr, &rHadrons, mergingHooksPtr,
    partonVertexPtr, stringInteractionsPtr, false) ) {
    logger.ABORT_MSG("partonLevel initialization failed");
    return false;
  }

  // Make pointer to shower available for merging machinery.
  if ( doMerging && mergingHooksPtr )
    mergingHooksPtr->setShowerPointer(&partonLevel);

  // Alternatively only initialize final-state showers in resonance decays.
  if ( (!doProcessLevel || !doPartonLevel)
    && (!doNonPert || doSoftQCD) ) partonLevel.init(
    timesDecPtr, nullptr, nullptr, &rHadrons, nullptr,
    partonVertexPtr, stringInteractionsPtr, false);

  // Set up shower variation groups if enabled
  bool doVariations = flag("UncertaintyBands:doVariations");
  if ( doVariations && showerModel==1 ) weightContainer.weightsShowerPtr->
    initWeightGroups(true);

  // Send info/pointers to parton level for trial shower initialization.
  if ( doMerging && !trialPartonLevel.init( timesDecPtr, timesPtr,
    spacePtr, &rHadrons, mergingHooksPtr, partonVertexPtr,
    stringInteractionsPtr, true) ) {
    logger.ABORT_MSG("trialPartonLevel initialization failed");
    return false;
  }

  // Initialise the merging wrapper class.
  if (doMerging && mergingPtr && mergingHooksPtr) {
    mergingPtr->initPtrs( mergingHooksPtr, &trialPartonLevel);
    mergingPtr->init();
  }

  // Send info/pointers to hadron level for initialization.
  // Note: forceHadronLevel() can come, so we must always initialize.
  if ( !hadronLevel.init( timesDecPtr, &rHadrons, decayHandlePtr,
    handledParticles, stringInteractionsPtr, partonVertexPtr,
    sigmaLowEnergy, nucleonExcitations) ) {
    logger.ABORT_MSG("hadronLevel initialization failed");
    return false;
  }

  // Optionally check particle data table for inconsistencies.
  if ( flag("Check:particleData") )
    particleData.checkTable( mode("Check:levelParticleData") );

  // Optionally show settings and particle data, changed or all.
  bool showCS  = flag("Init:showChangedSettings");
  bool showAS  = flag("Init:showAllSettings");
  bool showCPD = flag("Init:showChangedParticleData");
  bool showCRD = flag("Init:showChangedResonanceData");
  bool showAPD = flag("Init:showAllParticleData");
  int  show1PD = mode("Init:showOneParticleData");
  bool showPro = flag("Init:showProcesses");
  if (showCS)      settings.listChanged();
  if (showAS)      settings.listAll();
  if (show1PD > 0) particleData.list(show1PD);
  if (showCPD)     particleData.listChanged(showCRD);
  if (showAPD)     particleData.listAll();

  // Listing options for the next() routine.
  nCount       = mode("Next:numberCount");
  nShowLHA     = mode("Next:numberShowLHA");
  nShowInfo    = mode("Next:numberShowInfo");
  nShowProc    = mode("Next:numberShowProcess");
  nShowEvt     = mode("Next:numberShowEvent");
  showSaV      = flag("Next:showScaleAndVertex");
  showMaD      = flag("Next:showMothersAndDaughters");

  // Init junction splitting.
  junctionSplitting.init();

  // Flags for colour reconnection.
  doReconnect        = flag("ColourReconnection:reconnect");
  reconnectMode      = mode("ColourReconnection:mode");
  forceHadronLevelCR = flag("ColourReconnection:forceHadronLevelCR");
  if ( doReconnect ) colourReconnectionPtr =
    stringInteractionsPtr->getColourReconnections();

  // Succeeded.
  isInit = true;
  infoPrivate.addCounter(2);
  if (useNewLHA && showPro) lhaUpPtr->listInit();
  return true;

}

//--------------------------------------------------------------------------

// Check that combinations of settings are allowed; change if not.

void Pythia::checkSettings() {

  // Double rescattering not allowed if ISR or FSR.
  if ((flag("PartonLevel:ISR") || flag("PartonLevel:FSR"))
    && flag("MultipartonInteractions:allowDoubleRescatter")) {
    logger.WARNING_MSG(
      "double rescattering switched off since showering is on");
    settings.flag("MultipartonInteractions:allowDoubleRescatter", false);
  }

  // Optimize settings for collisions with direct photon(s).
  if ( beamSetup.beamA2gamma || beamSetup.beamB2gamma
    || (beamSetup.idA == 22) || (beamSetup.idB == 22) ) {
    if ( flag("PartonLevel:MPI") && (beamSetup.gammaMode > 1) ) {
      logger.WARNING_MSG(
        "MPIs turned off for collision with unresolved photon");
      settings.flag("PartonLevel:MPI", false);
    }
    if ( flag("SoftQCD:nonDiffractive")
      && (beamSetup.gammaMode > 1) ) {
      logger.WARNING_MSG(
        "soft QCD processes turned off for collision with unresolved photon");
      settings.flag("SoftQCD:nonDiffractive", false);
    }
  }

}

//--------------------------------------------------------------------------

// Main routine to generate the next event, using internal machinery.

bool Pythia::next(int procType) {

  // Check that constructor worked.
  if (!isConstructed) {
    endEvent(PhysicsBase::CONSTRUCTOR_FAILED);
    return false;
  }

  // Check that initialization worked.
  if (!isInit) {
    logger.ABORT_MSG("not properly initialized so cannot generate events");
    endEvent(PhysicsBase::INIT_FAILED);
    return false;
  }

  // Flexible-use call at the beginning of each new event.
  beginEvent();

  // Check if the generation is taken over by the HeavyIons object.
  // Allows HeavyIons::next to call next for this Pythia object
  // without going into a loop.
  if ( doHeavyIons ) {
    doHeavyIons = false;
    bool ok = heavyIonsPtr->next();
    doHeavyIons = true;
    endEvent(ok ? PhysicsBase::COMPLETE : PhysicsBase::HEAVYION_FAILED);
    return ok;
  }

  // Regularly print how many events have been generated.
  int nPrevious = infoPrivate.getCounter(3);
  if (nCount > 0 && nPrevious > 0 && nPrevious%nCount == 0)
    cout << "\n Pythia::next(): " << nPrevious
         << " events have been generated " << endl;

  // Set/reset info counters specific to each event. Also procType.
  infoPrivate.addCounter(3);
  for (int i = 10; i < 13; ++i) infoPrivate.setCounter(i);
  if (!beamSetup.doVarEcm) procType = 0;

  // Simpler option when no hard process, i.e. mainly hadron level.
  if (!doProcessLevel && !doNonPert) {
    // Optionally fetch in resonance decays from LHA interface.
    if (doLHA && !processLevel.nextLHAdec( event)) {
      if (infoPrivate.atEndOfFile())
        logger.ABORT_MSG("reached end of Les Houches Events File");
      endEvent(PhysicsBase::LHEF_END);
      return false;
    }

    // Reset info and partonSystems arrays (while event record contains data).
    infoPrivate.clear();
    weightContainer.clear();
    partonSystems.clear();

    // Set correct energy for system.
    Vec4 pSum = 0.;
    for (int i = 1; i < event.size(); ++i)
      if (event[i].isFinal()) pSum += event[i].p();
    event[0].p( pSum );
    event[0].m( pSum.mCalc() );

    // Generate parton showers where appropriate.
    if (doFSRinRes) {
      process = event;
      process.init("(hard process)", &particleData);
      partonLevel.setupShowerSys( process, event);
      partonLevel.resonanceShowers( process, event, true);
    }

    // Generate hadronization and decays.
    bool status = (doHadronLevel) ? forceHadronLevel() : true;
    if (status) infoPrivate.addCounter(4);
    if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
    if (doFSRinRes && nPrevious < nShowProc) process.list(showSaV, showMaD);
    if (status && nPrevious < nShowEvt) event.list(showSaV, showMaD);
    endEvent(status ? PhysicsBase::COMPLETE : PhysicsBase::HADRONLEVEL_FAILED);
    return status;
  }

  // Reset arrays.
  infoPrivate.clear();
  weightContainer.clear();
  process.clear();
  event.clear();
  partonSystems.clear();
  beamSetup.clear();

  // Pick current beam valence flavours (for pi0, K0S, K0L, Pomeron).
  beamSetup.newValenceContent();

  // Recalculate kinematics when beam momentum spread.
  if (beamSetup.doMomentumSpread || beamSetup.doVarEcm
    || beamSetup.doVertexSpread) beamSetup.nextKinematics();

  // Simplified special treatment for low-energy nonperturbative collisions.
  if (doLowEnergy) {
    double eMinPertNow = eMinPert
      + 2. * max( 0., particleData.m0(beamSetup.idA) - particleData.m0(2212))
      + 2. * max( 0., particleData.m0(beamSetup.idB) - particleData.m0(2212));
    double pertRate = (beamSetup.eCM - eMinPertNow) / eWidthPert;
    if ( (doNonPert && !doSoftQCD)
      || ( beamSetup.doVarEcm && pertRate < 10
        && (pertRate <= 0 || exp( -pertRate ) > rndm.flat())) ) {
      bool nextNP = nextNonPert();

      // Optionally check final event for problems.
      if (nextNP && checkEvent && !check()) {
        logger.ERROR_MSG("check of event revealed problems");
        endEvent(PhysicsBase::CHECK_FAILED);
        return false;
      }
      endEvent(nextNP ? PhysicsBase::COMPLETE : PhysicsBase::LOWENERGY_FAILED);
      return nextNP;
    }
  }

  // Outer loop over hard processes; only relevant for user-set vetoes.
  for ( ; ; ) {

    infoPrivate.addCounter(10);
    bool hasVetoed = false;
    bool hasVetoedDiff = false;

    // Provide the hard process that starts it off. Only one try.
    infoPrivate.clear();
    process.clear();
    partonSystems.clear();

    // Reset the event information. Necessary if the previous event was read
    // from LHEF, while the current event is not read from LHEF.
    infoPrivate.setLHEF3EventInfo();

    if ( !processLevel.next( process, procType) ) {
      if (doLHA && infoPrivate.atEndOfFile())
        logger.ABORT_MSG("reached end of Les Houches Events File");
      else
        logger.ABORT_MSG("processLevel failed; giving up");
      endEvent(PhysicsBase::PROCESSLEVEL_FAILED);
      return false;
    }

    infoPrivate.addCounter(11);

    // Update tried and selected events immediately after next event was
    // generated. Note: This does not accumulate cross section.
    processLevel.accumulate(false);

    // Possibility for a user veto of the process-level event.
    if (doVetoProcess) {
      hasVetoed = userHooksPtr->doVetoProcessLevel( process);
      if (hasVetoed) {
        if (abortIfVeto) {
          endEvent(PhysicsBase::PROCESSLEVEL_USERVETO);
          return false;
        }
        continue;
      }
    }

    // Possibility to perform matrix element merging for this event.
    if (doMerging && mergingPtr) {
      int veto = mergingPtr->mergeProcess( process );
      // Apply possible merging scale cut.
      if (veto == -1) {
        if (abortIfVeto) {
          endEvent(PhysicsBase::MERGING_FAILED);
          return false;
        }
        continue;
      // Exit because of vanishing no-emission probability.
      } else if (veto == 0) {
        event = process;
        break;
      }

      // Redo resonance decays after the merging, in case the resonance
      // structure has been changed because of reclusterings.
      if (veto == 2 && doResDec) processLevel.nextDecays( process);
    }

    // Possibility to stop the generation at this stage.
    if (!doPartonLevel) {
      beamSetup.boostAndVertex( process, event, true, true);
      processLevel.accumulate();
      event.scale( process.scale() );
      event.scaleSecond( process.scaleSecond() );
      infoPrivate.addCounter(4);
      if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
      if (nPrevious < nShowInfo) infoPrivate.list();
      if (nPrevious < nShowProc) process.list(showSaV, showMaD);
      endEvent(PhysicsBase::COMPLETE);
      return true;
    }

    // Save spare copy of process record in case of problems.
    Event processSave = process;
    int sizeMPI       = infoPrivate.sizeMPIarrays();
    infoPrivate.addCounter(12);
    for (int i = 14; i < 19; ++i) infoPrivate.setCounter(i);

    // Allow up to ten tries for parton- and hadron-level processing.
    bool physical   = true;
    for (int iTry = 0; iTry < NTRY; ++iTry) {

      infoPrivate.addCounter(14);
      physical  = true;
      hasVetoed = false;

      // Restore original process record if problems.
      if (iTry > 0) {
        process = processSave;
        infoPrivate.resizeMPIarrays( sizeMPI);
      }

      // Reset event record and (extracted partons from) beam remnants.
      event.clear();
      beamSetup.clear();
      partonSystems.clear();

      // Parton-level evolution: ISR, FSR, MPI.
      if ( !partonLevel.next( process, event) ) {

        // Abort event generation if parton level is set to abort.
        if (infoPrivate.getAbortPartonLevel()) {
          endEvent(PhysicsBase::PARTONLEVEL_FAILED);
          return false;
        }

        // Skip to next hard process for failure owing to veto in merging.
        if (partonLevel.hasVetoedMerging()) {
          event = process;
          break;
        }

        // Skip to next hard process for failure owing to deliberate veto,
        // or alternatively retry for the same hard process.
        hasVetoed = partonLevel.hasVetoed();
        if (hasVetoed) {
          if (retryPartonLevel) {
            --iTry;
            continue;
          }
          if (abortIfVeto) {
            endEvent(PhysicsBase::PARTONLEVEL_FAILED);
            return false;
          }
          break;
        }

        // If hard diffractive event has been discarded retry partonLevel.
        hasVetoedDiff = partonLevel.hasVetoedDiff();
        if (hasVetoedDiff) {
          logger.WARNING_MSG(
            "discarding hard diffractive event from partonLevel; try again");
          break;
        }

        // Else make a new try for other failures.
        logger.ERROR_MSG("partonLevel failed; try again");
        physical = false;
        continue;
      }
      infoPrivate.addCounter(15);

      // Possibility for a user veto of the parton-level event.
      if (doVetoPartons) {
        hasVetoed = userHooksPtr->doVetoPartonLevel( event);
        if (hasVetoed) {
          if (abortIfVeto) {
            endEvent(PhysicsBase::PARTONLEVEL_USERVETO);
            return false;
          }
          break;
        }
      }

      // Boost to lab frame (before decays, for vertices).
      beamSetup.boostAndVertex( process, event, true, true);

      // Possibility to stop the generation at this stage.
      if (!doHadronLevel) {
        processLevel.accumulate();
        partonLevel.accumulate();
        event.scale( process.scale() );
        event.scaleSecond( process.scaleSecond() );
        // Optionally check final event for problems.
        if (checkEvent && !check()) {
          logger.ABORT_MSG("check of event revealed problems");
          endEvent(PhysicsBase::CHECK_FAILED);
          return false;
        }
        infoPrivate.addCounter(4);
        if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
        if (nPrevious < nShowInfo) infoPrivate.list();
        if (nPrevious < nShowProc) process.list(showSaV, showMaD);
        if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);
        endEvent(PhysicsBase::COMPLETE);
        return true;
      }

      // Hadron-level: hadronization, decays.
      infoPrivate.addCounter(16);
      if ( !hadronLevel.next( event) ) {
        // Check if we aborted due to user intervention.
        if ( canVetoHadronization && hadronLevel.hasVetoedHadronize() ) {
          endEvent(PhysicsBase::HADRONLEVEL_USERVETO);
          return false;
        }
        logger.ERROR_MSG("hadronLevel failed; try again");
        physical = false;
        continue;
      }

      // If R-hadrons have been formed, then (optionally) let them decay.
      if (decayRHadrons && rHadrons.exist() && !doRHadronDecays()) {
        logger.ERROR_MSG("decayRHadrons failed; try again");
        physical = false;
        continue;
      }
      infoPrivate.addCounter(17);

      // Optionally check final event for problems.
      if (checkEvent && !check()) {
        logger.ERROR_MSG("check of event revealed problems");
        physical = false;
        continue;
      }

      // Stop parton- and hadron-level looping if you got this far.
      infoPrivate.addCounter(18);
      break;
    }

    // If event vetoed then to make a new try.
    if (hasVetoed || hasVetoedDiff)  {
      if (abortIfVeto) {
        endEvent(PhysicsBase::PARTONLEVEL_USERVETO);
        return false;
      }
      continue;
    }

    // If event failed any other way (after ten tries) then give up.
    if (!physical) {
      logger.ABORT_MSG("parton+hadronLevel failed; giving up");
      endEvent(PhysicsBase::OTHER_UNPHYSICAL);
      return false;
    }

    // Process- and parton-level statistics. Event scale.
    processLevel.accumulate();
    partonLevel.accumulate();
    event.scale( process.scale() );
    event.scaleSecond( process.scaleSecond() );

    // End of outer loop over hard processes. Done with normal option.
    infoPrivate.addCounter(13);
    break;
  }

  // List events.
  if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
  if (nPrevious < nShowInfo) infoPrivate.list();
  if (nPrevious < nShowProc) process.list(showSaV,showMaD);
  if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);

  // Done.
  infoPrivate.addCounter(4);
  endEvent(PhysicsBase::COMPLETE);
  return true;

}

//--------------------------------------------------------------------------

// Set beam CM energy.

bool Pythia::setKinematics(double eCMIn) {

  // Update heavy ions if necessary.
  if (doHeavyIons)
    if (!heavyIonsPtr->setKinematics(eCMIn))
      return false;

  // Pass to BeamSetup object.
  return beamSetup.setKinematics(eCMIn);

}

//--------------------------------------------------------------------------

// Set beam energies.

bool Pythia::setKinematics(double eAIn, double eBIn) {

  // Update heavy ions if necessary.
  if (doHeavyIons)
    if (!heavyIonsPtr->setKinematics(eAIn, eBIn))
      return false;

  // Pass to BeamSetup object.
  return beamSetup.setKinematics(eAIn, eBIn);

}

//--------------------------------------------------------------------------

// Set beam momenta.

bool Pythia::setKinematics(double pxAIn, double pyAIn, double pzAIn,
  double pxBIn, double pyBIn, double pzBIn) {

  // Update heavy ions if necessary.
  if (doHeavyIons)
    if (!heavyIonsPtr->setKinematics(pxAIn, pyAIn, pzAIn, pxBIn, pyBIn, pzBIn))
      return false;

  // Pass to BeamSetup object.
  return beamSetup.setKinematics(pxAIn, pyAIn, pzAIn, pxBIn, pyBIn, pzBIn);

}

//--------------------------------------------------------------------------

// Set beam momenta.

bool Pythia::setKinematics(Vec4 pAIn, Vec4 pBIn) {

  // Update heavy ions if necessary.
  if (doHeavyIons)
    if (!heavyIonsPtr->setKinematics(pAIn, pBIn))
      return false;

  // Pass to BeamSetup object.
  return beamSetup.setKinematics(pAIn, pBIn);

}

//--------------------------------------------------------------------------

// Flexible-use call at the beginning of each event in pythia.next().

void Pythia::beginEvent() {

  // Loop through all PhysicsBase-derived objects.
  for ( auto physicsPtr : physicsPtrs ) physicsPtr->beginEvent();

  // Done.
  return;

}

//--------------------------------------------------------------------------

// Flexible-use call at the end of each event in pythia.next().

void Pythia::endEvent(PhysicsBase::Status status) {

  // Loop through all PhysicsBase-derived objects.
  for ( auto physicsPtr : physicsPtrs ) physicsPtr->endEvent(status);

  // Update the event weight by the Dire shower weight when relevant.
  // Code to be moved to the Dire endEvent method.
  /*
  if (useNewDire) {
    // Retrieve the shower weight.
    direPtr->weightsPtr->calcWeight(0.);
    direPtr->weightsPtr->reset();
    double pswt = direPtr->weightsPtr->getShowerWeight();
    // Multiply the shower weight to the event weight.
    double wt = infoPrivate.weight();
    infoPrivate.updateWeight(wt * pswt);
  }
  */

  // Done.
  return;

}

//--------------------------------------------------------------------------

// Possibility to set a pointer to a new object.

void Pythia::pushInfo() {

  for ( auto physicsPtr : physicsPtrs ) physicsPtr->initInfoPtr(infoPrivate);
}

//--------------------------------------------------------------------------

// Generate only the hadronization/decay stage, using internal machinery.
// The "event" instance should already contain a parton-level configuration.

bool Pythia::forceHadronLevel(bool findJunctions) {

  // Can only generate event if initialization worked.
  if (!isInit) {
    logger.ABORT_MSG("not properly initialized so cannot generate events");
    return false;
  }

  // Check whether any junctions in system. (Normally done in ProcessLevel.)
  // Avoid it if there are no final-state coloured partons.
  if (findJunctions) {
    event.clearJunctions();
    for (int i = 0; i < event.size(); ++i)
    if (event[i].isFinal()
    && (event[i].col() != 0 || event[i].acol() != 0)) {
      processLevel.findJunctions( event);
      break;
    }
  }

  // Allow for CR before the hadronization.
  if (forceHadronLevelCR) {

    // Setup parton system for SK-I and SK-II colour reconnection.
    // Require all final state particles to have the Ws as mothers.
    if (reconnectMode == 3 || reconnectMode == 4) {
      partonSystems.clear();
      partonSystems.addSys();
      partonSystems.addSys();
      partonSystems.setInRes(0, 3);
      partonSystems.setInRes(1, 4);
      for (int i = 5; i < event.size(); ++i) {
        if (event[i].mother1() - 3 < 0 || event[i].mother1() - 3 > 1) {
          logger.ERROR_MSG("event not set up correctly for SK-I or SK-II CR");
          return false;
        }
        partonSystems.addOut(event[i].mother1() - 3,i);
      }
    }

    // save spare copy of event in case of failure.
    Event spareEvent = event;
    bool colCorrect = false;

    // Allow up to ten tries for CR.
    for (int iTry = 0; iTry < NTRY; ++ iTry) {
      if ( stringInteractionsPtr->getColourReconnections() )
         stringInteractionsPtr->getColourReconnections()->next(event, 0);
      if (junctionSplitting.checkColours(event)) {
        colCorrect = true;
        break;
      }
      else event = spareEvent;
    }

    if (!colCorrect) {
      logger.ERROR_MSG("colour reconnection failed");
      return false;
    }
  }

  // Save spare copy of event in case of failure.
  Event spareEvent = event;

  // Allow up to ten tries for hadron-level processing.
  bool physical = true;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {
    physical = true;

    // Check whether any resonances need to be handled at process level.
    if (doResDec) {
      process = event;
      processLevel.nextDecays( process);

      // Allow for showers if decays happened at process level.
      if (process.size() > event.size()) {
        if (doFSRinRes) {
          partonLevel.setupShowerSys( process, event);
          partonLevel.resonanceShowers( process, event, false);
        } else event = process;
      }
    }

    // Hadron-level: hadronization, decays.
    if (hadronLevel.next( event)) break;

    // Failure might be by user intervention. Then we should not try again.
    if (canVetoHadronization && hadronLevel.hasVetoedHadronize()) {
      endEvent(PhysicsBase::HADRONLEVEL_USERVETO);
      break;
    }

    // If failure then warn, restore original configuration and try again.
    logger.WARNING_MSG("hadronLevel failed; try again");
    physical = false;
    event    = spareEvent;
  }

  // Done for simpler option.
  if (!physical)  {
    logger.ERROR_MSG("hadronLevel failed; giving up");
    return false;
  }

  // Optionally check final event for problems.
  if (checkEvent && !check()) {
    logger.ERROR_MSG("check of event revealed problems");
    return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Simplified treatment for low-energy nonperturbative collisions.

bool Pythia::nextNonPert(int procType) {

  // Fill collision instate.
  process.append( 90, -11, 0, 0, 0, 0, 0, 0,
    Vec4(0., 0., 0., beamSetup.eCM), beamSetup.eCM, 0. );
  process.append(beamSetup.idA, -12, 0, 0, 0, 0, 0, 0,
    Vec4(0., 0., beamSetup.pzAcm, beamSetup.eA), beamSetup.mA, 0. );
  process.append(beamSetup.idB, -12, 0, 0, 0, 0, 0, 0,
    Vec4(0., 0., beamSetup.pzBcm, beamSetup.eB), beamSetup.mB, 0. );
  for (int i = 0; i < 3; ++i) event.append( process[i] );

  // Pick process type if it has not already been set.
  if (procType == 0) procType = hadronLevel.pickLowEnergyProcess(
    beamSetup.idA, beamSetup.idB, beamSetup.eCM, beamSetup.mA, beamSetup.mB);
  int procCode = 150 + min( 9, abs(procType));

  if (procType == 0) {
    logger.ERROR_MSG("unable to pick process");
    return false;
  }

  // Do a low-energy collision.
  if (!doLowEnergyProcess( 1, 2, procType)) {
    logger.ERROR_MSG("low energy process failed");
    return false;
  }

  // Boost back to frame of collision.
  beamSetup.boostAndVertex( process, event, true, true);

  // Do hadron level.
  if (doHadronLevel) {
    if (!hadronLevel.next(event)) {
      logger.ERROR_MSG("further hadron level processes failed");
      return false;
    }
  }

  // Set event info.
  string procName ="Low-energy ";
  if      (procCode == 151) procName += "nonDiffractive";
  else if (procCode == 152) procName += "elastic";
  else if (procCode == 153) procName += "single diffractive (XB)";
  else if (procCode == 154) procName += "single diffractive (AX)";
  else if (procCode == 155) procName += "double diffractive";
  else if (procCode == 157) procName += "excitation";
  else if (procCode == 158) procName += "annihilation";
  else if (procCode == 159) procName += "resonant";
  infoPrivate.setType( procName, procCode, 0, (procCode == 151), false,
    (procCode == 153 || procCode == 155),
    (procCode == 154 || procCode == 155));

  // List events.
  int nPrevious = infoPrivate.getCounter(3) - 1;
  if (doLHA && nPrevious < nShowLHA) lhaUpPtr->listEvent();
  if (nPrevious < nShowInfo) infoPrivate.list();
  if (nPrevious < nShowProc) process.list(showSaV,showMaD);
  if (nPrevious < nShowEvt)  event.list(showSaV, showMaD);

  // Done.
  infoPrivate.addCounter(4);
  return true;

}

//--------------------------------------------------------------------------

// Perform R-hadron decays, either as part of normal evolution or forced.

bool Pythia::doRHadronDecays( ) {

  // Check if R-hadrons exist to be processed.
  if ( !rHadrons.exist() ) return true;

  // Do the R-hadron decay itself.
  if ( !rHadrons.decay( event) ) return false;

  // Perform showers in resonance decay chains.
  if ( !partonLevel.resonanceShowers( process, event, false) ) return false;

  // Subsequent hadronization and decays.
  if ( !hadronLevel.next( event) ) return false;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Print statistics on event generation.

void Pythia::stat() {

  if ( doHeavyIons ) {
    heavyIonsPtr->stat();
    return;
  }

  // Read out settings for what to include.
  bool showPrL = flag("Stat:showProcessLevel");
  bool showPaL = flag("Stat:showPartonLevel");
  bool showErr = flag("Stat:showErrors");
  bool reset   = flag("Stat:reset");

  // Statistics on cross section and number of events.
  if (doProcessLevel) {
    if (showPrL) processLevel.statistics(false);
    if (reset)   processLevel.resetStatistics();
  }

  // Statistics from other classes, currently multiparton interactions.
  if (showPaL) partonLevel.statistics(false);
  if (reset)   partonLevel.resetStatistics();

  // Merging statistics.
  if (doMerging && mergingPtr) mergingPtr->statistics();

  // Summary of which and how many warnings/errors encountered.
  if (showErr) logger.errorStatistics();
  if (reset)   logger.errorReset();

  // Loop through all PhysicsBase-derived objects.
  for ( auto physicsPtr : physicsPtrs ) physicsPtr->stat();

}

//--------------------------------------------------------------------------

// Write the Pythia banner, with symbol and version information.

void Pythia::banner() {

  // Read in version number and last date of change.
  double versionNumber = parm("Pythia:versionNumber");
  int versionDate = mode("Pythia:versionDate");
  string month[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

  // Get date and time.
  time_t t = time(0);
  char dateNow[12];
  strftime(dateNow,12,"%d %b %Y",localtime(&t));
  char timeNow[9];
  strftime(timeNow,9,"%H:%M:%S",localtime(&t));

  cout << "\n"
       << " *-------------------------------------------"
       << "-----------------------------------------* \n"
       << " |                                           "
       << "                                         | \n"
       << " |  *----------------------------------------"
       << "--------------------------------------*  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   PPP   Y   Y  TTTTT  H   H  III    A  "
       << "    Welcome to the Lund Monte Carlo!  |  | \n"
       << " |  |   P  P   Y Y     T    H   H   I    A A "
       << "    This is PYTHIA version " << fixed << setprecision(3)
       << setw(5) << versionNumber << "      |  | \n"
       << " |  |   PPP     Y      T    HHHHH   I   AAAAA"
       << "    Last date of change: " << setw(2) << versionDate%100
       << " " << month[ min(11, (versionDate/100)%100 - 1) ]
       << " " << setw(4) << versionDate/10000 <<  "  |  | \n"
       << " |  |   P       Y      T    H   H   I   A   A"
       << "                                      |  | \n"
       << " |  |   P       Y      T    H   H  III  A   A"
       << "    Now is " << dateNow << " at " << timeNow << "    |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Program documentation and an archive "
       << "of historic versions is found on:     |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                               https://p"
       << "ythia.org/                            |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   PYTHIA is authored by a collaboration"
       << " consisting of:                       |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Javira Altmann, Christian Bierlich, N"
       << "aomi Cooke, Nishita Desai,            |  | \n"
       << " |  |   Leif Gellersen, Ilkka Helenius, Phili"
       << "p Ilten, Leif Lonnblad,               |  | \n"
       << " |  |   Stephen Mrenna, Christian Preuss, Tor"
       << "bjorn Sjostrand, Peter Skands,        |  | \n"
       << " |  |   Marius Utheim, and Rob Verheyen.     "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   The complete list of authors, includi"
       << "ng contact information and            |  | \n"
       << " |  |   affiliations, can be found on https:/"
       << "/pythia.org/.                         |  | \n"
       << " |  |   Problems or bugs should be reported "
       << "on email at authors@pythia.org.        |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   The main program reference is C. Bier"
       << "lich et al,                           |  | \n"
       << " |  |   'A comprehensive guide to the physics"
       << " and usage of Pythia 8.3',            |  | \n"
       << " |  |   SciPost Phys. Codebases 8-r8.3 (2022)"
       << " [arXiv:2203.11601 [hep-ph]]          |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   PYTHIA is released under the GNU Gene"
       << "ral Public Licence version 2 or later.|  | \n"
       << " |  |   Please respect the MCnet Guidelines f"
       << "or Event Generator Authors and Users. |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Disclaimer: this program comes withou"
       << "t any guarantees.                     |  | \n"
       << " |  |   Beware of errors and use common sense"
       << " when interpreting results.           |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   Copyright (C) 2024 Torbjorn Sjostrand"
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  *----------------------------------------"
       << "--------------------------------------*  | \n"
       << " |                                           "
       << "                                         | \n"
       << " *-------------------------------------------"
       << "-----------------------------------------* \n" << endl;

}

//--------------------------------------------------------------------------

// Check for lines in file that mark the beginning of new subrun.

int Pythia::readSubrun(string line, bool warn) {

  // If empty line then done.
  int subrunLine = SUBRUNDEFAULT;
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos)
    return subrunLine;

  // If first character is not a letter, then done.
  string lineNow = line;
  int firstChar = lineNow.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalpha(lineNow[firstChar])) return subrunLine;

  // Replace an equal sign by a blank to make parsing simpler.
  while (lineNow.find("=") != string::npos) {
    int firstEqual = lineNow.find_first_of("=");
    lineNow.replace(firstEqual, 1, " ");
  }

  // Get first word of a line.
  istringstream splitLine(lineNow);
  string name;
  splitLine >> name;

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");
  }

  // Convert to lowercase. If no match then done.
  if (toLower(name) != "main:subrun") return subrunLine;

  // Else find new subrun number and return it.
  splitLine >> subrunLine;
  if (!splitLine) {
    if (warn) cout << "\n PYTHIA Warning: Main:subrun number not"
        << " recognized; skip:\n   " << line << endl;
    subrunLine = SUBRUNDEFAULT;
  }
  return subrunLine;

}

//--------------------------------------------------------------------------

// Check for lines in file that mark the beginning or end of commented section.
// Return +1 for beginning, -1 for end, 0 else.

int Pythia::readCommented(string line) {

  // If less than two nontrivial characters on line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return 0;
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  if (int(line.size()) < firstChar + 2) return 0;

  // If first two nontrivial characters are /* or */ then done.
  if (line.substr(firstChar, 2) == "/*") return +1;
  if (line.substr(firstChar, 2) == "*/") return -1;

  // Else done.
  return 0;

}

//--------------------------------------------------------------------------

// Check that the final event makes sense: no unknown id codes;
// charge and energy-momentum conserved.

bool Pythia::check() {

  // Reset.
  bool physical     = true;
  bool listVertices = false;
  bool listHistory  = false;
  bool listSystems  = false;
  bool listBeams    = false;
  iErrId.resize(0);
  iErrCol.resize(0);
  iErrEpm.resize(0);
  iErrNan.resize(0);
  iErrNanVtx.resize(0);
  Vec4 pSum;
  double chargeSum  = 0.;

  // Incoming beams counted with negative momentum and charge.
  if (doProcessLevel || doNonPert) {
    // Incoming particles will be at position "1" and "2" in most cases.
    // However, need to be careful how to find incoming particles after
    // QED radiation in DIS-type collisions. Thus, first find both incoming
    // particles.
    int iA = 1;
    int iB = 2;
    if (!(beamSetup.beamA2gamma || beamSetup.beamB2gamma)) {
      if (beamA.isLepton() && beamB.isHadron())
        { iA = beamA[0].iPos(); iB = 2; }
      if (beamB.isLepton() && beamA.isHadron())
        { iB = beamB[0].iPos(); iA = 1; }
      int iPos = 0;
      while ( beamA.isHadron() && iPos < beamB.size()
           && beamA.id() == beamB[iPos].id() )
        { iA = beamA[iPos].iPos(); iPos++;}
      iPos = 0;
      while ( beamB.isHadron() && iPos < beamB.size()
           && beamB.id() == beamB[iPos].id() )
        { iB = beamB[iPos].iPos(); iPos++; }
    }
    // Count incoming momentum and charge.
    pSum      = - (event[iA].p() + event[iB].p());
    chargeSum = - (event[1].charge() + event[2].charge());

  // If no ProcessLevel then sum final state of process record.
  } else if (process.size() > 0) {
    pSum = - process[0].p();
    for (int i = 0; i < process.size(); ++i)
      if (process[i].isFinal()) chargeSum -= process[i].charge();

  // If process not filled, then use outgoing primary in event.
  } else {
    pSum = - event[0].p();
    for (int i = 1; i < event.size(); ++i)
      if (event[i].statusAbs() < 10 || event[i].statusAbs() == 23)
        chargeSum -= event[i].charge();
  }
  double eLab = abs(pSum.e());

  // Loop over particles in the event.
  for (int i = 0; i < event.size(); ++i) {

    // Look for any unrecognized particle codes.
    int id = event[i].id();
    if (id == 0 || !particleData.isParticle(id)) {
      ostringstream errCode;
      errCode << ", i = " << i << ", id = " << id;
      logger.ERROR_MSG("unknown particle code", errCode.str());
      physical = false;
      iErrId.push_back(i);

    // Check that colour assignments are the expected ones.
    } else {
      int colType = event[i].colType();
      int col     = event[i].col();
      int acol    = event[i].acol();
      if ( event[i].statusAbs() / 10 == 8 ) acol = col = 0;
      if ( (colType ==  0 && (col  > 0 || acol  > 0))
        || (colType ==  1 && (col <= 0 || acol  > 0))
        || (colType == -1 && (col  > 0 || acol <= 0))
        || (colType ==  2 && (col <= 0 || acol <= 0)) ) {
        ostringstream errCode;
        errCode << ", i = " << i << ", id = " << id << " cols = " << col
                << " " << acol;
        logger.ERROR_MSG("incorrect colours", errCode.str());
        physical = false;
        iErrCol.push_back(i);
      }
    }

    // Some intermediate shower partons excepted from (E, p, m) consistency.
    bool checkMass = event[i].statusAbs() != 49 && event[i].statusAbs() != 59;

    // Look for particles with mismatched or not-a-number energy/momentum/mass.
    if (isfinite(event[i].p()) && isfinite(event[i].m())) {
      double errMass = abs(event[i].mCalc() - event[i].m())
        / max( 1.0, event[i].e());
      if (checkMass && errMass > mTolErr) {
        logger.ERROR_MSG("unmatched particle energy/momentum/mass");
        physical = false;
        iErrEpm.push_back(i);
      } else if (checkMass && errMass > mTolWarn) {
        logger.WARNING_MSG("not quite matched particle energy/momentum/mass");
      }
    } else {
      logger.ERROR_MSG("not-a-number energy/momentum/mass");
      physical = false;
      iErrNan.push_back(i);
    }

    // Look for particles with not-a-number vertex/lifetime.
    if (isfinite(event[i].vProd()) && isfinite(event[i].tau())) ;
    else {
      logger.ERROR_MSG("not-a-number vertex/lifetime");
      physical     = false;
      listVertices = true;
      iErrNanVtx.push_back(i);
    }

    // Add final-state four-momentum and charge.
    if (event[i].isFinal()) {
      pSum      += event[i].p();
      chargeSum += event[i].charge();
    }

  // End of particle loop.
  }

  // Check energy-momentum/charge conservation.
  double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
    + abs(pSum.pz());
  if (epDev > epTolErr * eLab) {
    logger.ERROR_MSG("energy-momentum not conserved");
    physical = false;
  } else if (epDev > epTolWarn * eLab) {
    logger.WARNING_MSG("energy-momentum not quite conserved");
  }
  if (abs(chargeSum) > 0.1) {
    logger.ERROR_MSG("charge not conserved");
    physical = false;
  }

  // Check that beams and event records agree on incoming partons.
  // Only meaningful for resolved beams.
  if (infoPrivate.isResolved() && !info.hasUnresolvedBeams())
  for (int iSys = 0; iSys < beamA.sizeInit(); ++iSys) {
    int eventANw  = partonSystems.getInA(iSys);
    int eventBNw  = partonSystems.getInB(iSys);
    // For photon sub-beams make sure to use correct beams.
    int beamANw   = ( beamA.getGammaMode() == 0 || !beamSetup.beamA2gamma
                 || (beamA.getGammaMode() == 2 && beamB.getGammaMode() == 2)) ?
                 beamA[iSys].iPos() : beamSetup.beamGamA[iSys].iPos();
    int beamBNw   = ( beamB.getGammaMode() == 0 || !beamSetup.beamB2gamma
                 || (beamB.getGammaMode() == 2 && beamA.getGammaMode() == 2)) ?
                 beamB[iSys].iPos() : beamSetup.beamGamB[iSys].iPos();
    if (eventANw != beamANw || eventBNw != beamBNw) {
      logger.ERROR_MSG("event and beams records disagree");
      physical    = false;
      listSystems = true;
      listBeams   = true;
    }
  }

  // Check that mother and daughter information match for each particle.
  vector<int> noMot;
  vector<int> noDau;
  vector< pair<int,int> > noMotDau;
  if (checkHistory) {

    // Loop through the event and check that there are beam particles.
    bool hasBeams = false;
    for (int i = 0; i < event.size(); ++i) {
      int status = event[i].status();
      if (abs(status) == 12) hasBeams = true;

      // Check that mother and daughter lists not empty where not expected to.
      vector<int> mList = event[i].motherList();
      vector<int> dList = event[i].daughterList();
      if (mList.size() == 0 && abs(status) != 11 && abs(status) != 12)
        noMot.push_back(i);
      if (dList.size() == 0 && status < 0 && status != -11)
        noDau.push_back(i);

      // Check that the particle appears in the daughters list of each mother.
      for (int j = 0; j < int(mList.size()); ++j) {
        if ( event[mList[j]].daughter1() <= i
          && event[mList[j]].daughter2() >= i ) continue;
        vector<int> dmList = event[mList[j]].daughterList();
        bool foundMatch = false;
        for (int k = 0; k < int(dmList.size()); ++k)
        if (dmList[k] == i) {
          foundMatch = true;
          break;
        }
        if (!hasBeams && mList.size() == 1 && mList[0] == 0) foundMatch = true;
        if (!foundMatch) {
          bool oldPair = false;
          for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == mList[j] && noMotDau[k].second == i) {
            oldPair = true;
            break;
          }
          if (!oldPair) noMotDau.push_back( make_pair( mList[j], i) );
        }
      }

      // Check that the particle appears in the mothers list of each daughter.
      for (int j = 0; j < int(dList.size()); ++j) {
        if ( event[dList[j]].statusAbs() > 80
          && event[dList[j]].statusAbs() < 90
          && event[dList[j]].mother1() <= i
          && event[dList[j]].mother2() >= i) continue;
        vector<int> mdList = event[dList[j]].motherList();
        bool foundMatch = false;
        for (int k = 0; k < int(mdList.size()); ++k)
        if (mdList[k] == i) {
          foundMatch = true;
          break;
        }
        if (!foundMatch) {
          bool oldPair = false;
          for (int k = 0; k < int(noMotDau.size()); ++k)
          if (noMotDau[k].first == i && noMotDau[k].second == dList[j]) {
            oldPair = true;
            break;
          }
          if (!oldPair) noMotDau.push_back( make_pair( i, dList[j]) );
        }
      }
    }

    // Warn if any errors were found.
    if (noMot.size() > 0 || noDau.size() > 0 || noMotDau.size() > 0) {
      logger.ERROR_MSG("mismatch in daughter and mother lists");
      physical    = false;
      listHistory = true;
    }
  }

  // Done for sensible events.
  if (physical) return true;

  // Print (the first few) flawed events: local info.
  if (nErrEvent < nErrList) {
    cout << "\n PYTHIA erroneous event info: \n";
    if (iErrId.size() > 0) {
      cout << " unknown particle codes in lines ";
      for (int i = 0; i < int(iErrId.size()); ++i)
        cout << iErrId[i] << " ";
      cout << "\n";
    }
    if (iErrCol.size() > 0) {
      cout << " incorrect colour assignments in lines ";
      for (int i = 0; i < int(iErrCol.size()); ++i)
        cout << iErrCol[i] << " ";
      cout << "\n";
    }
    if (iErrEpm.size() > 0) {
      cout << " mismatch between energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrEpm.size()); ++i)
        cout << iErrEpm[i] << " ";
      cout << "\n";
    }
    if (iErrNan.size() > 0) {
      cout << " not-a-number energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrNan.size()); ++i)
        cout << iErrNan[i] << " ";
      cout << "\n";
    }
    if (iErrNanVtx.size() > 0) {
      cout << " not-a-number vertex/lifetime in lines ";
      for (int i = 0; i < int(iErrNanVtx.size()); ++i)
        cout << iErrNanVtx[i] << " ";
      cout << "\n";
    }
    if (epDev > epTolErr * eLab) cout << scientific << setprecision(3)
      << " total energy-momentum non-conservation = " << epDev << "\n";
    if (abs(chargeSum) > 0.1) cout << fixed << setprecision(2)
      << " total charge non-conservation = " << chargeSum << "\n";
    if (noMot.size() > 0) {
      cout << " missing mothers for particles ";
      for (int i = 0; i < int(noMot.size()); ++i) cout << noMot[i] << " ";
      cout << "\n";
    }
    if (noDau.size() > 0) {
      cout << " missing daughters for particles ";
      for (int i = 0; i < int(noDau.size()); ++i) cout << noDau[i] << " ";
      cout << "\n";
    }
    if (noMotDau.size() > 0) {
      cout << " inconsistent history for (mother,daughter) pairs ";
      for (int i = 0; i < int(noMotDau.size()); ++i)
        cout << "(" << noMotDau[i].first << "," << noMotDau[i].second << ") ";
      cout << "\n";
    }

    // Print (the first few) flawed events: standard listings.
    infoPrivate.list();
    event.list(listVertices, listHistory);
    if (listSystems) partonSystems.list();
    if (listBeams) beamSetup.list();
  }

  // Update error counter. Done also for flawed event.
  ++nErrEvent;
  return false;

}

//==========================================================================

} // end namespace Pythia8
