// HeavyIons.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HeavyIons.h header) for the
// heavy ion classes classes, and some related global functions.

#include "Pythia8/BeamShape.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8/HINucleusModel.h"
#include "Pythia8/HISubCollisionModel.h"

namespace Pythia8 {

//==========================================================================

// The abstract HeavyIons class

//--------------------------------------------------------------------------

// Before doing anything Pythia should add special heavy ion versions
// for some groups of settings.

void HeavyIons::addSpecialSettings(Settings & settings) {
  setupSpecials(settings, "Diffraction:");
  setupSpecials(settings, "MultipartonInteractions:");
  setupSpecials(settings, "PDF:");
  setupSpecials(settings, "SigmaDiffractive:");
  setupSpecials(settings, "BeamRemnants:");
}

void HeavyIons::setupSpecials(Settings & settings, string match) {
  map<string,Flag> flags = settings.getFlagMap(match);
  for ( map<string,Flag>::iterator it = flags.begin();
        it != flags.end(); ++it )
    settings.addFlag("HI" + it->second.name, it->second.valDefault);
  map<string,Mode> modes = settings.getModeMap(match);
  for ( map<string,Mode>::iterator it = modes.begin();
        it != modes.end(); ++it )
    settings.addMode("HI" + it->second.name, it->second.valDefault,
                     it->second.hasMin, it->second.hasMax,
                     it->second.valMin, it->second.valMax, it->second.optOnly);
  map<string,Parm> parms = settings.getParmMap(match);
  for ( map<string,Parm>::iterator it = parms.begin();
        it != parms.end(); ++it )
    settings.addParm("HI" + it->second.name, it->second.valDefault,
                 it->second.hasMin, it->second.hasMax,
                 it->second.valMin, it->second.valMax);
  map<string,Word> words = settings.getWordMap(match);
  for ( map<string,Word>::iterator it = words.begin();
        it != words.end(); ++it )
    settings.addWord("HI" + it->second.name, it->second.valDefault);
  map<string,FVec> fvecs = settings.getFVecMap(match);
  for ( map<string, FVec>::iterator it = fvecs.begin();
        it != fvecs.end(); ++it )
    settings.addFVec("HI" + it->second.name, it->second.valDefault);
  map<string,MVec> mvecs = settings.getMVecMap(match);
  for ( map<string,MVec>::iterator it = mvecs.begin();
        it != mvecs.end(); ++it )
    settings.addMVec("HI" + it->second.name, it->second.valDefault,
                 it->second.hasMin, it->second.hasMax,
                 it->second.valMin, it->second.valMax);
  map<string,PVec> pvecs = settings.getPVecMap(match);
  for ( map<string,PVec>::iterator it = pvecs.begin();
        it != pvecs.end(); ++it )
    settings.addPVec("HI" + it->second.name, it->second.valDefault,
                 it->second.hasMin, it->second.hasMax,
                 it->second.valMin, it->second.valMax);
  map<string,WVec> wvecs = settings.getWVecMap(match);
  for ( map<string,WVec>::iterator it = wvecs.begin();
        it != wvecs.end(); ++it )
    settings.addWVec("HI" + it->second.name, it->second.valDefault);
}

void HeavyIons::setupSpecials(Pythia & p, string match) {
  Settings & opts = p.settings;
  map<string, Flag> flags = opts.getFlagMap(match);
  for ( map<string, Flag>::iterator it = flags.begin();
        it != flags.end(); ++it )
    opts.flag(it->second.name.substr(2), it->second.valNow, true);
  map<string, Mode> modes = opts.getModeMap(match);
  for ( map<string, Mode>::iterator it = modes.begin();
        it != modes.end(); ++it )
    opts.mode(it->second.name.substr(2), it->second.valNow, true);
  map<string, Parm> parms = opts.getParmMap(match);
  for ( map<string, Parm>::iterator it = parms.begin();
        it != parms.end(); ++it )
    opts.parm(it->second.name.substr(2), it->second.valNow, true);
  map<string, Word> words = opts.getWordMap(match);
  for ( map<string, Word>::iterator it = words.begin();
       it != words.end(); ++it )
    opts.word(it->second.name.substr(2), it->second.valNow, true);
  map<string, FVec> fvecs = opts.getFVecMap(match);
  for ( map<string, FVec>::iterator it = fvecs.begin();
        it != fvecs.end(); ++it )
    opts.fvec(it->second.name.substr(2), it->second.valNow, true);
  map<string, MVec> mvecs = opts.getMVecMap(match);
  for ( map<string, MVec>::iterator it = mvecs.begin();
        it != mvecs.end(); ++it )
    opts.mvec(it->second.name.substr(2), it->second.valNow, true);
  map<string, PVec> pvecs = opts.getPVecMap(match);
  for ( map<string, PVec>::iterator it = pvecs.begin();
        it != pvecs.end(); ++it )
    opts.pvec(it->second.name.substr(2), it->second.valNow, true);
  map<string, WVec> wvecs = opts.getWVecMap(match);
  for ( map<string, WVec>::iterator it = wvecs.begin();
        it != wvecs.end(); ++it )
    opts.wvec(it->second.name.substr(2), it->second.valNow, true);
}

//--------------------------------------------------------------------------

// Reset all process level settings in the given Pythia object. NOTE
// must be expanded if new process groups are included in Pythia.

void HeavyIons::clearProcessLevel(Pythia & pyt) {
  string path = pyt.settings.word("xmlPath");
  pyt.settings.mode("Tune:ee", 0);
  pyt.settings.mode("Tune:pp", 0);
  pyt.settings.init(path + "QCDSoftProcesses.xml", true);
  pyt.settings.init(path + "QCDHardProcesses.xml", true);
  pyt.settings.init(path + "ElectroweakProcesses.xml", true);
  pyt.settings.init(path + "OniaProcesses.xml", true);
  pyt.settings.init(path + "TopProcesses.xml", true);
  pyt.settings.init(path + "FourthGenerationProcesses.xml", true);
  pyt.settings.init(path + "HiggsProcesses.xml", true);
  pyt.settings.init(path + "SUSYProcesses.xml", true);
  pyt.settings.init(path + "NewGaugeBosonProcesses.xml", true);
  pyt.settings.init(path + "LeftRightSymmetryProcesses.xml", true);
  pyt.settings.init(path + "LeptoquarkProcesses.xml", true);
  pyt.settings.init(path + "CompositenessProcesses.xml", true);
  pyt.settings.init(path + "HiddenValleyProcesses.xml", true);
  pyt.settings.init(path + "ExtraDimensionalProcesses.xml", true);
  pyt.settings.init(path + "DarkMatterProcesses.xml", true);
  pyt.settings.init(path + "SecondHardProcess.xml", true);
  pyt.settings.init(path + "PhaseSpaceCuts.xml", true);
  // NOTE! if new processes are added in separate xml files these have
  // to be added here.
}

//--------------------------------------------------------------------------

// Update the Info object in the main Pythia object.

void HeavyIons::updateInfo() {
  *infoPtr =  hiInfo.primInfo;
  infoPtr->hiInfo = &hiInfo;
  infoPtr->weightContainerPtr->setWeightNominal(hiInfo.weight());
  infoPtr->sigmaReset();
  double norm = 1.0/double(hiInfo.NSave);
  int Nall = 0;
  double wall = 0.0;
  double w2all = 0.0;
  for ( map<int,int>::iterator ip = hiInfo.NPrim.begin();
        ip != hiInfo.NPrim.end(); ++ip ) {
    int N = ip->second;
    if ( !N ) continue;
    int pc = ip->first;
    double w = hiInfo.sumPrimW[pc]/millibarn;
    double w2 = hiInfo.sumPrimW2[pc]/pow2(millibarn);
    infoPtr->setSigma(pc, hiInfo.NamePrim[pc], N, N, N,
                      w*norm, sqrt(w2*norm)/N, w * millibarn);
    Nall += N;
    wall += w;
    w2all += w2;
  }
  infoPtr->setSigma(0, "sum", hiInfo.NSave, Nall, Nall,
                    wall*norm, sqrt(w2all*norm)/Nall, wall * millibarn);
}

//--------------------------------------------------------------------------

// Print out statistics from a HeavyIons run.

void HeavyIons::stat() {
  bool showPrL = flag("Stat:showProcessLevel");
  //  bool showPaL = settings.flag("Stat:showPartonLevel");
  bool showErr = flag("Stat:showErrors");
  bool reset   = flag("Stat:reset");
  Info & in = *infoPtr;
  // Header.
  if ( showPrL ) {
    cout << "\n *-----  HeavyIon Event and Cross Section Statistics  ------"
         << "-------------------------------------------------------*\n"
         << " |                                                            "
         << "                                                     |\n"
         << " | Subprocess                                    Code |       "
         << "     Number of events       |      sigma +- delta    |\n"
         << " |                                                    |       "
         << "Tried   Selected   Accepted |     (estimated) (mb)   |\n"
         << " |                                                    |       "
         << "                            |                        |\n"
         << " |------------------------------------------------------------"
         << "-----------------------------------------------------|\n"
         << " |                                                    |       "
         << "                            |                        |\n";

    vector<int> pc = in.codesHard();
    for ( int i = 0, N = pc.size(); i < N; ++i ) {
      cout << " | " << left << setw(45) << in.nameProc(pc[i])
           << right << setw(5) << pc[i] << " | "
           << setw(11) << in.nTried(pc[i]) << " "
           << setw(10) << in.nSelected(pc[i]) << " "
           << setw(10) << in.nAccepted(pc[i]) << " | "
           << scientific << setprecision(3)
           << setw(11) << in.sigmaGen(pc[i])
           << setw(11) << in.sigmaErr(pc[i]) << " |\n";
    }
    if ( pc.empty() ) in.setSigma(0, "sum", hiInfo.NSave, 0, 0, 0.0, 0.0, 0.0);

    cout << " |                                                    |       "
         << "                            |                        |\n"
         << " | " << left << setw(50) << "sum" << right << " | " << setw(11)
         << in.nTried(0) << " " << setw(10) << in.nSelected(0) << " "
         << setw(10) << in.nAccepted(0) << " | " << scientific
         << setprecision(3) << setw(11)
         << in.sigmaGen(0) << setw(11) << in.sigmaErr(0) << " |\n";
    cout << " | " << left << setw(50) << "(Estimated total cross section)"
         << right << " | " << setw(11)
         << hiInfo.nAttempts() << " " << setw(10) << 0 << " " << setw(10)
         << 0 << " | " << scientific << setprecision(3) << setw(11)
         << hiInfo.sigmaTot() << setw(11) << hiInfo.sigmaTotErr() << " |\n";
    cout << " | " << left << setw(50)
         << "(Estimated non-diffractive cross section)"
         << right << " | " << setw(11)
         << hiInfo.nAttempts() << " " << setw(10) << 0 << " " << setw(10)
         << 0 << " | " << scientific << setprecision(3) << setw(11)
         << hiInfo.sigmaND() << setw(11) << hiInfo.sigmaNDErr() << " |\n";
    // Listing finished.
    cout << " |                                                            "
         << "                                                     |\n"
         << " *-----  End HeavyIon Event and Cross Section Statistics -----"
         << "-----------------------------------------------------*" << endl;
  }
  if ( reset ) hiInfo = HIInfo();
  if ( showErr ) {
    for ( int i = 1, np = pythia.size(); i < np; ++i )
      loggerPtr->errorCombine(pythia[i]->logger, "(" + pythiaNames[i] + ")");
    loggerPtr->errorStatistics();
  }
  if ( reset ) loggerPtr->errorReset();

}

//--------------------------------------------------------------------------

// Check the settings and return false of there are no heavy ion beams.

bool HeavyIons::isHeavyIon(Settings & settings) {
  int idProj = settings.mode("Beams:idA");
  int idTarg = settings.mode("Beams:idB");
  return ( abs(idProj/100000000) == 10 ||abs(idTarg/100000000) == 10 );
}

//==========================================================================

// Angantyr is the main HeavyIons model in Pythia.

//--------------------------------------------------------------------------

// Constructor.

Angantyr::Angantyr(Pythia & mainPythiaIn)
  : HeavyIons(mainPythiaIn), hasSignal(true),
    collPtr(0), bGenPtr(0), projPtr(0), targPtr(0), recoilerMode(1), bMode(0),
    doAbort(false) {
  selectMB = make_shared<ProcessSelectorHook>();
  selectSASD = make_shared<ProcessSelectorHook>();
  pythia.resize(ALL);
  info.resize(ALL);
  pythiaNames.resize(ALL);
  pythiaNames[HADRON] = "HADRON";
  pythiaNames[MBIAS] = "MBIAS";
  pythiaNames[SASD] = "SASD";
  pythiaNames[SIGPP] = "SIGPP";
  pythiaNames[SIGPN] = "SIGPN";
  pythiaNames[SIGNP] = "SIGNP";
  pythiaNames[SIGNN] = "SIGNN";

}

//--------------------------------------------------------------------------

// Destructor deleting model objects that are not provided from the
// outside (via HIUserHooks).

Angantyr::~Angantyr() {
  for ( int i = MBIAS; i < ALL; ++i ) if ( pythia[i] ) delete pythia[i];
}

//--------------------------------------------------------------------------

// Add a HIUserHooks object to customise the Angantyr model.

bool Angantyr::setUserHooksPtr(PythiaObject sel, shared_ptr<UserHooks> uhook) {
  for ( int i = HADRON; i < ALL; ++i )
    if ( ( i == sel || ALL == sel ) && !pythia[i]->setUserHooksPtr(uhook) )
      return false;
  return true;
}

//--------------------------------------------------------------------------

// Figure out what beams the user wants.

void Angantyr::setBeamKinematics(int idA, int idB) {
  // We will use the MBIAS BeamSetup object to figure out what is
  // happening. Whatever we do here will be overridden when we do the
  // proper init().
  beamSetupPtr = pythia[MBIAS]->info.beamSetupPtr;
  pythia[MBIAS]->settings.mode("Beams:idA", idA);
  pythia[MBIAS]->settings.mode("Beams:idB", idB);
  beamSetupPtr->mA = particleDataPtr->m0(idA);
  beamSetupPtr->mB = particleDataPtr->m0(idB);
  beamSetupPtr->initFrame();
  unifyFrames();
}

//--------------------------------------------------------------------------

// Create an EventInfo object connected to a SubCollision from the
// last event generated by the given PythiaObject.

EventInfo Angantyr::mkEventInfo(Pythia & pyt, Info & infoIn,
                                const SubCollision * coll) {
    EventInfo ei;
    ei.coll = coll;
    ei.event = pyt.event;
    ei.info = infoIn;
    ei.code =  pyt.info.code();
    ei.ordering = ( ( HIHooksPtr && HIHooksPtr->hasEventOrdering() )?
                    HIHooksPtr->eventOrdering(ei.event, infoIn):
                    pyt.info.bMPI() );
    if ( coll ) {
      ei.projs[coll->proj] = make_pair(1, ei.event.size());
      ei.targs[coll->targ] = make_pair(2, ei.event.size());
    }

    ei.ok = true;
    return ei;
  }

//--------------------------------------------------------------------------

void Angantyr::banner(int idProj, int idTarg) const {

  string colOut = "              ";
  string cols = particleDataPtr->name(idProj)+" on "+
    particleDataPtr->name(idTarg);
  colOut.replace(colOut.begin(), colOut.begin() + cols.size(), cols);

  cout << " *----------------------  Initializing Angantyr  ----------------"
        << "------*\n"
        << " |                    We collide: " + colOut + "                 "
        << "      |\n"
        << " |                                                               "
        << "      |\n"
        << " |                    Below follows initialization               "
        << "      |\n"
        << " |                    of sub-collisions.                         "
        << "      |\n"
        << " |                                                               "
        << "      |\n"
        << " |                   //>________________________________         "
        << "      |\n"
        << " |          [########[]_________________________________>        "
        << "      |\n"
        << " |                   \\\\>                                       "
        << "        |\n";
  if (!settingsPtr->flag("HeavyIon:SigFitPrint"))
    cout << " *-------------------------------------------------------------"
          << "--------*" << endl;
  else
    cout << " |                                                             "
          << "        |" << endl;

}

//--------------------------------------------------------------------------

// Initialise Angantyr. Called from within Pythia::init().

bool Angantyr::init() {

  // Read settings.
  int idProj = mode("Beams:idA");
  int idTarg = mode("Beams:idB");
  doSDTest = flag("Angantyr:SDTest");
  glauberOnly = flag("Angantyr:GlauberOnly");
  recoilerMode = mode("Angantyr:SDRecoil");
  bMode = mode("Angantyr:impactMode");
  doVarECM = flag("Beams:allowVariableEnergy");
  doHadronLevel = flag("HadronLevel:all");

  int idProjP = idProj;
  int idProjN = 0;
  int idTargP = idTarg;
  int idTargN = 0;
  bool isHIProj = ( abs(idProj/100000000) == 10 );
  bool isHITarg = ( abs(idTarg/100000000) == 10 );
  bool isHI = isHIProj || isHITarg || mode("HeavyIon:mode") > 1;
  if ( isHIProj ) {
    idProjP = idProj > 0? 2212: -2212;
    idProjN = idProj > 0? 2112: -2112;
  }
  if ( isHITarg ) {
    idTargP = idTarg > 0? 2212: -2212;
    idTargN = idTarg > 0? 2112: -2112;
  }
  if ( mode("HeavyIon:mode") == 1 && !isHI ) {
    loggerPtr->ABORT_MSG("no heavy ions requested");
    return false;
  }

  bool print = flag("HeavyIon:showInit") && !settingsPtr->flag("Print:quiet");
  if ( print ) banner(idProj, idTarg);

  // Fix settings to be used for subobjects.
  settingsPtr->mode("Next:numberCount", 0);
  settingsPtr->mode("Next:numberShowLHA", 0);
  settingsPtr->mode("Next:numberShowInfo", 0);
  settingsPtr->mode("Next:numberShowProcess", 0);
  settingsPtr->mode("Next:numberShowEvent", 0);
  settingsPtr->flag("HadronLevel:all", false);
  settingsPtr->flag("SoftQCD:all", false);
  settingsPtr->flag("SoftQCD:elastic", false);
  settingsPtr->flag("SoftQCD:nonDiffractive", false);
  settingsPtr->flag("SoftQCD:singleDiffractive", false);
  settingsPtr->flag("SoftQCD:doubleDiffractive", false);
  settingsPtr->flag("SoftQCD:centralDiffractive", false);

  // Create Pythia subobjects.
  for ( int i = MBIAS; i < ALL; ++i ) {
    pythia[i] = new Pythia(*settingsPtr, *particleDataPtr, false);
    pythia[i]->settings.mode("HeavyIon:mode", 1);
    pythia[i]->settings.flag("Beams:allowVertexSpread", false);
    if (i != MBIAS)
      pythia[i]->settings.mode("MultipartonInteractions:reuseInit", 0);
  }

  // Allow for user to override with a custom HIUserHooks.
  if ( HIHooksPtr ) HIHooksPtr->init(idProj, idTarg);

  // Initialize kinematics and cross sections.
  setBeamKinematics(idProjP, idTargP);
  sigTotNN.init();
  sigTotNN.calc(idProjP, idTargP, beamSetupPtr->eCM);

  // Set up nucleus geometry.
  if (HIHooksPtr && HIHooksPtr->hasProjectileModel())
    projPtr = HIHooksPtr->projectileModel();
  else
    projPtr = NucleusModel::create(mode("Angantyr:NucleusModelA"));
  if (!projPtr) {
    loggerPtr->ABORT_MSG("nucleus model not found for projectile");
    return false;
  }
  projPtr->initPtr(idProj, true, *infoPtr);
  if (!projPtr->init()) {
    loggerPtr->ABORT_MSG("projectile nucleus model failed to initialize");
    return false;
  }
  projPtr->setPN(beamSetupPtr->pAinit);

  if (HIHooksPtr && HIHooksPtr->hasTargetModel())
    targPtr = HIHooksPtr->targetModel();
  else
    targPtr = NucleusModel::create(mode("Angantyr:NucleusModelB"));
  if (!targPtr) {
    loggerPtr->ABORT_MSG("nucleus model not found for target");
    return false;
  }
  targPtr->initPtr(idTarg, false, *infoPtr);
  if (!targPtr->init()) {
    loggerPtr->ABORT_MSG("target nucleus model failed to initialize");
    return false;
  }
  targPtr->setPN(beamSetupPtr->pBinit);

  // Set up subcollision model.
  if ( HIHooksPtr && HIHooksPtr->hasSubCollisionModel() )
    collPtr = HIHooksPtr->subCollisionModel();
  else
    collPtr = SubCollisionModel::create(mode("Angantyr:CollisionModel"));
  if (!collPtr) {
    loggerPtr->ABORT_MSG("subcollision model not found");
    return false;
  }
  collPtr->initPtr(*projPtr, *targPtr, sigTotNN, *settingsPtr,
                   *infoPtr, *rndmPtr);
  if (!collPtr->init(idProjP, idTargP, beamSetupPtr->eCM)) {
    loggerPtr->ABORT_MSG("subcollision model failed to initialize");
    return false;
  }
  hiInfo.avNDbSave = collPtr->avNDB();

  // Set up impact parameter generator.
  if ( HIHooksPtr && HIHooksPtr->hasImpactParameterGenerator() )
    bGenPtr = HIHooksPtr->impactParameterGenerator();
  else
    bGenPtr = make_shared<ImpactParameterGenerator>();
  bGenPtr->initPtr(*infoPtr, *collPtr, *projPtr, *targPtr);
  if ( !bGenPtr->init() ) {
    loggerPtr->ABORT_MSG("impact parameter generator failed to initialize");
    return false;
  }

  // Initialize subobject for minimum bias processes.
  clearProcessLevel(*pythia[MBIAS]);
  pythia[MBIAS]->settings.flag("SoftQCD:all", true);
  pythia[MBIAS]->settings.mode("Beams:idA", idProjP);
  pythia[MBIAS]->settings.mode("Beams:idB", idTargP);
  if ( beamSetupPtr->frameType > 3 ) {
    pythia[MBIAS]->settings.mode("Beams:eA", beamSetupPtr->eA);
    pythia[MBIAS]->settings.mode("Beams:eB", beamSetupPtr->eB);
    pythia[MBIAS]->settings.mode("Beams:frameType", 2);
  }

  pythia[MBIAS]->addUserHooksPtr(selectMB);
  init(MBIAS, "minimum bias processes");

  // Initialize subobject for secondary absorptive processes.
  clearProcessLevel(*pythia[SASD]);
  Settings & sdabsopts = pythia[SASD]->settings;
  sdabsopts.flag("SoftQCD:singleDiffractive", true);

  setupSpecials(*pythia[SASD], "HIDiffraction:");
  setupSpecials(*pythia[SASD], "HIMultipartonInteractions:");
  setupSpecials(*pythia[SASD], "HIPDF:");
  setupSpecials(*pythia[SASD], "HISigmaDiffractive:");
  setupSpecials(*pythia[SASD], "HIBeamRemnants:");
  if ( sdabsopts.mode("Angantyr:SASDmode") > 0 ) {
    double pT0Ref = sdabsopts.parm("MultipartonInteractions:pT0Ref");
    double ecmRef = sdabsopts.parm("MultipartonInteractions:ecmRef");
    double ecmPow = sdabsopts.parm("MultipartonInteractions:ecmPow");
    double ecm = beamSetupPtr->eCM;
    sdabsopts.parm("Beams:eCM", ecm);
    double pT0     = pT0Ref * pow(ecm / ecmRef, ecmPow);
    sdabsopts.parm("MultipartonInteractions:pT0Ref", pT0);
    sdabsopts.parm("MultipartonInteractions:ecmRef", ecm);
    sdabsopts.parm("MultipartonInteractions:ecmPow", 0.0);
    sdabsopts.word("PDF:PomSet", "11");
    int reuseMpi = settingsPtr->mode("HeavyIon:SasdMpiReuseInit");
    if (reuseMpi != 0) {
      string initFile = settingsPtr->word("HeavyIon:SasdMpiInitFile");
      sdabsopts.mode("MultipartonInteractions:reuseInit", reuseMpi);
      sdabsopts.word("MultipartonInteractions:initFile", initFile);
    }
    if ( sdabsopts.mode("Angantyr:SASDmode") == 2 ) {
      sdabsopts.parm("Diffraction:mRefPomP", ecm);
      double sigND = sigTotNN.sigmaND();
      double mmin = sdabsopts.parm("Diffraction:mMinPert");
      double powp = sdabsopts.parm("HIDiffraction:mPowPomP");
      sdabsopts.parm("Diffraction:mPowPomP", powp, true);
      if ( powp > 0.0 ) sigND /= ((1.0 - pow(mmin/ecm, powp))/powp);
      else sigND /= log(ecm/mmin);
      sdabsopts.parm("Diffraction:sigmaRefPomP", sigND, true);
    }
    if ( sdabsopts.mode("Angantyr:SASDmode") >= 3 ) {
      sdabsopts.parm("Diffraction:mRefPomP", ecm);
      double sigND = sigTotNN.sigmaND();
      sdabsopts.parm("Diffraction:sigmaRefPomP", sigND, true);
      sdabsopts.parm("Diffraction:mPowPomP", 0.0);
    }
  }
  sdabsopts.mode("Beams:idA", idProjP);
  sdabsopts.mode("Beams:idB", idTargP);
  if ( beamSetupPtr->frameType > 3 ) {
    sdabsopts.mode("Beams:eA", beamSetupPtr->eA);
    sdabsopts.mode("Beams:eB", beamSetupPtr->eB);
    sdabsopts.mode("Beams:frameType", 2);
  }

  pythia[SASD]->addUserHooksPtr(selectSASD);
  init(SASD, "secondary absorptive processes as single diffraction.");

  // Initialize subobject for hadronization.
  clearProcessLevel(*pythia[HADRON]);
  pythia[HADRON]->settings.flag("ProcessLevel:all", false);
  pythia[HADRON]->settings.flag("PartonLevel:all", false);
  pythia[HADRON]->settings.flag("HadronLevel:all", doHadronLevel);
  pythia[HADRON]->settings.mode("Beams:idA", idProj);
  pythia[HADRON]->settings.mode("Beams:idB", idTarg);

  // Initialize subobjects for signal processes.
  pythia[SIGPP]->settings.mode("Beams:idA", idProjP);
  pythia[SIGPP]->settings.mode("Beams:idB", idTargP);
  if ( idTargN ) {
    pythia[SIGPN]->settings.mode("Beams:idA", idProjP);
    pythia[SIGPN]->settings.mode("Beams:idB", idTargN);
  }
  if ( idProjN ) {
    pythia[SIGNP]->settings.mode("Beams:idA", idProjN);
    pythia[SIGNP]->settings.mode("Beams:idB", idTargP);
  }
  if ( idProjN && idTargN ) {
    pythia[SIGNN]->settings.mode("Beams:idA", idProjN);
    pythia[SIGNN]->settings.mode("Beams:idB", idTargN);
  }

  if ( hasSignal )
    hasSignal = pythia[SIGPP]->settings.hasHardProc() ||
      pythia[SIGPP]->settings.mode("Beams:frameType") >= 4;
  if ( hasSignal ) {
    init(SIGPP, "signal process (pp)", 10);
    if ( idTargN ) init(SIGPN, "signal process (pn)", 10);
    if ( idProjN ) init(SIGNP, "signal process (np)", 10);
    if ( idProjN && idTargN ) init(SIGNN, "signal process (nn)", 10);
  }

  if (doHadronLevel) {
    if ( print )
      cout << " Angantyr Info: Initializing hadronisation processes." << endl;
  }
  settingsPtr->flag("ProcessLevel:all", false);
  return true;

}

//--------------------------------------------------------------------------

// Initialize a specific Pythia object and optionally run a number
// of events to get a handle of the cross section.

bool Angantyr::init(PythiaObject sel, string name, int n) {
  bool print = flag("HeavyIon:showInit") && !flag("Print:quiet");
  shared_ptr<InfoGrabber> ihg = make_shared<InfoGrabber>();
  pythia[sel]->addUserHooksPtr(ihg);
  if ( print ) cout << " Angantyr Info: Initializing " << name << "." << endl;
  if ( !pythia[sel]->init() ) return false;
  info[sel] = ihg->getInfo();
  if ( n <= 0 ) return true;
  if ( print ) cout << "Generating a few signal events for " << name
                    << " to build up statistics" << endl;
  for ( int i = 0; i < 10; ++i ) pythia[sel]->next();
  return true;
}


//--------------------------------------------------------------------------

// Generate events and return EventInfo objects for different process
// types.

EventInfo Angantyr::getSignal(const SubCollision & coll) {
  if ( !hasSignal ) return EventInfo();
  int pytsel = SIGPP + coll.nucleons();
  int itry = MAXTRY;
  while ( itry-- ) {
    if ( pythia[pytsel]->next() )
      return mkEventInfo(*pythia[pytsel], *info[pytsel], &coll);
  }
  loggerPtr->WARNING_MSG("could not setup signal sub-collision");
  return EventInfo();
}

EventInfo Angantyr::getMBIAS(const SubCollision * coll, int procid) {
  int itry = MAXTRY;
  double bp = -1.0;
  if ( bMode > 0 && procid == 101 ) bp = coll->bp;
  HoldProcess hold(selectMB, procid, bp);
  while ( --itry ) {
    if ( !pythia[MBIAS]->next() ) continue;
    if (pythia[MBIAS]->info.code() != procid) {
      loggerPtr->ERROR_MSG("MBIAS info code not equal to set procid",
                          "contact the authors");
      doAbort = true;
    }
    return mkEventInfo(*pythia[MBIAS], *info[MBIAS], coll);
  }
  return EventInfo();
}

EventInfo Angantyr::getSASD(const SubCollision * coll, int procid) {
  int itry = MAXTRY;
  double bp = -1.0;
  if ( bMode > 1 ) bp = coll->bp;
  HoldProcess hold(selectSASD, procid, bp);
  while ( --itry ) {
    if ( !pythia[SASD]->next() ) continue;
    if (pythia[SASD]->info.code() != procid) {
      loggerPtr->ERROR_MSG("SASD info code not equal to set procid",
                          "contact the authors");
      doAbort = true;
    }
    return mkEventInfo(*pythia[SASD], *info[SASD], coll);
  }
  return EventInfo();
}

//--------------------------------------------------------------------------

// Generate primary absorptive (non-diffractive) nucleon-nucleon
// sub-collisions.

bool Angantyr::genAbs(SubCollisionSet& subCollsIn,
  list<EventInfo>& subEventsIn) {
  // The fully absorptive
  vector<const SubCollision*> abscoll;
   // The partly absorptive
  vector<const SubCollision*> abspart;
  // The non-diffractive and signal events
  multiset<EventInfo> ndeve, sigeve;

  // Select the primary absorptive sub collisions.
  for (const SubCollision& subColl : subCollsIn) {

    if ( subColl.type != SubCollision::ABS ) continue;
    if (!subColl.proj->done() && !subColl.targ->done() ) {
      abscoll.push_back(&subColl);
      if ( bMode > 0 ) {
        EventInfo ie = getND(subColl);
        if (ie.code != 101) {
          loggerPtr->ERROR_MSG("ND code not equal to 101",
                            "contact the authors");
          doAbort = true;
        }
        ndeve.insert(ie);
      }
      subColl.proj->select();
      subColl.targ->select();
    } else
      abspart.push_back(&subColl);
  }

  if ( abscoll.empty() ) return true;

  int Nabs = abscoll.size();
  int Nadd = abspart.size();

  if ( bMode == 0 ) {
    for ( int i = 0; i < Nabs + Nadd; ++i ) {
      EventInfo ie = getND();
      if (ie.code != 101) {
        loggerPtr->ERROR_MSG("ND code not equal to 101",
                            "contact the authors");
        doAbort = true;
      }
      ndeve.insert(ie);
    }
  }
  vector<int> Nii(4, 0);
  vector<double> w(4, 0.0);
  double wsum = 0.0;
  double P1 = 1.0;
  if ( hasSignal ) {

    // Count how many potential absorpitve collisions there are for
    // each iso-spin combination.
    for ( int i = 0, N = abscoll.size(); i < N; ++i )
      ++Nii[abscoll[i]->nucleons()];
    for ( int i = 0, N = abspart.size(); i < N; ++i )
      ++Nii[abspart[i]->nucleons()];

    if ( Nii[0] )
      w[0] = pythia[SIGPP]->info.sigmaGen()*millibarn/collPtr->sigND();
    if ( Nii[1] )
      w[1] = pythia[SIGPN]->info.sigmaGen()*millibarn/collPtr->sigND();
    if ( Nii[2] )
      w[2] = pythia[SIGNP]->info.sigmaGen()*millibarn/collPtr->sigND();
    if ( Nii[3] )
      w[3] = pythia[SIGNN]->info.sigmaGen()*millibarn/collPtr->sigND();

    wsum = Nii[0]*w[0] + Nii[1]*w[1] + Nii[2]*w[2] + Nii[3]*w[3];
    P1 = 1.0 - pow(1.0 - w[0], Nii[0])*pow(1.0 - w[1], Nii[1])*
               pow(1.0 - w[2], Nii[2])*pow(1.0 - w[3], Nii[3]);

  }

  bool noSignal = hasSignal;

  // *** THINK *** Is it ok to always pair the hardest events with the
  // *** most central sub-collisions, or will this introduce a strange
  // *** bias?
  multiset<EventInfo>::iterator it = ndeve.begin();
  EventInfo ei;
  for ( int i = 0, N = abscoll.size(); i < N; ++i ) {
    int b = abscoll[i]->nucleons();
    if ( Nii[b]
         && ( noSignal || w[b]*(wsum/P1 - 1.0)/(wsum - w[b]) > rndmPtr->flat())
         && (ei = getSignal(*abscoll[i])).ok ) {
      noSignal = false;
    }
    else
      ei =*it++;
    subEventsIn.push_back(ei);
    if ( !setupFullCollision(subEventsIn.back(), *abscoll[i],
                             Nucleon::ABS, Nucleon::ABS) )
      return false;
  }

  if ( noSignal ) return false;

  hiInfo.reweight(P1);

  return true;

}

//--------------------------------------------------------------------------

// Add secondary absorptive sub-collisions to the primary ones.

void Angantyr::addSASD(const SubCollisionSet& subCollsIn) {
  // Collect absorptively wounded nucleons in secondary
  // sub-collisions.
  int ntry = mode("Angantyr:SDTries");
  if ( settingsPtr->isMode("HI:SDTries") )
    ntry = mode("HI:SDTries");
  for (const SubCollision& subColl : subCollsIn)
    if ( subColl.type == SubCollision::ABS ) {
      if ( subColl.targ->done() && !subColl.proj->done() ) {
        EventInfo * evp = subColl.targ->event();
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getSDabsP(subColl);
          if ( addNucleonExcitation(*evp, add, true) ) {
            subColl.proj->select(*evp, Nucleon::ABS);
            break;
          }
          if ( itry == ntry - 1 ) hiInfo.failedExcitation();
        }
      } else if ( subColl.proj->done() && !subColl.targ->done() ) {
        EventInfo * evp = subColl.proj->event();
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getSDabsT(subColl);
          if ( addNucleonExcitation(*evp, add, true) ) {
            subColl.targ->select(*evp, Nucleon::ABS);
            break;
          }
          if ( itry == ntry - 1 ) hiInfo.failedExcitation();
        }
      }
    }
}

//--------------------------------------------------------------------------

// Add primary double diffraction sub-collisions.

bool Angantyr::addDD(const SubCollisionSet& subCollsIn,
  list<EventInfo>& subEventsIn) {
  // Collect full double diffraction collisions.
  for (const SubCollision& subColl : subCollsIn)
    if ( subColl.type == SubCollision::DDE &&
         !subColl.proj->done() && !subColl.targ->done() ) {
      subEventsIn.push_back(getDD(subColl));
      if ( !setupFullCollision(subEventsIn.back(), subColl,
                               Nucleon::DIFF, Nucleon::DIFF) )
        return false;
    }
  return true;
}

//--------------------------------------------------------------------------

// Add primary single diffraction sub-collisions.

bool Angantyr::addSD(const SubCollisionSet& subCollsIn,
  list<EventInfo> & subEventsIn) {
  // Collect full single diffraction collisions.
  for (const SubCollision& subColl : subCollsIn)
    if ( !subColl.proj->done() && !subColl.targ->done() ) {
      if ( subColl.type == SubCollision::SDEP ) {
        subEventsIn.push_back(getSDP(subColl));
        if ( !setupFullCollision(subEventsIn.back(), subColl,
                                 Nucleon::DIFF, Nucleon::ELASTIC) )
          return false;
      }
      if ( subColl.type == SubCollision::SDET ) {
        subEventsIn.push_back(getSDT(subColl));
        if ( !setupFullCollision(subEventsIn.back(), subColl,
                                 Nucleon::ELASTIC, Nucleon::DIFF) )
          return false;
      }
    }
  return true;
}

//--------------------------------------------------------------------------

// Add all secondary single diffractive sub-collisions to primary
// ones.

void Angantyr::addSDsecond(const SubCollisionSet& subCollsIn) {
  // Collect secondary single diffractive sub-collisions.
  int ntry = mode("Angantyr:SDTries");
  if ( settingsPtr->isMode("HI:SDTries") )  ntry = mode("HI:SDTries");
  for (const SubCollision& subColl : subCollsIn) {
    if ( !subColl.proj->done() &&
         ( subColl.type == SubCollision::SDEP ||
           subColl.type == SubCollision::DDE ) ) {
      EventInfo * evp = subColl.targ->event();
      for ( int itry = 0; itry < ntry; ++itry ) {
        EventInfo add = getSDP(subColl);
        if ( addNucleonExcitation(*evp, add, false) ) {
          subColl.proj->select(*evp, Nucleon::DIFF);
          break;
        }
        if ( itry == ntry - 1 ) hiInfo.failedExcitation();
      }
    }
    if ( !subColl.targ->done() &&
         ( subColl.type == SubCollision::SDET ||
           subColl.type == SubCollision::DDE ) ) {
      EventInfo * evp = subColl.proj->event();
      for ( int itry = 0; itry < ntry; ++itry ) {
        EventInfo add = getSDT(subColl);
        if ( addNucleonExcitation(*evp, add, false) ) {
          subColl.targ->select(*evp, Nucleon::DIFF);
          break;
        }
        if ( itry == ntry - 1 ) hiInfo.failedExcitation();
      }
    }
  }
}

//--------------------------------------------------------------------------

// Add all primary central diffraction sub-colliions

bool Angantyr::addCD(const SubCollisionSet& subCollsIn,
  list<EventInfo>& subEventsIn) {
  // Collect full central diffraction collisions.
  for (const SubCollision& subColl : subCollsIn)
    if ( subColl.type == SubCollision::CDE &&
         !subColl.proj->done() && !subColl.targ->done() ) {
      subEventsIn.push_back(getCD(subColl));
      if ( !setupFullCollision(subEventsIn.back(), subColl,
                               Nucleon::ELASTIC, Nucleon::ELASTIC) )
        return false;
    }
  return true;
}

//--------------------------------------------------------------------------

// Add all secondary central diffraction sub-colliions to primary
// ones.

void Angantyr::addCDsecond(const SubCollisionSet& subCollsIn) {
  // Collect secondary central diffractive sub-collisions.
  for (const SubCollision& subColl : subCollsIn) {
    if ( !subColl.proj->done() && subColl.type == SubCollision::CDE ) {
      EventInfo * evp = subColl.targ->event();
      EventInfo add = getCD(subColl);
      if ( addNucleonExcitation(*evp, add, false) ) {
        subColl.proj->select(*evp, Nucleon::ELASTIC);
      }
    }
    if ( !subColl.targ->done() && subColl.type == SubCollision::CDE ) {
      EventInfo * evp = subColl.proj->event();
      EventInfo add = getCD(subColl);
      if ( addNucleonExcitation(*evp, add, false) ) {
        subColl.targ->select(*evp, Nucleon::ELASTIC);
      }
    }
  }
}

//--------------------------------------------------------------------------

// Add all primary elastic sub-colliions

bool Angantyr::addEL(const SubCollisionSet& subCollsIn,
  list<EventInfo>& subEventsIn) {
  // Collect full elastic collisions.
  for (const SubCollision& subColl : subCollsIn)
    if ( subColl.type == SubCollision::ELASTIC &&
         !subColl.proj->done() && !subColl.targ->done() ) {
      subEventsIn.push_back(getEl(subColl));
      if (!setupFullCollision(subEventsIn.back(), subColl,
                              Nucleon::ELASTIC, Nucleon::ELASTIC))
        return false;
    }
  return true;
}

//--------------------------------------------------------------------------

// Add all secondary elastic sub-colliions to primary ones.

void Angantyr::addELsecond(const SubCollisionSet& subCollsIn) {
    // Collect secondary elastic sub-collisions.
  for (const SubCollision& subColl : subCollsIn) {
    if ( !subColl.proj->done() && subColl.type == SubCollision::ELASTIC ) {
      EventInfo * evp = subColl.targ->event();
      EventInfo add = getEl(subColl);
      if ( addNucleonExcitation(*evp, add, false) ) {
        subColl.proj->select(*evp, Nucleon::ELASTIC);
      }
    }
    if ( !subColl.targ->done() && subColl.type == SubCollision::ELASTIC ) {
      EventInfo * evp = subColl.proj->event();
      EventInfo add = getEl(subColl);
      if ( addNucleonExcitation(*evp, add, false) ) {
        subColl.targ->select(*evp, Nucleon::ELASTIC);
      }
    }
  }
}

//--------------------------------------------------------------------------

// Shift an event in impact parameter from the nucleon-nucleon
// sub-collision to the overall nucleus-nucleus frame. It is assumed
// that all partonic vertices are given in units of femtometers.

EventInfo & Angantyr::shiftEvent(EventInfo & ei) {
  if ( HIHooksPtr && HIHooksPtr->canShiftEvent() )
    return HIHooksPtr->shiftEvent(ei);

  double ymax = ei.event[1].y();
  Vec4 bmax = ei.coll->proj->bPos();
  double ymin = ei.event[2].y();
  Vec4 bmin = ei.coll->targ->bPos();
  for ( int i = 0, N = ei.event.size(); i < N; ++i ) {
    Vec4 shift = bmin + (bmax - bmin)*(ei.event[i].y() - ymin)/(ymax - ymin);
    ei.event[i].vProdAdd( shift * FM2MM);
  }
  return ei;
}

//--------------------------------------------------------------------------

// Prepare a primary sub-collision.

bool Angantyr::
setupFullCollision(EventInfo & ei, const SubCollision & coll,
                   Nucleon::Status projStatus, Nucleon::Status targStatus) {
  if ( !ei.ok ) return false;
  coll.proj->select(ei, projStatus);
  coll.targ->select(ei, targStatus);
  ei.coll = &coll;
  ei.projs.clear();
  ei.projs[coll.proj] = make_pair(1, ei.event.size());
  ei.targs.clear();
  ei.targs[coll.targ] = make_pair(2, ei.event.size());
  shiftEvent(ei);
  ei.event[1].status(-203);
  ei.event[1].mother1(1);
  ei.event[1].mother2(0);
  ei.event[2].status(-203);
  ei.event[2].mother1(2);
  ei.event[2].mother2(0);
  return fixIsoSpin(ei);
}

//--------------------------------------------------------------------------

// Trace a particle back to one of the beams in an event.

int Angantyr::getBeam(Event & ev, int i) {
  if ( int mom = ev[i].mother1() ) {
    if ( ev[mom].status() != -203 && ev[mom].mother1() < mom )
      return getBeam(ev, mom);
    else
      return mom;
  }
  else
    return i;
}

//--------------------------------------------------------------------------

// Minimum-bias sub-collisions are always generated as p-p events, and
// it is assumed to be safe to be assumed that they are iso-spin
// invariant so we can just modify the quark content in the remnants to
// get p-n, n-p, and n-n collisions.

bool Angantyr::fixIsoSpin(EventInfo & ei) {
  if ( HIHooksPtr && HIHooksPtr->canFixIsoSpin() )
    return HIHooksPtr->fixIsoSpin(ei);

  // Check if isospin needs fixing.
  int pshift = 0, tshift = 0;
  if ( ei.event[1].id() == 2212 && ei.coll->proj->id() == 2112 )
    pshift = 1;
  if ( ei.event[1].id() == -2212 && ei.coll->proj->id() == -2112 )
    pshift = -1;
  if ( pshift )
    ei.event[1].id(pshift*2112);
  if ( ei.event[2].id() == 2212 && ei.coll->targ->id() == 2112 )
    tshift = 1;
  if ( ei.event[2].id() == -2212 && ei.coll->targ->id() == -2112 )
    tshift = -1;
  if ( tshift )
    ei.event[2].id(tshift*2112);

  if ( !pshift && !tshift ) return true;

  // Try to find corresponding remnants that change flavour
  for ( int i = ei.event.size()  - 1; i > 2 && ( pshift || tshift ); --i ) {
    if ( pshift && ( isRemnant(ei, i) || ei.event[i].status() == 14 )
         &&  getBeam(ei.event, i) == 1 ) {
      int newid = 0;
      if ( ei.event[i].id() == 2*pshift ) newid = 1*pshift;
      if ( ei.event[i].id() == 2101*pshift ) newid = 1103*pshift;
      if ( ei.event[i].id() == 2103*pshift ) newid = 1103*pshift;
      if ( ei.event[i].id() == 2203*pshift ) newid = 2103*pshift;
      if ( ei.event[i].id() == 2212*pshift ) newid = 2112*pshift;
      if ( newid ) {
        ei.event[i].id(newid);
        pshift = 0;
        continue;
      }
    }
    if ( tshift && ( isRemnant(ei, i) || ei.event[i].status() == 14 )
         &&  getBeam(ei.event, i) == 2 ) {
      int newid = 0;
      if ( ei.event[i].id() ==    2*tshift ) newid =    1*tshift;
      if ( ei.event[i].id() == 2101*tshift ) newid = 1103*tshift;
      if ( ei.event[i].id() == 2103*tshift ) newid = 1103*tshift;
      if ( ei.event[i].id() == 2203*tshift ) newid = 2103*tshift;
      if ( ei.event[i].id() == 2212*tshift ) newid = 2112*tshift;
      if ( newid ) {
        ei.event[i].id(newid);
        tshift = 0;
        continue;
      }
    }
  }

  if ( !pshift && !tshift ) return true;

  // Try to find any final state quark that we modify, preferably far
  // in the beam direction.
  int qselp = 0;
  int qselt = 0;
  double yselp = 0.0;
  double yselt = 0.0;
  for ( int i = ei.event.size()  - 1; i > 2 && ( pshift || tshift ); --i ) {
    if ( pshift && ei.event[i].isFinal() && ei.event[i].id() == 2*pshift) {
      if ( ei.event[i].y() > yselp ) {
        qselp = i;
        yselp = ei.event[i].y();
      }
    }
    if ( tshift && ei.event[i].isFinal() && ei.event[i].id() == 2*tshift) {
      if ( ei.event[i].y() < yselt ) {
        qselt = i;
        yselt = ei.event[i].y();
      }
    }
  }
  if ( qselp ) {
    ei.event[qselp].id(1*pshift);
    pshift = 0;
  }
  if ( qselt ) {
    ei.event[qselt].id(1*tshift);
    tshift = 0;
  }

  return !pshift && !tshift;

}

//--------------------------------------------------------------------------

// Find recoilers in a primary sub-collisions to conserve energy and
// momentum when adding a secondary one. Not actually used yet.

vector<int> Angantyr::
findRecoilers(const Event & e, bool tside, int beam, int end,
              const Vec4 & pdiff, const Vec4 & pbeam) {
  vector<int> ret;
  multimap<double,int> ordered;
  double mtd2 = pdiff.m2Calc() + pdiff.pT2();
  int dir = tside? -1: 1;
  double ymax = -log(pdiff.pNeg());
  if ( tside ) ymax = -log(pdiff.pPos());
  for ( int i = beam, N = end; i < N; ++i )
    if ( e[i].status() > 0 )
      ordered.insert(make_pair(e[i].y()*dir, i));
  Vec4 prec;
  double pzfree2 = 0.0;
  multimap<double,int>::iterator it = ordered.begin();
  // cout << "--- find recoilers ----" << endl;
  // cout << setw(10) << setprecision(4)
  //      << log(pdiff.pPos())
  //      << setw(10) << pdiff.rap()
  //      << " diffractive system" << pdiff;
  while ( it != ordered.end() ) {
    if ( it->first > ymax ) break;
    int i = (*it++).second;
    Vec4 test = prec + e[i].p();
    double mtr2 = test.m2Calc() + test.pT2();
    double S = (pbeam + test).m2Calc();
    double pz2 = 0.25*(pow2(S - mtr2 - mtd2) - 4.0*mtr2*mtd2)/S;
    if ( pz2 < pzfree2 ) {
      // cout << setw(10) << setprecision(4) << it->first*dir << " failed "
      //      << (pz2 < 0.0? -sqrt(-pz2): sqrt(pz2)) << endl;
      break;
    }
    // cout << setw(10) << setprecision(4) << it->first*dir << " accept "
    //      << sqrt(pz2) << " (" << sqrt(S) << ")" << test;
    prec = test;
    pzfree2 = pz2;
    ret.push_back(i);
  }
  //  cout << "--- found recoilers ---" << endl;

  // *** THINK! *** Is this the best way?
  return ret;

}

//--------------------------------------------------------------------------

// Add a secondary sub-collision to a primary one.

bool Angantyr::addNucleonExcitation(EventInfo & ei, EventInfo & sub,
                                    bool colConnect) {
  fixIsoSpin(sub);
  shiftEvent(sub);
  if ( HIHooksPtr && HIHooksPtr->canAddNucleonExcitation() )
    return HIHooksPtr->addNucleonExcitation(ei, sub, colConnect);

  typedef map<Nucleon *, pair<int, int> >::iterator NucPos;
  bool tside = false;
  NucPos recnuc = ei.projs.find(sub.coll->proj);
  if ( recnuc != ei.projs.end() ) tside = true;
  NucPos rectarg = ei.targs.find(sub.coll->targ);
  if ( rectarg != ei.targs.end() ) {
    if ( tside ) loggerPtr->WARNING_MSG("nucleon already added");
    tside = false;
    recnuc = rectarg;
  }

  // First get the projectile system to take recoil and their momentum.
  int olddiff = tside? 4: 3;
  int beam = tside? 2: 1;
  int recbeam = recnuc->second.first;
  int recend = recnuc->second.second;
  Vec4 pbeam = sub.event[beam].p();
  Vec4 pdiff = sub.event[olddiff].p();
  if ( sub.code == 106 ) pdiff += sub.event[5].p();
  vector<int> rec;
  Vec4 prec;
  if ( HIHooksPtr && HIHooksPtr->canFindRecoilers() )
    rec = HIHooksPtr->findRecoilers(ei.event, tside, recbeam, recend,
                                    pdiff, pbeam);
  else if ( recoilerMode == 2 )
    rec = findRecoilers(ei.event, tside, recbeam, recend, pdiff, pbeam);
  else {
    if ( tside && ei.code == 104 && ei.event[4].status() > 0 )
      rec.push_back(4);
    else if ( !tside && ei.code == 103 && ei.event[3].status() > 0 )
      rec.push_back(3);
    else if ( tside && ei.event[3].status() > 0 &&
              ( ei.code == 102 || ei.code == 106 ) )
      rec.push_back(3);
    else if ( !tside && ei.event[4].status() > 0 &&
              ( ei.code == 102 || ei.code == 106 ) )
      rec.push_back(4);
    else
      for ( int i = recbeam, N = recend; i < N; ++i )
        if ( isRemnant(ei, i) && getBeam(ei.event, i) == recbeam )
          rec.push_back(i);
  }
  if ( rec.empty() ) return false;
  for ( int i = 0, N = rec.size(); i < N; ++i ) prec += ei.event[rec[i]].p();

  // Find the transform to the recoilers and the diffractive combined cms.
  pair<RotBstMatrix,RotBstMatrix> R12;
  if ( !getTransforms(prec, pdiff, pbeam, R12) )
    return false;

  // Transform the recoilers.
  for ( int i = 0, N = rec.size(); i < N; ++i )
    ei.event[rec[i]].rotbst(R12.first);

  // Copy the event and transform and offset the particles appropriately.

  int newbeam = ei.event.size();
  ei.event.append(sub.event[beam]);
  ei.event.back().status(-203);
  ei.event.back().mother1(beam);
  ei.event.back().mother2(0);
  ei.event.back().daughter1(ei.event.size());
  int newdiff = ei.event.size();
  int nextpos = 5;
  ei.event.append(sub.event[olddiff]);
  ei.event.back().rotbst(R12.second);
  ei.event.back().mother1(newbeam);
  ei.event.back().mother2(0);
  if ( sub.code == 102 ) {
    if ( tside )
      ei.targs[sub.coll->targ] = make_pair(newbeam, ei.event.size());
    else
      ei.projs[sub.coll->proj] = make_pair(newbeam, ei.event.size());
    return true;
  }

  int idoff = tside? newdiff - olddiff: newdiff - olddiff - 1;
  if ( sub.code == 106 ) {
    // Special handling of central diffraction.
    ++newdiff;
    ++nextpos;
    idoff = newdiff - 5;
    ei.event.append(sub.event[5]);
    ei.event.back().rotbst(R12.second);
    ei.event.back().mother1(newbeam);
    ei.event.back().mother2(0);
  }
  ei.event.back().daughter1(sub.event[olddiff].daughter1() + idoff);
  ei.event.back().daughter2(sub.event[olddiff].daughter2() + idoff);
  int coloff = ei.event.lastColTag();
  // Add energy to zeroth line and calculate new invariant mass.
  ei.event[0].p( ei.event[0].p() + pbeam );
  ei.event[0].m( ei.event[0].mCalc() );
  for (int i = nextpos; i < sub.event.size(); ++i) {
    Particle temp = sub.event[i];

    // Add offset to nonzero mother, daughter and colour indices.
    if ( temp.mother1() == olddiff ) temp.mother1(newdiff);
    else if ( temp.mother1() > 0 ) temp.mother1(temp.mother1() + idoff );
    if ( temp.mother2() == olddiff ) temp.mother2(newdiff);
    else if ( temp.mother2() > 0 ) temp.mother2(temp.mother2() + idoff );
    if ( temp.daughter1() > 0 ) temp.daughter1( temp.daughter1() + idoff );
    if ( temp.daughter2() > 0 ) temp.daughter2( temp.daughter2() + idoff );
    if ( temp.col() > 0 ) temp.col( temp.col() + coloff );
    if ( temp.acol() > 0 ) temp.acol( temp.acol() + coloff );
    temp.rotbst(R12.second);
    // Append particle to summed event.
    ei.event.append( temp );
  }

  addJunctions(ei.event, sub.event, coloff);

  if ( tside )
    ei.targs[sub.coll->targ] = make_pair(newbeam, ei.event.size());
  else
    ei.projs[sub.coll->proj] = make_pair(newbeam, ei.event.size());
  return true;

}

//--------------------------------------------------------------------------

// Calculate boosts to shuffle momenta when adding secondary
// sub-collisions.

bool
Angantyr::getTransforms(Vec4 prec, Vec4 pdiff, const Vec4 & pbeam,
                      pair<RotBstMatrix,RotBstMatrix> & R12) {
  RotBstMatrix Ri;
  Ri.toCMframe(pbeam, prec);
  Vec4 pr1 = prec;
  Vec4 pb1 = pbeam;
  Vec4 pd1 = pdiff;
  pr1.rotbst(Ri);
  pb1.rotbst(Ri);
  pd1.rotbst(Ri);
  Vec4 pr2 = pr1;
  if ( pd1.pT() >= abs(pr2.pz()) ) {
    return false;
  }
  double the = asin(pd1.pT()/abs(pr2.pz()));
  RotBstMatrix R1;
  R1.rot(the, pd1.phi());
  pr2.rotbst(R1);

  double S = (prec + pbeam).m2Calc();
  double mtr2 = pr2.pT2() + pr2.m2Calc();
  double mtd2 = pd1.pT2() + pd1.m2Calc();
  if ( sqrt(S) <= sqrt(mtr2) + sqrt(mtd2) ) {
    return false;
  }
  double z2 = 0.25*(mtr2*mtr2 + (mtd2 - S)*(mtd2 - S) - 2.0*mtr2*(mtd2 + S))/S;
  if ( z2 <= 0.0 ) {
    return false;
  }
  double z = sqrt(z2);
  double ppo2 = pow2(pr2.pNeg());
  double ppn2 = pow2(z + sqrt(z2 + mtr2));
  R1.bst(0.0, 0.0, -(ppn2 - ppo2)/(ppn2 + ppo2));

  ppo2 = pow2(pd1.pPos());
  ppn2 = pow2(z + sqrt(z2 + mtd2));
  RotBstMatrix R2;
  R2.bst(0.0, 0.0, (ppn2 - ppo2)/(ppn2 + ppo2));
  Vec4 pr3 = pr1;
  pr3.rotbst(R1);
  Vec4 pd3 = pd1;
  pd3.rotbst(R2);

  RotBstMatrix Rf = Ri;
  Rf.invert();
  Vec4 pr4 = pr3;
  pr4.rotbst(Rf);
  Vec4 pd4 = pd3;
  pd4.rotbst(Rf);

  R12.first = R12.second = Ri;
  R12.first.rotbst(R1);
  R12.second.rotbst(R2);
  R12.first.rotbst(Rf);
  R12.second.rotbst(Rf);
  prec.rotbst(R12.first);
  pdiff.rotbst(R12.second);

  return true;

}

//--------------------------------------------------------------------------

// Add sub-events together taking special care with the status of the
// incoming nucleons, and also handle the junctions correctly.

void Angantyr::addSubEvent(Event & evnt, Event & sub) {

  int idoff = evnt.size() - 1;
  int coloff = evnt.lastColTag();

  for (int i = 1; i < sub.size(); ++i) {
    Particle temp = sub[i];

    // Add offset to nonzero mother, daughter and colour indices.
    if ( temp.status() == -203 )
      temp.status(-13);
    else {
      if ( temp.mother1() > 0 ) temp.mother1(temp.mother1() + idoff );
      if ( temp.mother2() > 0 ) temp.mother2( temp.mother2() + idoff );
    }
    if ( temp.daughter1() > 0 ) temp.daughter1( temp.daughter1() + idoff );
    if ( temp.daughter2() > 0 ) temp.daughter2( temp.daughter2() + idoff );
    if ( temp.col() > 0 ) temp.col( temp.col() + coloff );
    if ( temp.acol() > 0 ) temp.acol( temp.acol() + coloff );
    // Append particle to summed event.
    evnt.append( temp );
  }

  addJunctions(evnt, sub, coloff);

}

void Angantyr::addJunctions(Event & ev, Event & addev, int coloff) {

  // Read out junctions one by one.
  Junction tempJ;
  int begCol, endCol;
  for (int i = 0; i < addev.sizeJunction(); ++i) {
    tempJ = addev.getJunction(i);

    // Add colour offsets to all three legs.
    for (int  j = 0; j < 3; ++j) {
      begCol = tempJ.col(j);
      endCol = tempJ.endCol(j);
      if (begCol > 0) begCol += coloff;
      if (endCol > 0) endCol += coloff;
      tempJ.cols( j, begCol, endCol);
    }
    // Append junction to summed event.
    ev.appendJunction( tempJ );
  }
}

//--------------------------------------------------------------------------

// Special function to generate secondary absorptive events as single
// diffraction. Called from Angantyr::next() and used for debugging
// and tuning purposes.

bool Angantyr::nextSASD(int procid) {
  Nucleon dummy;
  double bp = pythia[SASD]->parm("Angantyr:SDTestB");
  SubCollision coll(dummy, dummy, bp*collPtr->avNDB(), bp, SubCollision::ABS);
  EventInfo ei = getSASD(&coll, procid);
  if ( !ei.ok ) return false;
  pythia[HADRON]->event = ei.event;
  updateInfo();
  if (doHadronLevel) {
    if ( HIHooksPtr && HIHooksPtr->canForceHadronLevel() ) {
      if ( !HIHooksPtr->forceHadronLevel(*pythia[HADRON]) ) return false;
    } else {
      if ( !pythia[HADRON]->forceHadronLevel(false) ) return false;
    }
  }
  return true;
}

//--------------------------------------------------------------------------

// Take all sub-events and merge them together.

bool Angantyr::buildEvent(list<EventInfo> & subEventsIn) {
    Event & etmp = pythia[HADRON]->event;
    etmp.reset();
    etmp.append(projPtr->produceIon());
    etmp.append(targPtr->produceIon());
    etmp[0].p(etmp[1].p() + etmp[2].p());
    etmp[0].m(etmp[0].mCalc());
    double bx = 0.5*FM2MM*hiInfo.b()*cos(hiInfo.phi());
    double by = 0.5*FM2MM*hiInfo.b()*sin(hiInfo.phi());
    etmp[1].vProd( bx,  by, 0.0, 0.0);
    etmp[2].vProd(-bx, -by, 0.0, 0.0);

    // Start with the signal event(s)
    if ( hasSignal ) {
      bool found = false;
      for ( list<EventInfo>::iterator sit = subEventsIn.begin();
            sit != subEventsIn.end(); ++sit  ) {
        if ( sit->code >= 101 && sit->code <= 106 ) continue;
        addSubEvent(etmp, sit->event);
        hiInfo.select(sit->info);
        hiInfo.addSubCollision(*sit->coll);
        subEventsIn.erase(sit);
        found = true;
        break;
      }
      if ( !found ) {
        loggerPtr->ERROR_MSG("failed to generate signal event");
        return false;
      }
    } else
      hiInfo.select(subEventsIn.begin()->info);

    // Then all the others
    for ( list<EventInfo>::iterator sit = subEventsIn.begin();
          sit != subEventsIn.end(); ++sit  ) {
      addSubEvent(etmp, sit->event);
      hiInfo.addSubCollision(*sit->coll);
    }

    // Finally add all nucleon remnants.
    return addNucleusRemnants();

}

//--------------------------------------------------------------------------

// Construct nucleus remnants fron all non-interacting nucleons and
// add them to the main event.

bool Angantyr::addNucleusRemnants() {
  Event & etmp = pythia[HADRON]->event;
  int npp = 0;
  int nnp = 0;
  Vec4 ppsum;
  for (const Nucleon& nucleon : proj) {
    if (nucleon.event())
      hiInfo.addProjectileNucleon(nucleon);
    else {
      double e = pythia[HADRON]->parm("Beams:eA");
      double m = pythia[HADRON]->particleData.m0(nucleon.id());
      double pz = sqrt(max(e*e - m*m, 0.0));
      if ( nucleon.id() == 2212 ) {
        ++npp;
        ppsum += Vec4(0.0, 0.0, pz, e);
      } else if ( nucleon.id() == 2112 ) {
        ++nnp;
        ppsum += Vec4(0.0, 0.0, pz, e);
      } else
        etmp.append(nucleon.id(), 14, 1, 0, 0, 0, 0, 0, 0.0, 0.0, pz, e, m);
    }
  }
  int npt = 0;
  int nnt = 0;
  Vec4 tpsum;
  for (const Nucleon& nucleon : targ) {
    if (nucleon.event())
      hiInfo.addTargetNucleon(nucleon);
    else {
      double e = pythia[HADRON]->parm("Beams:eB");
      double m = pythia[HADRON]->particleData.m0(nucleon.id());
      double pz = -sqrt(max(e*e - m*m, 0.0));
      if ( nucleon.id() == 2212 ) {
        ++npt;
        tpsum += Vec4(0.0, 0.0, pz, e);
      } else if ( nucleon.id() == 2112 ) {
        ++nnt;
        tpsum += Vec4(0.0, 0.0, pz, e);
      } else
        etmp.append(nucleon.id(), 14, 2, 0, 0, 0, 0, 0, 0.0, 0.0, pz, e, m);
    }
  }

  Vec4 ptot = etmp[0].p();
  for ( int i = 0, N = etmp.size(); i < N; ++i )
    if ( etmp[i].status() > 0 ) ptot -= etmp[i].p();

  if ( npp + nnp +npt + nnt  == 0 ) return true;
  ParticleData & pdt = pythia[HADRON]->particleData;
  int idp = 0;
  if ( npp + nnp > 1 ) {
    idp = 1000000009 + 10000*npp + 10*(nnp + npp);
    pdt.addParticle(idp, "NucRem", 0, 3*npp, 0, ppsum.mCalc());
    pdt.particleDataEntryPtr(idp)->setHasChanged(false);
  }
  else if ( npp == 1 ) idp = 2212;
  else if ( nnp == 1 ) idp = 2112;
  int idt = 0;
  if ( npt + nnt > 1 ) {
    idt = 1000000009 + 10000*npt + 10*(nnt + npt);
    pdt.addParticle(idt, "NucRem", 0, 3*npt, 0, tpsum.mCalc());
    pdt.particleDataEntryPtr(idt)->setHasChanged(false);
  }
  else if ( npt == 1 ) idt = 2212;
  else if ( nnt == 1 ) idt = 2112;

  if ( npp + nnp > npt + nnt ) {
    if ( npt + nnt > 0 ) {
      etmp.append(idt, 14, 2, 0, 0, 0, 0, 0, tpsum, tpsum.mCalc());
      ptot -= tpsum;
    }
    etmp.append(idp, 14, 1, 0, 0, 0, 0, 0, ptot, ptot.mCalc());
  } else {
    if ( npp + nnp > 0 ) {
      etmp.append(idp, 14, 1, 0, 0, 0, 0, 0, ppsum, ppsum.mCalc());
      ptot -= ppsum;
    }
    etmp.append(idt, 14, 2, 0, 0, 0, 0, 0, ptot, ptot.mCalc());
  }
  return true;
}

//--------------------------------------------------------------------------

// Set beam kinematics.

bool Angantyr::setKinematics(){
  unifyFrames();
  if (!sigTotNN.calc(beamSetupPtr->idA, beamSetupPtr->idB, beamSetupPtr->eCM))
    return false;
  collPtr->updateSig();
  hiInfo.avNDbSave = collPtr->avNDB();
  collPtr->setKinematics(beamSetupPtr->eCM);
  bGenPtr->updateWidth();
  projPtr->setPN(beamSetupPtr->pAinit);
  targPtr->setPN(beamSetupPtr->pBinit);
  return true;
}

bool Angantyr::setKinematics(double eCMIn) {
  pythia[MBIAS]->setKinematics(eCMIn);
  if (!glauberOnly)
    pythia[SASD]->setKinematics(eCMIn);
  return setKinematics();
}

bool Angantyr::setKinematics(double eAIn, double eBIn) {
  pythia[MBIAS]->setKinematics(eAIn, eBIn);
  if (!glauberOnly)
    pythia[SASD]->setKinematics(eAIn, eBIn);
  return setKinematics();
}

bool Angantyr::setKinematics(double pxAIn, double pyAIn, double pzAIn,
  double pxBIn, double pyBIn, double pzBIn) {
  pythia[MBIAS]->setKinematics(pxAIn, pyAIn, pzAIn, pxBIn, pyBIn, pzBIn);
  if (!glauberOnly)
    pythia[SASD]->setKinematics(pxAIn, pyAIn, pzAIn, pxBIn, pyBIn, pzBIn);
  return setKinematics();
}

bool Angantyr::setKinematics(Vec4 pAIn, Vec4 pBIn) {
  pythia[MBIAS]->setKinematics(pAIn, pBIn);
  if (!glauberOnly)
    pythia[SASD]->setKinematics(pAIn, pBIn);
  return setKinematics();
}

//--------------------------------------------------------------------------

// Make sure the correct information is available irrespective of frame type.

void Angantyr::unifyFrames() {
  BeamSetup &bs = *beamSetupPtr;

  if ( bs.frameType == 1 ) {
    bs.eA     = bs.eB = bs.eCM/2;
    bs.pzA    =  sqrt(pow2(bs.eA) - pow2(bs.mA));
    bs.pzB    = -sqrt(pow2(bs.eB) - pow2(bs.mB));
    bs.pxA    = bs.pyA = bs.pxB = bs.pyB = 0.0;
    bs.pAinit = Vec4(bs.pxA, bs.pyA, bs.pzA, bs.eA);
    bs.pBinit = Vec4(bs.pxB, bs.pyB, bs.pzB, bs.eB);
  } else if ( bs.frameType == 3 ) {
    bs.eA     = sqrt(pow2(bs.pxA) + pow2(bs.pyA) + pow2(bs.pzA) + pow2(bs.mA));
    bs.eB     = sqrt(pow2(bs.pxB) + pow2(bs.pyB) + pow2(bs.pzB) + pow2(bs.mB));
    bs.pAinit = Vec4(bs.pxA, bs.pyA, bs.pzA, bs.eA);
    bs.pBinit = Vec4(bs.pxB, bs.pyB, bs.pzB, bs.eB);
    bs.eCM    = (bs.pAinit + bs.pBinit).mCalc();
  } else {
    bs.pzA    =  sqrt(pow2(bs.eA) - pow2(bs.mA));
    bs.pzB    = -sqrt(pow2(bs.eB) - pow2(bs.mB));
    bs.pxA    = bs.pyA = bs.pxB = bs.pyB = 0.0;
    bs.pAinit = Vec4(bs.pxA, bs.pyA, bs.pzA, bs.eA);
    bs.pBinit = Vec4(bs.pxB, bs.pyB, bs.pzB, bs.eB);
    bs.eCM    = (bs.pAinit + bs.pBinit).mCalc();
  }

  if ( !bs.doMomentumSpread ) {
    bs.pAnow = bs.pAinit;
    bs.pBnow = bs.pBinit;
  }

}

//--------------------------------------------------------------------------

// The main method called from Pythia::next().

bool Angantyr::next() {

  if (doSDTest)
    return nextSASD(104);

  int itry = MAXTRY;

  while ( itry-- && !doAbort) {

    // Generate impact parameter, nuclei, and sub-collisions.
    double bweight = 0.0;
    Vec4 bvec = bGenPtr->generate(bweight);
    proj = Nucleus(projPtr->generate(), bvec / 2.);
    targ = Nucleus(targPtr->generate(), -bvec / 2.);

    subColls = collPtr->getCollisions(proj, targ);
    hiInfo.addAttempt(subColls.T(), bvec.pT(), bvec.phi(), bweight);

    if ( subColls.empty() ) continue;
    if ( glauberOnly ) return true;

    list<EventInfo> subEvents;

    if ( !genAbs(subColls, subEvents) ) {
      loggerPtr->WARNING_MSG("could not setup signal or ND collisions");
      continue;
    }
    if ( hasSignal && subEvents.empty() ) continue;

    // Collect absorptively wounded nucleons in secondary sub-collisions.
    addSASD(subColls);

    // Collect full double diffraction collisions.
    if ( !addDD(subColls, subEvents) ) {
      loggerPtr->ERROR_MSG("could not setup DD sub-collision");
      continue;
    }

    // Collect full single diffraction collisions.
    if ( !addSD(subColls, subEvents) ) {
      loggerPtr->ERROR_MSG("could not setup SD sub-collision");
      continue;
    }

    // Collect secondary single diffractive sub-collisions.
    addSDsecond(subColls);

    // Collect full central diffraction collisions.
    if ( !addCD(subColls, subEvents) ) {
      loggerPtr->ERROR_MSG("could not setup CD sub-collision");
      continue;
    }

    // Collect secondary central diffractive sub-collisions.
    addCDsecond(subColls);

    // Collect full elastic collisions.
    if ( !addEL(subColls, subEvents) ) {
      loggerPtr->ERROR_MSG("could not setup elastic sub-collision");
      continue;
    }

    // Collect secondary elastic sub-collisions.
    addELsecond(subColls);

    // Finally bunch all events together.
    if ( subEvents.empty() ) continue;
    if ( !buildEvent(subEvents) ) continue;

    // Finally we hadronise everything, if requested.
    if (doHadronLevel) {
      if ( HIHooksPtr && HIHooksPtr->canForceHadronLevel() ) {
        if ( !HIHooksPtr->forceHadronLevel(*pythia[HADRON]) ) continue;
      } else {
        if ( !pythia[HADRON]->forceHadronLevel(false) ) continue;
      }
    }

    if ( settingsPtr->flag("Beams:allowVertexSpread") ) {
      pythia[HADRON]->getBeamShapePtr()->pick();
      Vec4 vertex = pythia[HADRON]->getBeamShapePtr()->vertex();
      for ( Particle & p : pythia[HADRON]->event ) p.vProdAdd( vertex);
    }

    hiInfo.accept();

    updateInfo();


    return true;

  }
  if (doAbort)
    loggerPtr->ABORT_MSG("Angantyr was aborted due to a critical error");
  else
    loggerPtr->ABORT_MSG("too many attempts to generate a working impact "
      "parameter point", "consider reducing HeavyIon:bWidth");
  hiInfo.reject();
  return false;

}

//==========================================================================

} // end namespace Pythia8
