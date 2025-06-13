// HeavyIons.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
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
    double w = hiInfo.sumPrimW[pc]*FMSQ2MB;
    double w2 = hiInfo.sumPrimW2[pc]*pow2(FMSQ2MB);
    infoPtr->setSigma(pc, hiInfo.NamePrim[pc], N, N, N,
                      w*norm, sqrt(w2*norm)/N, w * MB2FMSQ);
    Nall += N;
    wall += w;
    w2all += w2;
  }
  infoPtr->setSigma(0, "sum", hiInfo.NSave, Nall, Nall,
                    wall*norm, sqrt(w2all*norm)/Nall, wall * MB2FMSQ);
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
         << " | Primary NN sub-collision subprocess           Code |       "
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
    bool caveat = false;
    for ( int i = 0, N = pc.size(); i < N; ++i ) {
      string pname = in.nameProc(pc[i]);
      if ( pc[i] == 102 ) {
        pname += " (*)";
        caveat = true;
      }
      cout << " | " << left << setw(45) << pname
           << right << setw(5) << pc[i] << " | "
           << setw(11) << in.nTried(pc[i]) << " "
           << setw(10) << in.nSelected(pc[i]) << " "
           << setw(10) << in.nAccepted(pc[i]) << " | "
           << scientific << setprecision(3)
           << setw(11) << in.sigmaGen(pc[i])
           << setw(11) << in.sigmaErr(pc[i]) << " |\n";
    }
    if ( pc.empty() ) in.setSigma(0, "sum", hiInfo.NSave, 0, 0, 0.0, 0.0, 0.0);

    cout << " | " << left << setw(50) << "sum" << right << " | " << setw(11)
         << in.nTried(0) << " " << setw(10) << in.nSelected(0) << " "
         << setw(10) << in.nAccepted(0) << " | " << scientific
         << setprecision(3) << setw(11)
         << in.sigmaGen(0) << setw(11) << in.sigmaErr(0) << " |\n"
         << " |                                                    |       "
         << "                            |                        |\n";
    if ( caveat )
      cout << " | (*) Note: elastic events are not correctly treated |       "
           << "                            |                        |\n";
    cout << " |------------------------------------------------------------"
         << "-----------------------------------------------------|\n"
         << " |                                                            "
         << "                            |                        |\n";
    string line = "Semi-inclusive " + particleDataPtr->name(idProj) + " on " +
      particleDataPtr->name(idTarg) +
      " cross sections from the Glauber calculation:";
    cout << " | " << left << setw(86) << line
         << " |                        |\n"
         << " |                                                            "
         << "                            |                        |\n";
    cout << " | " << left << setw(86)
         << "Total" << " | "
         << right << scientific << setprecision(3)
         << setw(11) << hiInfo.glauberTot()
         << setw(11) << hiInfo.glauberTotErr() << " |\n";
    cout << " | " << left << setw(86)
         << "Non-Diffractive" << " | "
         << right << scientific << setprecision(3)
         << setw(11) << hiInfo.glauberND()
         << setw(11) << hiInfo.glauberNDErr() << " |\n";
    cout << " | " << left << setw(86)
         << "Total inelastic" << " | "
         << right << scientific << setprecision(3) << setw(11)
         << hiInfo.glauberINEL() << setw(11)
         << hiInfo.glauberINELErr() << " |\n";
    cout << " | " << left << setw(86)
         << "Elastic" << " | "
         << right << scientific << setprecision(3) << setw(11)
         << hiInfo.glauberEL() << setw(11)
         << hiInfo.glauberELErr() << " |\n";
    cout << " | " << left << setw(86)
         << "Diffractive target excitation" << " | "
         << right << scientific << setprecision(3) << setw(11)
         << hiInfo.glauberDiffT() << setw(11)
         << hiInfo.glauberDiffTErr() << " |\n";
    cout << " | " << left << setw(86)
         << "Diffractive projectile excitation" << " | "
         << right << scientific << setprecision(3) << setw(11)
         << hiInfo.glauberDiffP() << setw(11)
         << hiInfo.glauberDiffPErr() << " |\n";
    cout << " | " << left << setw(86)
         << "Double diffractive excitation" << " | "
         << right << scientific << setprecision(3) << setw(11)
         << hiInfo.glauberDDiff() << setw(11)
         << hiInfo.glauberDDiffErr() << " |\n";
     cout << " | " << left << setw(86)
         << "Elastic b-slope (GeV^-2)" << " | "
         << right << scientific << setprecision(3) << setw(11)
         << hiInfo.glauberBSlope() << setw(11)
         << hiInfo.glauberBSlopeErr() << " |\n";
    // Listing finished.
    cout << " |                                                            "
         << "                            |                        |\n"
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
  pythia[MBIAS]->settings.mode("Beams:idA", idA);
  pythia[MBIAS]->settings.mode("Beams:idB", idB);
  beamSetupPtr->mA = particleDataPtr->m0(idA);
  beamSetupPtr->mB = particleDataPtr->m0(idB);
  if ( idProj != idA ) {
    int A = (idProj/10)%1000;
    beamSetupPtr->mA = particleDataPtr->m0(idProj)/A;
  }
  if ( idTarg != idB ) {
    int A = (idTarg/10)%1000;
    beamSetupPtr->mB = particleDataPtr->m0(idTarg)/A;
  }
  beamSetupPtr->initFrame();
  unifyFrames();
}

//--------------------------------------------------------------------------

// Switch to new beam particle identities.
bool Angantyr::setBeamIDs(int idAIn, int idBIn) {

  if ( idAIn == projPtr->id() && ( idBIn == 0 || idBIn == targPtr->id() ) )
    return true;

  // Reset the statistics.
  hiInfo.glauberReset();

  // Set the projectile and target IDs.
  projPtr->setParticle(idAIn);
  if ( idBIn != 0 ) targPtr->setParticle(idBIn);

  // Set the beam IDs in minimum bias.
  if (!pythia[MBIAS]->setBeamIDs(projPtr->idN(), targPtr->idN()))
    return false;
  if (!pythia[SASD]->setBeamIDs(projPtr->idN(), targPtr->idN()))
    return false;

  // Calculate the total cross-section.
  sigTotNN.calc(projPtr->idN(), targPtr->idN(), beamSetupPtr->eCM);

  // Set masses and IDs.
  beamSetupPtr->mA = projPtr->mN();
  beamSetupPtr->mB = targPtr->mN();
  beamSetupPtr->idA = idAIn;
  beamSetupPtr->idB = idBIn;

  collPtr->setIDA(beamSetupPtr->represent(projPtr->idN()));
  bGenPtr->updateWidth();
  unifyFrames();

  idProj = idAIn;
  idTarg = idBIn;

  return true;
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

void Angantyr::banner() const {

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
  if (!settingsPtr->flag("HeavyIon:SigFitPrint") ||
       settingsPtr->mode("HeavyIon:SigFitNGen") <= 0 )
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
  idProj = mode("Beams:idA");
  idTarg = mode("Beams:idB");
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
  if ( print ) banner();

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
  for ( int i = MBIAS; i < ALL; ++i ) {
    pythia[i]->settings.mode("Beams:frameType", 1);
    pythia[i]->settings.parm("Beams:eCM", beamSetupPtr->eCM);
  }
  sigTotNN.init();
  if (!sigTotNN.calc(idProjP, idTargP, beamSetupPtr->eCM))
    return false;

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
    pythia[MBIAS]->settings.parm("Beams:eA", beamSetupPtr->eA);
    pythia[MBIAS]->settings.parm("Beams:eB", beamSetupPtr->eB);
    pythia[MBIAS]->settings.mode("Beams:frameType", 2);
  }

  pythia[MBIAS]->addUserHooksPtr(selectMB);
  init(MBIAS, "minimum bias processes");

  settingsPtr->wvec("Init:reuseMPIiDiffSys0",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys0"));
  settingsPtr->wvec("Init:reuseMPIiDiffSys1",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys1"));
  settingsPtr->wvec("Init:reuseMPIiDiffSys2",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys2"));
  settingsPtr->wvec("Init:reuseMPIiDiffSys3",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys3"));

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
    sdabsopts.parm("MultipartonInteractions:pT0Ref", pT0, true);
    sdabsopts.parm("MultipartonInteractions:ecmRef", ecm, true);
    sdabsopts.parm("MultipartonInteractions:ecmPow", 0.0, true);
    sdabsopts.word("PDF:PomSet", "11");
    int reuseMpi = settingsPtr->mode("HeavyIon:SasdMpiReuseInit");
    if (reuseMpi != 0) {
      string initFile = settingsPtr->word("HeavyIon:SasdMpiInitFile");
      sdabsopts.mode("MultipartonInteractions:reuseInit", reuseMpi);
      sdabsopts.word("MultipartonInteractions:initFile", initFile);
      sdabsopts.wvec("Init:reuseMPIiDiffSys0",
                     settingsPtr->wvec("Init:reuseSasdMPIiDiffSys0"));
      sdabsopts.wvec("Init:reuseMPIiDiffSys1",
                     settingsPtr->wvec("Init:reuseSasdMPIiDiffSys1"));
      sdabsopts.wvec("Init:reuseMPIiDiffSys2",
                     settingsPtr->wvec("Init:reuseSasdMPIiDiffSys2"));
      sdabsopts.wvec("Init:reuseMPIiDiffSys3",
                     settingsPtr->wvec("Init:reuseSasdMPIiDiffSys3"));
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

  settingsPtr->wvec("Init:reuseSasdMPIiDiffSys0",
                    sdabsopts.wvec("Init:reuseMPIiDiffSys0"));
  settingsPtr->wvec("Init:reuseSasdMPIiDiffSys1",
                    sdabsopts.wvec("Init:reuseMPIiDiffSys1"));
  settingsPtr->wvec("Init:reuseSasdMPIiDiffSys2",
                    sdabsopts.wvec("Init:reuseMPIiDiffSys2"));
  settingsPtr->wvec("Init:reuseSasdMPIiDiffSys3",
                    sdabsopts.wvec("Init:reuseMPIiDiffSys3"));


  // Initialize subobject for hadronization.
  clearProcessLevel(*pythia[HADRON]);
  pythia[HADRON]->settings.flag("ProcessLevel:all", false);
  pythia[HADRON]->settings.flag("PartonLevel:all", false);
  pythia[HADRON]->settings.flag("HadronLevel:all", doHadronLevel);
  pythia[HADRON]->settings.mode("Beams:idA", idProj);
  pythia[HADRON]->settings.mode("Beams:idB", idTarg);
  pythia[HADRON]->settings.flag("LowEnergyQCD:all", false);

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
    if ( pythia[pytsel]->next() ) {
      if ( pythia[pytsel]->event[0].pAbs2() != 0.0 )
        pythia[pytsel]->event.rotbst(toCMframe(pythia[pytsel]->event[1].p(),
                                               pythia[pytsel]->event[2].p()));
      return mkEventInfo(*pythia[pytsel], *info[pytsel], &coll);
    }
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
      w[0] = pythia[SIGPP]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();
    if ( Nii[1] )
      w[1] = pythia[SIGPN]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();
    if ( Nii[2] )
      w[2] = pythia[SIGNP]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();
    if ( Nii[3] )
      w[3] = pythia[SIGNN]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();

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
          if ( itry == ntry - 1 ) hiInfo.failedExcitation(subColl);
        }
      } else if ( subColl.proj->done() && !subColl.targ->done() ) {
        EventInfo * evp = subColl.proj->event();
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getSDabsT(subColl);
          if ( addNucleonExcitation(*evp, add, true) ) {
            subColl.targ->select(*evp, Nucleon::ABS);
            break;
          }
          if ( itry == ntry - 1 ) hiInfo.failedExcitation(subColl);
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
        if ( itry == ntry - 1 ) hiInfo.failedExcitation(subColl);
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
        if ( itry == ntry - 1 ) hiInfo.failedExcitation(subColl);
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

// Reset the main event.

void Angantyr::resetEvent() {

  Event & etmp = pythia[HADRON]->event;
  unifyFrames();
  etmp.reset();
  etmp.append(projPtr->produceIon());
  etmp.append(targPtr->produceIon());
  double mA = projPtr->mN();
  double mB = targPtr->mN();
  double eCM = beamSetupPtr->eCM;
  double pz = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
                           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;

  etmp[1].p(max(projPtr->A(), 1)*Vec4(0.0, 0.0, pz,
                                      sqrt(pow2(pz) + pow2(mA))));
  etmp[1].m(particleDataPtr->m0(idProj));
  etmp[2].p(max(targPtr->A(), 1)*Vec4(0.0, 0.0, -pz,
                                      sqrt(pow2(pz) + pow2(mB))));
  etmp[2].m(particleDataPtr->m0(idTarg));
  etmp[0].p(etmp[1].p() + etmp[2].p());
  etmp[0].m(etmp[0].mCalc());

}

//--------------------------------------------------------------------------

// Take all sub-events and merge them together.

bool Angantyr::buildEvent(list<EventInfo> & subEventsIn) {

  resetEvent();
  Event & etmp = pythia[HADRON]->event;
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
  }
  // Add statistics about participating nucleons and subcollisions.
  hiInfo.glauberStatistics();

  // Finally add all nucleon remnants.
  return addNucleusRemnants();
}

//--------------------------------------------------------------------------

// Construct nucleus remnants fron all non-interacting nucleons and
// add them to the main event.

bool Angantyr::addNucleusRemnants() {

  Event & etmp = pythia[HADRON]->event;
  BeamSetup & bs = *beamSetupPtr;
  ParticleData & pdt = pythia[HADRON]->particleData;

  // Get beam particle energies in rest frame.
  double eA = 0.5*(pow2(bs.eCM) + pow2(bs.mA) - pow2(bs.mB))/bs.eCM;
  double eB = bs.eCM - eA;


  // Sum up number of ineracted nucleons in the projectile.
  int npp = 0;
  int nnp = 0;
  for (const Nucleon& nucleon : proj)
    if (!nucleon.event()) {
      if ( abs(nucleon.id()) == 2212 ) ++npp;
      else if ( abs(nucleon.id()) == 2112 ) ++nnp;
      else etmp.append(nucleon.id(), 14, 1, 0, 0, 0, 0, 0,
                       0.0, 0.0, sqrt(pow2(eA) - pow2(bs.mA)), eA, bs.mA);
      }

  // Sum up number of ineracted nucleons in the target.
  int npt = 0;
  int nnt = 0;
  for (const Nucleon& nucleon : targ)
    if (!nucleon.event()) {
      if ( abs(nucleon.id()) == 2212 ) ++npt;
      else if ( abs(nucleon.id()) == 2112 ) ++nnt;
      else etmp.append(nucleon.id(), 14, 2, 0, 0, 0, 0, 0,
                       0.0, 0.0, -sqrt(pow2(eB) - pow2(bs.mB)), eB, bs.mB);
    }

  // sum up the missing momentum. Also sum up all remnant momenta.
  Vec4 ptot = etmp[0].p();
  vector<int> iRemP, iRemT;
  for ( int i = 0, N = etmp.size(); i < N; ++i ) {
    if ( etmp[i].status() <= 0 ) continue;
    ptot -= etmp[i].p();
    if ( etmp[i].status() != 63 ) continue;
    switch ( getBeam(etmp, i) ) {
    case 1: iRemP.push_back(i); break;
    case 2: iRemT.push_back(i); break;
    default: break;
    }
  }

  // Zero, one or two nucleus remnants?
  Vec4 pAr, pBr;
  int idAr = 0, idBr = 0;

  if ( npp + nnp > 1 ) idAr = 1000000009 + 10000*npp + 10*(nnp + npp);
  else if ( npp == 1 ) idAr = 2212;
  else if ( nnp == 1 ) idAr = 2112;
  if ( bs.idA < 0 ) idAr = -idAr;
  if ( npt + nnt > 1 ) idBr = 1000000009 + 10000*npt + 10*(nnt + npt);
  else if ( npt == 1 ) idBr = 2212;
  else if ( nnt == 1 ) idBr = 2112;
  if ( bs.idB < 0 ) idBr = -idBr;

  // Get masses and momenta of remnants
  double mAr = 0.0, mBr = 0.0;
  double mABr = ptot.mCalc();
  if ( idAr && idBr ) {
    mAr = (npp + nnp)*bs.mA;
    mBr = (npt + nnt)*bs.mB;
    double eAr = 0.5*(pow2(mABr) + pow2(mAr) - pow2(mBr))/mABr;
    double eBr = mABr - eAr;
    if ( eAr < mAr || eBr < mBr ) return false;
    auto M = fromCMframe(ptot);
    pAr = M*Vec4(0.0, 0.0,  sqrt(pow2(eAr) - pow2(mAr)), eAr);
    pBr = M*Vec4(0.0, 0.0, -sqrt(pow2(eBr) - pow2(mBr)), eBr);
  } else if ( idAr ) {
    pAr = ptot;
    mAr = pAr.mCalc();
  } else if ( idBr ) {
    pBr = ptot;
    mBr = pBr.mCalc();
  }


  // Add remnants.
  if ( idAr ) {
    if ( npp + nnp > 1 ) {
      pdt.addParticle(idAr, "NucRem", 0, 3*npp, 0, mAr);
      pdt.particleDataEntryPtr(idAr)->setHasChanged(false);
    }
    etmp.append(idAr, 14, 1, 0, 0, 0, 0, 0, pAr, mAr);
    ptot = Vec4();
  }
  if ( idBr ) {
    if ( npt + nnt > 1 ) {
      pdt.addParticle(idBr, "NucRem", 0, 3*npt, 0, mBr);
      pdt.particleDataEntryPtr(idBr)->setHasChanged(false);
    }
    etmp.append(idBr, 14, 2, 0, 0, 0, 0, 0, pBr, mBr);
    ptot = Vec4();
  }

  // Special case for ions, where all nucleons have interacted.
  // Select the most energetic remnant and dump excess momentum there.
  int irmax = 0;
  double ermax = 0.0;
  if ( projPtr->A() > 1 && idAr == 0 && iRemP.size() )
    for ( int i : iRemP )
      if ( etmp[i].e() > ermax ) {
        irmax = i;
        ermax = etmp[i].e();
      }
  if ( targPtr->A() > 1 && idBr == 0 && iRemT.size() )
    for ( int i : iRemT )
      if ( etmp[i].e() > ermax ) {
        irmax = i;
        ermax = etmp[i].e();
      }
  if ( irmax ) {
    etmp[irmax].p(ptot + etmp[irmax].p());
    etmp[irmax].m(etmp[irmax].mCalc());
    ptot = Vec4();
  }

  // Return successful.
  return true;

}

//--------------------------------------------------------------------------

// Set beam kinematics.

bool Angantyr::setKinematics(){
  unifyFrames();
  if (!sigTotNN.calc(projPtr->idN(), targPtr->idN(), beamSetupPtr->eCM))
    return false;
  collPtr->updateSig();
  hiInfo.avNDbSave = collPtr->avNDB();
  collPtr->setKinematics(beamSetupPtr->eCM);
  bGenPtr->updateWidth();
  projPtr->setPN(beamSetupPtr->pAinit);
  targPtr->setPN(beamSetupPtr->pBinit);
  return true;
}

bool Angantyr::setKinematicsCM() {
  hiInfo.glauberReset();
  if ( !setKinematics() ) return false;
  if (!glauberOnly && !pythia[SASD]->setKinematics(beamSetupPtr->eCM) )
    return false;
  return pythia[MBIAS]->setKinematics(beamSetupPtr->eCM);
}


bool Angantyr::setKinematics(double eCMIn) {
  if ( eCMIn == beamSetupPtr->eCM ) return true;
  if ( !beamSetupPtr->setKinematics(eCMIn) ) return false;
  return setKinematicsCM();
}

bool Angantyr::setKinematics(double eAIn, double eBIn) {
  if ( eAIn == beamSetupPtr->eA && eBIn == beamSetupPtr->eB )
    return true;
  if ( !beamSetupPtr->setKinematics(eAIn, eBIn) ) return false;
  return setKinematicsCM();
}

bool Angantyr::setKinematics(double pxAIn, double pyAIn, double pzAIn,
  double pxBIn, double pyBIn, double pzBIn) {
  if ( pxAIn == beamSetupPtr->pxA && pyAIn == beamSetupPtr->pyA &&
       pzAIn == beamSetupPtr->pzA && pxBIn == beamSetupPtr->pxB &&
       pyBIn == beamSetupPtr->pyB && pzBIn == beamSetupPtr->pzB )
    return true;

  if ( !beamSetupPtr->setKinematics(pxAIn, pyAIn, pzAIn,
                                    pxBIn, pyBIn, pzBIn) ) return false;
  return setKinematicsCM();
}

bool Angantyr::setKinematics(Vec4 pA, Vec4 pB) {
  return setKinematics(pA.px(), pA.py(), pA.pz(), pB.px(), pB.py(), pB.pz());
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
    // If beam energy is set to less than the mass, it is assumed at rest.
    if (bs.eA < bs.mA ||
        ( projPtr && projPtr->A() > 1 &&
          bs.eA <= particleDataPtr->m0(2112) ) ) {
      bs.pzA = 0.;
      bs.eA = bs.mA;
    }
    else {
      bs.pzA = sqrt(pow2(bs.eA) - pow2(bs.mA));
    }
    if ( bs.eB <= bs.mB ||
         ( targPtr && targPtr->A() > 1 &&
           bs.eB <= particleDataPtr->m0(2112) ) ) {
      bs.pzB = 0.;
      bs.eB = bs.mB;
    }
    else {
      bs.pzB = -sqrt(pow2(bs.eB) - pow2(bs.mB));
    }

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

  resetEvent();
  int itry = MAXTRY;

  while ( itry-- && !doAbort) {

    // Generate impact parameter, nuclei, and sub-collisions.
    double bweight = 0.0;
    Vec4 bvec = bGenPtr->generate(bweight);
    proj = Nucleus(projPtr->generate(), bvec / 2.);
    targ = Nucleus(targPtr->generate(), -bvec / 2.);

    subColls = collPtr->getCollisions(proj, targ);
    hiInfo.setSubCollisions(&subColls);
    hiInfo.addAttempt(subColls.T(), bvec.pT(), bvec.phi(),
                      bweight, bGenPtr->xSecScale());

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

    // Bunch all events together.
    if ( subEvents.empty() || !buildEvent(subEvents) ) {
      loggerPtr->ERROR_MSG("failed to build full event");
      continue;
    }

    // Hadronise everything, if requested.
    if (doHadronLevel) {
     if ( HIHooksPtr && HIHooksPtr->canForceHadronLevel() ) {
        if ( !HIHooksPtr->forceHadronLevel(*pythia[HADRON]) ) continue;
      } else {
        if ( !pythia[HADRON]->forceHadronLevel(false) ) continue;
      }
    }

    // Finally, boost to the requested frame and optionally do vertex
    // spreading.
    pythia[HADRON]->event.rotbst(
      fromCMframe(beamSetupPtr->pAnow, beamSetupPtr->pBnow));

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
      "parameter point", "consider reducing "
                         "HeavyIon:bWidth or HeavyIon:bWidthCut ");
  hiInfo.reject();
  return false;

}

//==========================================================================

} // end namespace Pythia8
