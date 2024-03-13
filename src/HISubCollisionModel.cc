// HISubCollisionModel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HISubCollisionModel.h header) for
// the built-in heavy ion subcollision models.

#include "Pythia8/Pythia.h"
#include "Pythia8/HISubCollisionModel.h"

namespace Pythia8 {

//==========================================================================

// ImpactParameterGenerator samples the impact parameter space.

//--------------------------------------------------------------------------

// Initialise base class, passing pointers to important objects.

void ImpactParameterGenerator::initPtr(Info & infoIn,
  SubCollisionModel & collIn, NucleusModel & projIn, NucleusModel & targIn) {
  infoPtr = &infoIn;
  settingsPtr = infoIn.settingsPtr;
  rndmPtr = infoIn.rndmPtr;
  loggerPtr = infoIn.loggerPtr;
  collPtr = &collIn;
  projPtr = &projIn;
  targPtr = &targIn;

}

//--------------------------------------------------------------------------

// Initialise base class, bay be overridden by subclasses.

bool ImpactParameterGenerator::init() {
  if ( settingsPtr->isParm("HI:bWidth") )
    widthSave = settingsPtr->parm("HI:bWidth")*femtometer;
  else
    widthSave = settingsPtr->parm("HeavyIon:bWidth")*femtometer;

  if ( widthSave <= 0.0 )
    updateWidth();

  return true;
}

//--------------------------------------------------------------------------

// Set width based on the associated subcollision and nucleus models.

void ImpactParameterGenerator::updateWidth() {
  double Rp = sqrt(collPtr->sigTot()/M_PI)/2.0;
  double RA = max(Rp, projPtr->R());
  double RB = max(Rp, targPtr->R());
  widthSave = RA + RB + 2.0*Rp;
}

//--------------------------------------------------------------------------

// Generate an impact parameter according to a gaussian distribution.

Vec4 ImpactParameterGenerator::generate(double & weight) const {
  double b = sqrt(-2.0*log(rndmPtr->flat()))*width();
  double phi = 2.0*M_PI*rndmPtr->flat();
  weight = 2.0*M_PI*width()*width()*exp(0.5*b*b/(width()*width()));
  return Vec4(b*sin(phi), b*cos(phi), 0.0, 0.0);
}

//==========================================================================

// The SubCollisionModel base class for modeling the collision between
// two nucleons to tell which type of collision has occurred. The
// model may manipulate the corresponing state of the nucleons.

//--------------------------------------------------------------------------

shared_ptr<SubCollisionModel> SubCollisionModel::create(int model) {
  switch (model) {
    case 0: return make_shared<NaiveSubCollisionModel>();
    case 1: return make_shared<DoubleStrikmanSubCollisionModel>();
    case 2: return make_shared<DoubleStrikmanSubCollisionModel>(1);
    case 3: return make_shared<BlackSubCollisionModel>();
    case 4: return make_shared<LogNormalSubCollisionModel>();
    case 5: return make_shared<LogNormalSubCollisionModel>(1);
    default: return nullptr;
  }
}

//--------------------------------------------------------------------------

// Initialize the base class. Subclasses should consider calling this
// in overriding functions.

bool SubCollisionModel::init(int idAIn, int idBIn, double eCMIn) {

  // Store input.
  idASave = idAIn;
  idBSave = idBIn;

  // Read basic settings.
  NInt = settingsPtr->mode("HeavyIon:SigFitNInt");
  NPop = settingsPtr->mode("HeavyIon:SigFitNPop");
  sigErr = settingsPtr->pvec("HeavyIon:SigFitErr");
  sigFuzz = settingsPtr->parm("HeavyIon:SigFitFuzz");
  fitPrint = settingsPtr->flag("HeavyIon:SigFitPrint");
  impactFudge = settingsPtr->parm("Angantyr:impactFudge");
  doVarECM = settingsPtr->flag("Beams:allowVariableEnergy");
  if (doVarECM) {
    eMin = settingsPtr->parm("HeavyIon:varECMMin");
    eMax = settingsPtr->parm("HeavyIon:varECMMax");
    if (eMax == 0)
      eMax = eCMIn;
    else if (eMax < eCMIn) {
      loggerPtr->ERROR_MSG("maximum energy is lower than requested eCM");
      return false;
    }
  }
  else
    eMin = eMax = eCMIn;
  updateSig();

  // If there are parameters, no further initialization is necessary.
  if (nParms() == 0) return true;

  // First try to load configuration from file, if requested.
  int    reuseInitMode = settingsPtr->mode("HeavyIon:SigFitReuseInit");
  string reuseInitFile = settingsPtr->word("HeavyIon:SigFitInitFile");
  bool   reuseWorked   = (reuseInitMode == 2 || reuseInitMode == 3)
    && loadParms(reuseInitFile);

  if (!reuseWorked) {
    if (reuseInitMode == 2) {
      loggerPtr->ABORT_MSG("unable to load parameter data");
      return false;
    }

    // If parameters were not loaded, generate from scratch.
    if (!genParms()) {
      loggerPtr->ABORT_MSG("evolutionary algorithm failed");
      return false;
    }
  }

  // Set parameters at the correct kinematics.
  setKinematics(eCMIn);

  // Set initial avNDb
  avNDb = getSig().avNDb * impactFudge;

  // Save parameters to disk, if requested.
  if (reuseInitMode == 1 || (reuseInitMode == 3 && !reuseWorked) ) {
    if (saveParms(reuseInitFile)) loggerPtr->INFO_MSG(
      "wrote initialization configuration to file", reuseInitFile);
    else loggerPtr->WARNING_MSG("couldn't save initialization configuration");
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Generate parameters based on run settings and the evolutionary algorithm.

bool SubCollisionModel::genParms() {

  // Initialize with default parameters.
  int nGen = settingsPtr->mode("HeavyIon:SigFitNGen");
  vector<double> defaultParms = settingsPtr->pvec("HeavyIon:SigFitDefPar");
  if ( settingsPtr->isPVec("HI:SigFitDefPar") )
    defaultParms = settingsPtr->pvec("HI:SigFitDefPar");
  if (defaultParms.size() == 0)
    defaultParms = defParm();
  if (int(defaultParms.size()) < nParms()) {
    loggerPtr->ERROR_MSG("too few parameters have been specified",
      "(expected " + to_string(nParms())
      + ", got " + to_string(defaultParms.size()) + ")");
    return false;
  }
  if (int(defaultParms.size()) > nParms()) {
    loggerPtr->WARNING_MSG("too many parameters have been specified",
      "(expected " + to_string(nParms())
      + ", got " + to_string(defaultParms.size()) + ")");
    defaultParms.resize(nParms());
  }
  setParm(defaultParms);

  // If nGen is zero, there is nothing to do, just use the default parameters.
  if (nGen == 0) {
    subCollParms = vector<LogInterpolator>(nParms());
    for (int iParm = 0; iParm < nParms(); ++iParm)
      subCollParms[iParm] = LogInterpolator(eMin, eMax, {defaultParms[iParm]});
    return true;
  }

  // Run evolutionary algorithm.
  if ( fitPrint ) {
    cout << " *------ HeavyIon fitting of SubCollisionModel to "
         << "cross sections ------* " << endl;
    flush(cout);
  }
  if (!evolve(nGen, eMax)) {
    loggerPtr->ERROR_MSG("evolutionary algorithm failed");
    return false;
  }
  defaultParms = getParm();

  // If we don't care about varECM, we are done.
  if (!doVarECM) {
    if (fitPrint) {
      cout << " *--- End HeavyIon fitting of parameters in "
        << "nucleon collision model ---* "
        << endl << endl;
      cout << " To avoid refitting, add the following lines to your "
             "configuration file: " << endl;
      cout << "  HeavyIon:SigFitNGen = 0" << endl;
      cout << "  HeavyIon:SigFitDefPar = ";
      for (int iParm = 0; iParm < nParms(); ++iParm) {
        if (iParm > 0) cout << ",";
        cout << defaultParms[iParm];
      }
      cout << endl << endl;
    }
    subCollParms = vector<LogInterpolator>(nParms());
    for (int iParm = 0; iParm < nParms(); ++iParm)
      subCollParms[iParm] = LogInterpolator(eMin, eMax, {defaultParms[iParm]});
    return true;
  }

  // Read settings for varECM evolution.
  eCMPts = settingsPtr->mode("HeavyIon:varECMSigFitNPts");
  bool doStepwiseEvolve = settingsPtr->flag("HeavyIon:varECMStepwiseEvolve");
  int nGenNext = settingsPtr->mode("HeavyIon:varECMSigFitNGen");

  // Vector of size nParms, each entry contains the parameter values.
  vector<vector<double>> parmsByECM(nParms(), vector<double>(eCMPts));

  // Write parameters at original eCM.
  for (int iParm = 0; iParm < nParms(); ++iParm)
    parmsByECM[iParm].back() = defaultParms[iParm];

  // Evolve down to eMin.
  vector<double> eCMs = logSpace(eCMPts, eMin, eMax);
  for (int i = eCMPts - 2; i >= 0; --i) {
    // Update to correct eCM.
    double eNow = eCMs[i];
    sigTotPtr->calc(idASave, idBSave, eNow);
    updateSig();

    // Alternatively reset to default parameters (mostly for debug purposes).
    if (!doStepwiseEvolve)
      setParm(defaultParms);

    // Evolve and get next set of parameters.
    if (fitPrint)
      cout << " *------------------------------------------"
             "---------------------------* "
           << endl;

    if (!evolve(nGenNext, eNow)) {
      loggerPtr->ERROR_MSG("evolutionary algorithm failed");
      return false;
    }
    vector<double> parmsNow = getParm();
    for (int iParm = 0; iParm < nParms(); ++iParm)
      parmsByECM[iParm][i] = parmsNow[iParm];
  }
  if (fitPrint){
    cout << " *--- End HeavyIon fitting of parameters in "
         << "nucleon collision model ---* "
         << endl << endl;
    cout << " To avoid refitting, you may use the HeavyIon:SigFitReuseInit"
            " parameter \n to store the configuration to disk."
         << endl << endl;
  }
  // Reset cross section and parameters to their eCM values.
  sigTotPtr->calc(idASave, idBSave, eMax);
  updateSig();
  setParm(defaultParms);

  // Store parameter values as logarithmic interpolators.
  subCollParms = vector<LogInterpolator>(nParms());
  for (int iParm = 0; iParm < nParms(); ++iParm)
    subCollParms[iParm] = LogInterpolator(eMin, eMax, parmsByECM[iParm]);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Save parameter configuration to disk.

bool SubCollisionModel::saveParms(string fileName) const {

  if (nParms() == 0) {
    loggerPtr->ERROR_MSG("model does not have any parameters");
    return true;
  }

  ofstream stream(fileName);
  if (!stream.good()) {
    loggerPtr->ERROR_MSG("unable to open file for writing", fileName);
    return false;
  }

  // Write energy range
  stream << subCollParms.front().data().size()
         << " " << eMin << " " << eMax << endl;

  // Each line corresponds to one parameter.
  for (int iParm = 0; iParm < nParms(); ++iParm) {
    stream << setprecision(14);
    for (double val : subCollParms[iParm].data())
      stream << val << " ";
    stream << endl;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Load parameter configuration from disk.

bool SubCollisionModel::loadParms(string fileName) {

  if (nParms() == 0) {
    loggerPtr->ERROR_MSG("model does not have any parameters");
    return true;
  }

  ifstream stream(fileName);
  if (!stream.good()) {
    loggerPtr->ERROR_MSG("unable to open file for reading", fileName);
    return false;
  }

  // Lambda function as a shorthand for error message.
  auto formatError = [this]() {
    loggerPtr->ERROR_MSG("invalid format");
    return false;
  };

  // Read first line
  string line;
  if (!getline(stream, line)) return formatError();
  if (!(istringstream(line) >> eCMPts >> eMin >> eMax)) return formatError();

  // Read each line and use the data to define an interpolator.
  subCollParms = vector<LogInterpolator>(nParms());
  vector<double> defaultParms(nParms());
  for (int iParm = 0; iParm < nParms(); ++iParm) {
    if (!getline(stream, line)) return formatError();

    istringstream lineStream(line);
    vector<double> parmData(eCMPts);
    for (int iPt = 0; iPt < eCMPts; ++iPt) {
      if (!(lineStream >> parmData[iPt]))
        return formatError();
    }

    subCollParms[iParm] = LogInterpolator(eMin, eMax, parmData);
    defaultParms[iParm] = parmData.back();
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Update the parameters to the interpolated value at the given eCM.

void SubCollisionModel::setKinematics(double eCMIn) {
  if (nParms() > 0) {
    vector<double> parmsNow(subCollParms.size());
    for (size_t iParm = 0; iParm < parmsNow.size(); ++iParm)
      parmsNow[iParm] = subCollParms[iParm](eCMIn);
    setParm(parmsNow);
    avNDb = getSig().avNDb * impactFudge;
  }
}

//--------------------------------------------------------------------------

// Update internally stored cross sections.

void SubCollisionModel::updateSig() {
  sigTarg[0] = sigTotPtr->sigmaTot()*millibarn;
  sigTarg[1] = sigTotPtr->sigmaND()*millibarn;
  sigTarg[2] = sigTotPtr->sigmaXX()*millibarn;
  sigTarg[3] = sigTotPtr->sigmaAX()*millibarn + sigTarg[1] + sigTarg[2];
  sigTarg[4] = sigTotPtr->sigmaXB()*millibarn + sigTarg[1] + sigTarg[2];
  sigTarg[5] = sigTotPtr->sigmaAXB()*millibarn;
  sigTarg[6] = sigTotPtr->sigmaEl()*millibarn;
  sigTarg[7] = sigTotPtr->bSlopeEl();
  // preliminarily set average ND impact parameter as if black disk.
  avNDb = 2.0 * sqrt(sigTarg[1]/M_PI) * impactFudge / 3.0;
}

//--------------------------------------------------------------------------

// Calculate the Chi^2 for the cross section that model in a subclass
// tries to model.

double SubCollisionModel::Chi2(const SigEst & se, int npar) const {

  double chi2 = 0.0;
  int nval = 0;
  for ( int i = 0, Nval = se.sig.size(); i < Nval; ++i ) {
    if ( sigErr[i] == 0.0 ) continue;
    ++nval;
    chi2 += pow2(se.sig[i] - sigTarg[i])/
      (se.dsig2[i] + pow2(sigTarg[i]*sigErr[i]));
  }
  return chi2/double(max(nval - npar, 1));
}


//--------------------------------------------------------------------------

// Anonymous helper function to print out stuff.

namespace {

void printFit(string name, double fit, double sig, double sigerr,
                 string unit = "mb    ") {
  cout << " |" << setw(25) << name << ": "
       << setw(8) << fit
       << (sigerr > 0.0? " *(": "  (")
       << setw(6) << sig;
  if ( sigerr > 0.0 )
    cout << " +- " << setw(2) << int(100.0*sigerr)  << "%";
  else
    cout << "       ";
  cout << ") " << unit << "          | " << endl;
}

}
//--------------------------------------------------------------------------

// A simple genetic algorithm for fitting the parameters in a subclass
// to reproduce desired cross sections.

bool SubCollisionModel::evolve(int nGenerations, double eCM) {

  if (nParms() == 0)
    return true;
  if (nGenerations <= 0)
    return true;

  if ( fitPrint ) {
    cout << " |                                      "
         << "                               | \n"
         << " |   Fitting parameters for " << setprecision(1) << setw(8)
         << eCM << " GeV                               | \n |   ";
    flush(cout);
  }

  // We're going to use a home-made genetic algorithm. We start by
  // creating a population of random parameter points.
  typedef vector<double> Parms;
  Parms minp = minParm();
  Parms maxp = maxParm();
  Parms defp = getParm();
  int dim = nParms();

  // Population of parameter sets. The most accurate sets will propagate
  // to the next generation.
  vector<Parms> pop(NPop, Parms(dim));
  for ( int j = 0; j < dim; ++j )
    pop[0][j] = clamp(defp[j], minp[j], maxp[j]);
  for ( int i = 1; i < NPop; ++i )
    for ( int j = 0; j < dim; ++j )
      pop[i][j] = minp[j] + rndmPtr->flat()*(maxp[j] - minp[j]);

  // Now we evolve our population for a number of generations.
  for ( int iGen = 0; iGen < nGenerations; ++iGen ) {

    // Calculate Chi2 for each parameter set and order them.
    multimap<double, Parms> chi2map;
    double chi2max = 0.0;
    for ( int i = 0; i < NPop; ++i ) {
      setParm(pop[i]);
      double chi2 = Chi2(getSig(), dim);
      chi2map.insert(make_pair(chi2, pop[i]));
      chi2max = max(chi2max, chi2);
    }

    if (fitPrint) {
      if ( iGen >= nGenerations - 20 ) cout << ".";
      flush(cout);
    }

    // Keep the best one, and move the other closer to a better one or
    // kill them if they are too bad.
    multimap<double, Parms>::iterator it = chi2map.begin();
    pop[0] = it->second;
    for ( int i = 1; i < NPop; ++i ) {
      ++it;
      double chi2Now = it->first;
      const Parms& parmsNow = it->second;
      pop[i] = it->second;
      if ( chi2Now > rndmPtr->flat()*chi2max ) {
        // Kill this individual and create a new one.
        for ( int j = 0; j < dim; ++j )
          pop[i][j] = minp[j] + rndmPtr->flat()*(maxp[j] - minp[j]);
      } else {
        // Pick one of the better parameter sets and move this closer.
        int ii = int(rndmPtr->flat()*i);
        for ( int j = 0; j < dim; ++j ) {
          double d = pop[ii][j] - parmsNow[j];
          double pl = clamp(parmsNow[j] - sigFuzz*d, minp[j], maxp[j]);
          double pu = clamp(parmsNow[j] + (1.0 + sigFuzz)*d, minp[j], maxp[j]);
          pop[i][j] = pl + rndmPtr->flat()*(pu - pl);
        }
      }
    }
  }

  // Update resulting parameter set.
  setParm(pop[0]);
  SigEst se = getSig();

  // Output information.
  double chi2 = Chi2(se, dim);
  avNDb = se.avNDb*impactFudge;
  if ( fitPrint ) {
    for ( int i = nGenerations; i < 20; ++i ) cout << " ";
    cout << "                                              | \n";
    cout << " |                                      "
         << "                               | "
         << endl;
    if ( nGenerations > 0 ) {
      cout << " |     Resulting parameters:         "
           << "                                  | "
           << endl;
      for (int iParm = 0; iParm < this->nParms(); ++iParm) {
        cout << " |" << setw(25) << "[" + to_string(iParm) << "]: "
             << setprecision(2) << setw(7) << pop[0][iParm]
             << setw(36) << "| " << endl;
      }
      cout << " |                                      "
           << "                               | "
           << endl;
    }
    cout << " |     Resulting cross sections        (target value) "
         << "                 | "
         << endl;
    printFit("Total", se.sig[0]/millibarn,
             sigTarg[0]/millibarn, sigErr[0]);
    printFit("non-diffractive", se.sig[1]/millibarn,
             sigTarg[1]/millibarn, sigErr[1]);
    printFit("XX diffractive", se.sig[2]/millibarn,
             sigTarg[2]/millibarn, sigErr[2]);
    printFit("wounded target (B)", se.sig[3]/millibarn,
             sigTarg[3]/millibarn, sigErr[3]);
    printFit("wounded projectile (A)", se.sig[4]/millibarn,
             sigTarg[4]/millibarn, sigErr[4]);
    printFit("AXB diffractive", se.sig[5]/millibarn,
             sigTarg[5]/millibarn, sigErr[5]);
    printFit("elastic", se.sig[6]/millibarn,
             sigTarg[6]/millibarn, sigErr[6]);
    printFit("elastic b-slope", se.sig[7], sigTarg[7], sigErr[7], "GeV^-2");
    cout << " |                                   "
         << "                                  | "
         << endl;
    cout << " |" << setw(25) << "Chi2/Ndf" << ": ";
    cout << fixed << setprecision(2);
    cout << setw(8) << chi2 << "                                  | \n";

    cout << " |                                      "
         << "                               | "
         << endl;
  }

  // Done.
  return true;

}

//==========================================================================

// The BlackSubCollisionModel uses fixed size, black-disk
// nucleon-nucleon cross section, equal to the total inelastic pp cross
// section. Everything else is elastic -- Diffraction not included.

//--------------------------------------------------------------------------

SubCollisionSet BlackSubCollisionModel::
getCollisions(Nucleus& proj, Nucleus& targ) {

  multiset<SubCollision> ret;

  // Go through all pairs of nucleons
  for (Nucleon& p : proj)
    for (Nucleon& t : targ) {
      double b = (p.bPos() - t.bPos()).pT();
      if ( b > sqrt(sigTot()/M_PI) ) continue;
      if ( b < sqrt((sigTot() - sigEl())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ABS));
      }
      else {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ELASTIC));
      }
    }

  return SubCollisionSet(ret, 0.5);
}

//==========================================================================

// The NaiveSubCollisionModel uses a fixed size, black-disk-like
// nucleon-nucleon cross section where. Central collisions will always
// be absorptive, less central will be doubly diffractive, more
// peripheral will be single diffractive and the most peripheral will
// be elastic.

//--------------------------------------------------------------------------

SubCollisionSet NaiveSubCollisionModel::
getCollisions(Nucleus& proj, Nucleus& targ) {

  multiset<SubCollision> ret;

  // Go through all pairs of nucleons
  for (Nucleon& p : proj)
    for (Nucleon& t : targ) {
      double b = (p.bPos() - t.bPos()).pT();
      if ( b > sqrt(sigTot()/M_PI) ) continue;
      if ( b < sqrt(sigND()/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ABS));
      }
      else if ( b < sqrt((sigND() + sigDDE())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::DDE));
      }
      else if ( b < sqrt((sigND() + sigSDE() + sigDDE())/M_PI) ) {
         if ( sigSDEP() > rndmPtr->flat()*sigSDE() ) {
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDEP));
        } else {
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDET));
        }
      }
      else if ( b < sqrt((sigND() + sigSDE() + sigDDE() + sigCDE())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::CDE));
      }
      else {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ELASTIC));
      }
    }

  return SubCollisionSet(ret, 0.5);
}

//==========================================================================

// DoubleStrikman uses a fluctuating and semi-transparent disk for
// each Nucleon in a sub-collision resulting in a fluctuating
// interaction probability. To assess the fluctuation each Nucleon has
// two random states in each collision, one main state and one helper
// state to assess the frlutuations.

//--------------------------------------------------------------------------

// Helper functions to get the correct average elastic and wounded
// cross sections for fluctuating models.

static void shuffle(double PND1, double PND2, double & PW1, double & PW2) {
  if ( PND1 > PW1 ) {
    PW2 += PW1 - PND1;
    PW1 = PND1;
    return;
  }
  if ( PND2 > PW2 ) {
    PW1 += PW2 - PND2;
    PW2 = PND2;
    return;
  }
}

static void shuffle(double & PEL11, double P11,
                    double P12, double P21, double P22) {
  double PEL12 = PEL11, PEL21 = PEL11, PEL22 = PEL11;
  map<double, double *> ord;
  ord[P11] = &PEL11;
  ord[P12] = &PEL12;
  ord[P21] = &PEL21;
  ord[P22] = &PEL22;
  map<double, double *>::iterator next = ord.begin();
  map<double, double *>::iterator prev = next++;
  while ( next != ord.end() ) {
    if ( *prev->second > prev->first ) {
      *next->second += *prev->second - prev->first;
      *prev->second = prev->first;
    }
    prev = next++;
  }
}

static double pnw(double PWp, double PWt, double PND) {
  return ( 1.0 - PWp <= 0.0 || 1.0 - PWt <= 0.0 )?
    0.0: (1.0 - PWp)*(1.0 - PWt)/(1.0 - PND);
}

static double el(double s1, double s2, double u1, double u2) {
  return s1/u1 > s2/u2? s2*u1: s1*u2;
}

//--------------------------------------------------------------------------

// Numerically estimate the cross sections corresponding to the
// current parameter setting.

SubCollisionModel::SigEst FluctuatingSubCollisionModel::getSig() const {

  SigEst s;
  for ( int n = 0; n < NInt; ++n ) {
    double rp1 = pickRadiusProj();
    double rp2 = pickRadiusProj();
    double rt1 = pickRadiusTarg();
    double rt2 = pickRadiusTarg();
    double s11 = pow2(rp1 + rt1)*M_PI;
    double s12 = pow2(rp1 + rt2)*M_PI;
    double s21 = pow2(rp2 + rt1)*M_PI;
    double s22 = pow2(rp2 + rt2)*M_PI;

    double stot = (s11 + s12 + s21 + s22)/4.0;
    s.sig[0] += stot;
    s.dsig2[0] += pow2(stot);

    double u11 = opacity(s11)/2.0;
    double u12 = opacity(s12)/2.0;
    double u21 = opacity(s21)/2.0;
    double u22 = opacity(s22)/2.0;

    double avb = sqrt(2.0/M_PI)*(s11*sqrt(s11/(2.0*u11))*(1.0 - u11) +
                                 s12*sqrt(s12/(2.0*u12))*(1.0 - u12) +
                                 s21*sqrt(s21/(2.0*u21))*(1.0 - u21) +
                                 s22*sqrt(s22/(2.0*u22))*(1.0 - u22))/12.0;
    s.avNDb += avb;
    s.davNDb2 += pow2(avb);

    double snd = (s11 - s11*u11 + s12 - s12*u12 +
                  s21 - s21*u21 + s22 - s22*u22)/4.0;
    s.sig[1] += snd;
    s.dsig2[1] += pow2(snd);

    double sel = (el(s11, s22, u11, u22) + el(s12, s21, u12, u21))/2.0;
    s.sig[6] += sel;
    s.dsig2[6] += pow2(sel);

    double swt = stot - (el(s11, s12, u11, u12) + el(s21, s22, u21, u22))/2.0;
    double swp = stot - (el(s11, s21, u11, u21) + el(s12, s22, u12, u22))/2.0;
    s.sig[4] += swp;
    s.dsig2[4] += pow2(swp);
    s.sig[3] += swt;
    s.dsig2[3] += pow2(swt);

    s.sig[2] += swt + swp - snd  + sel - stot;
    s.dsig2[2] += pow2(swt + swp - snd  + sel - stot);

    s.sig[5] += s11;
    s.dsig2[5] += pow2(s11);

    s.sig[7] += pow2(s11)/u11;
    s.dsig2[7] += pow2(pow2(s11)/u11);

  }

  s.sig[0] /= double(NInt);
  s.dsig2[0] = (s.dsig2[0]/double(NInt) - pow2(s.sig[0]))/double(NInt);

  s.sig[1] /= double(NInt);
  s.dsig2[1] = (s.dsig2[1]/double(NInt) - pow2(s.sig[1]))/double(NInt);

  s.sig[2] /= double(NInt);
  s.dsig2[2] = (s.dsig2[2]/double(NInt) - pow2(s.sig[2]))/double(NInt);

  s.sig[3] /= double(NInt);
  s.dsig2[3] = (s.dsig2[3]/double(NInt) - pow2(s.sig[3]))/double(NInt);

  s.sig[4] /= double(NInt);
  s.dsig2[4] = (s.dsig2[4]/double(NInt) - pow2(s.sig[4]))/double(NInt);

  s.sig[6] /= double(NInt);
  s.dsig2[6] = (s.dsig2[6]/double(NInt) - pow2(s.sig[6]))/double(NInt);

  s.sig[5] /= double(NInt);
  s.dsig2[5] /= double(NInt);

  s.sig[7] /= double(NInt);
  s.dsig2[7] /= double(NInt);
  double bS = (s.sig[7]/s.sig[5])/(16.0*M_PI*pow2(0.19732697));
  double b2S = pow2(bS)*(s.dsig2[7]/pow2(s.sig[7]) - 1.0 +
                        s.dsig2[5]/pow2(s.sig[5]) - 1.0)/double(NInt);
  s.sig[5] = 0.0;
  s.dsig2[5] = 0.0;
  s.sig[7] = bS;
  s.dsig2[7] = b2S;

  s.avNDb /= double(NInt);
  s.davNDb2 = (s.davNDb2/double(NInt) - pow2(s.avNDb))/double(NInt);
  s.avNDb /= s.sig[1];
  s.davNDb2 /= pow2(s.sig[1]);

  return s;

}

//--------------------------------------------------------------------------

// Main function returning the possible sub-collisions.

SubCollisionSet FluctuatingSubCollisionModel::
getCollisions(Nucleus& proj, Nucleus& targ) {

  multiset<SubCollision> ret;

  // Assign two states to each nucleon.
  for (Nucleon& p : proj) {
    p.state({ pickRadiusProj() });
    p.addAltState({ pickRadiusProj() });
  }
  for (Nucleon& t : targ) {
    t.state({ pickRadiusTarg() });
    t.addAltState({ pickRadiusTarg() });
  }

  // The factorising S-matrix.
  double S = 1.0;

  // Go through all pairs of nucleons
  for (Nucleon& p : proj)
    for (Nucleon& t : targ) {
      double b = (p.bPos() - t.bPos()).pT();

      double T11 = Tpt(p.state(), t.state(), b);
      double T12 = Tpt(p.state(), t.altState(), b);
      double T21 = Tpt(p.altState(), t.state(), b);
      double T22 = Tpt(p.altState(), t.altState(), b);
      double S11 = 1.0 - T11;
      double S12 = 1.0 - T12;
      double S21 = 1.0 - T21;
      double S22 = 1.0 - T22;
      S *= S11;
      double PND11 = 1.0 - pow2(S11);
      // First and most important, check if this is an absorptive
      // scattering.
      if ( PND11 > rndmPtr->flat() ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ABS));
        continue;
      }

      // Now set up calculation for probability of diffractively
      // wounded nucleons.
      double PND12 = 1.0 - pow2(S12);
      double PND21 = 1.0 - pow2(S21);
      double PWp11 = 1.0 - S11*S21;
      double PWp21 = 1.0 - S11*S21;
      shuffle(PND11, PND21, PWp11, PWp21);
      double PWt11 = 1.0 - S11*S12;
      double PWt12 = 1.0 - S11*S12;
      shuffle(PND11, PND12, PWt11, PWt12);

      bool wt = ( PWt11 - PND11 > (1.0 - PND11)*rndmPtr->flat() );
      bool wp = ( PWp11 - PND11 > (1.0 - PND11)*rndmPtr->flat() );
      if ( wt && wp ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::DDE));
        continue;
      }
      if ( wt ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDET));
        continue;
      }
      if ( wp ) {
        ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::SDEP));
        continue;
      }

      // Finally set up calculation for elastic scattering. This can
      // never be exact, but let's do as well as we can.

      double PND22 = 1.0 - pow2(S22);
      double PWp12 = 1.0 - S12*S22;
      double PWp22 = 1.0 - S12*S22;
      shuffle(PND12, PND22, PWp12, PWp22);
      double PWt21 = 1.0 - S21*S22;
      double PWt22 = 1.0 - S21*S22;
      shuffle(PND21, PND22, PWt21, PWt22);

      double PNW11 = pnw(PWp11, PWt11, PND11);
      double PNW12 = pnw(PWp12, PWt12, PND12);
      double PNW21 = pnw(PWp21, PWt21, PND21);
      double PNW22 = pnw(PWp22, PWt22, PND22);

      double PEL = (T12*T21 + T11*T22)/2.0;
      shuffle(PEL, PNW11, PNW12, PNW21, PNW22);
      if ( PEL > PNW11*rndmPtr->flat() ) {
        if ( sigCDE() > rndmPtr->flat()*(sigCDE() + sigEl()) )
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::CDE));
        else
          ret.insert(SubCollision(p, t, b, b/avNDb, SubCollision::ELASTIC));
      }
    }

  return SubCollisionSet(ret, 1.0 - S);
}

//==========================================================================

} // end namespace Pythia8
