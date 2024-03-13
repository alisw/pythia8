// HINucleusModel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HINucleusModel.h header) for
// the built-in heavy ion nucleus models.

#include "Pythia8/Pythia.h"
#include "Pythia8/HINucleusModel.h"

namespace Pythia8 {

//==========================================================================

// The Nucleon class describing a Nucleon in a Nucleus.

//--------------------------------------------------------------------------

// Print out information about a Nucleon. To be called from inside a
// debugger.

void Nucleon::debug() {
  cout << "Nucleon id: " << id() << endl;
  cout << "index:      " << index() << endl;
  cout << "b(rel):     " << nPos().px() << " " << nPos().py() << endl;
  cout << "b(abs):     " << bPos().px() << " " << bPos().py() << endl;
  cout << "status:     " << status() << (done()? " done": "     ") << endl;
  cout << "state:      ";
  for ( int i = 0, N = state().size(); i < N; ++i )
    cout << state()[i] << " ";
  cout << endl;
  for ( int j = 0, M = altStatesSave.size(); j < M; ++j ) {
    cout << "state " << j+1 << ":    ";
    for ( int i = 0, N = altStatesSave[j].size(); i < N; ++i )
      cout << altStatesSave[j][i] << " ";
    cout << endl;
  }
}

//==========================================================================

// NucleusModel base class.

//--------------------------------------------------------------------------

shared_ptr<NucleusModel> NucleusModel::create(int model) {

  switch (model) {
    case 1: return make_shared<GLISSANDOModel>();
    case 2: return make_shared<WoodsSaxonModel>();
    case 3: return make_shared<HOShellModel>();
    case 4: return make_shared<GaussianModel>();
    case 5: return make_shared<HulthenModel>();
    default: return nullptr;
  }
}

//--------------------------------------------------------------------------

// Initialise base class, passing pointers to important objects.

void NucleusModel::initPtr(int idIn, bool isProjIn, Info& infoIn) {
  isProj = isProjIn;
  idSave = idIn;
  infoPtr = &infoIn;
  settingsPtr = infoIn.settingsPtr;
  loggerPtr = infoIn.loggerPtr;
  rndmPtr = infoIn.rndmPtr;
  mSave  = infoIn.particleDataPtr->m0(idSave);
  int decomp = abs(idSave);
  ISave = decomp%10;
  decomp /= 10;
  ASave = decomp%1000;
  decomp /= 1000;
  ZSave = decomp%1000;
  decomp /= 1000;
  LSave = decomp%10;
  decomp /= 10;

  if ( decomp != 10 ) {
    LSave = 0;
    ISave = 0;
    ASave = 0;
    ZSave = 0;
  }
}

//--------------------------------------------------------------------------

// Produce a proper particle corresponding to the nucleus handled by
// this NucleusModel.

Particle NucleusModel::produceIon() {

  Particle p(id(), -12);
  p.p(fromCMframe(pNSave)*Vec4(0.0, 0.0, 0.0, mSave));
  p.m(mSave);
  p.daughter1(isProj ? 3 : 4);
  return p;
}

//==========================================================================

// ExternalNucleusModel reads in the configuration from a file.

//--------------------------------------------------------------------------

// Initialize, read in the file, shuffle the configurations.

bool ExternalNucleusModel::init() {
  // Read settings.

  // File to read.
  fName = isProj ? settingsPtr->word("HeavyIonA:NucleusFile")
                 : settingsPtr->word("HeavyIonB:NucleusFile");

  // Do shuffling.
  doShuffle = isProj ? settingsPtr->flag("HeavyIonA:Shuffle")
                     : settingsPtr->flag("HeavyIonB:Shuffle");

  // Read the file.
  ifstream ifs(fName);
  if (!ifs.is_open()) {
    loggerPtr->ABORT_MSG("could not open file", fName);
    return false;
  }
  for (string line; getline(ifs, line); ) {
    // Remove comment lines.
    if (line.find("#") != std::string::npos) continue;
    stringstream ss(line);
    vector<double> positions;
    for (double coord; ss >> coord; )
      positions.push_back(coord);
    // Check malformed line.
    if (positions.size() != size_t(3 * A())) {
      loggerPtr->ABORT_MSG(
        "number of entries on each line must be 3 x A", fName);
      return false;
    }
    else {
      vector<Vec4> tmp;
      for (int i = 0; i < A(); ++i)
        tmp.push_back(Vec4(positions[3 * i],
          positions[3 * i + 1], positions[3 * i + 2]));
      nucleonPositions.push_back(tmp);
    }
  }
  ifs.close();

  if (nucleonPositions.size() == 0) {
    loggerPtr->ABORT_MSG("no entries found");
    return false;
  }

  if (doShuffle) rndmPtr->shuffle(nucleonPositions);
  return NucleusModel::init();
}

//--------------------------------------------------------------------------

// Generate a vector of nucleons.

vector<Nucleon> ExternalNucleusModel::generate() const {
  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;

  // Get the next nucleon
  vector<Vec4> nucleus = nucleonPositions[nUsed];
  int Np = Z();
  int Nn = A() - Z();
  for (int i = 0; i < A(); ++i) {
    if ( int(rndmPtr->flat()*(Np + Nn)) >= Np ) {
      --Nn;
      nucleons[i] = Nucleon(nid, i, nucleus[i]);
    } else {
      --Np;
      nucleons[i] = Nucleon(pid, i, nucleus[i]);
    }
  }

  // Update position in list of nucleons.
  if (++nUsed == nucleonPositions.size()) {
    nUsed = 0;
    if (doShuffle) rndmPtr->shuffle(nucleonPositions);
  }

  return nucleons;
}


//==========================================================================

// HardCoreModel is a base class for models implementing a hard core.

//--------------------------------------------------------------------------

// Init the hard core parameters. To be called in init() in derived classes.
void HardCoreModel::initHardCore() {
  useHardCore = (isProj ? settingsPtr->flag("HeavyIonA:HardCore")
                        : settingsPtr->flag("HeavyIonB:HardCore"));
  hardCoreRadius = (isProj ? settingsPtr->parm("HeavyIonA:HardCoreRadius")
                           : settingsPtr->parm("HeavyIonB:HardCoreRadius"));
  gaussHardCore = (isProj ? settingsPtr->flag("HeavyIonA:GaussHardCore")
                          : settingsPtr->flag("HeavyIonB:GaussHardCore"));
}

//==========================================================================

// WoodsSaxonModel is a subclass of NucleusModel and implements a
// general Wood-Saxon distributed nucleus.

//--------------------------------------------------------------------------

// Initialize.
bool WoodsSaxonModel::init() {
  if (A() == 0) return true;

  // Initialize hard core.
  initHardCore();

  // In the basic Woods-Saxon model we get parameters directly from settings.
  RSave = settingsPtr->parm(isProj ? "HeavyIonA:WSR" : "HeavyIonB:WSR");
  aSave = settingsPtr->parm(isProj ? "HeavyIonA:WSa" : "HeavyIonB:WSa");

  // Calculate the overestimates.
  overestimates();
  return NucleusModel::init();
}

// Place a nucleon inside a nucleus.
Vec4 WoodsSaxonModel::generateNucleon() const {

  while ( true ) {
    double r = R();
    double sel = rndmPtr->flat()*(intlo + inthi0 + inthi1 + inthi2);
    if ( sel > intlo ) r -= a()*log(rndmPtr->flat());
    if ( sel > intlo + inthi0 ) r -= a()*log(rndmPtr->flat());
    if ( sel > intlo + inthi0 + inthi1 )  r -= a()*log(rndmPtr->flat());
    if ( sel <= intlo ) {
      r = R()*pow(rndmPtr->flat(), 1.0/3.0);
      if ( rndmPtr->flat()*(1.0 + exp((r - R())/a())) > 1.0 ) continue;
    } else
      if ( rndmPtr->flat()*(1.0 + exp((r - R())/a())) > exp((r - R())/a()) )
        continue;

    double costhe = 2.0*rndmPtr->flat() - 1.0;
    double sinthe = sqrt(max(1.0 - costhe*costhe, 0.0));
    double phi = 2.0*M_PI*rndmPtr->flat();

    return Vec4(r*sinthe*cos(phi), r*sinthe*sin(phi), r*costhe);
  }
  return Vec4();
}

//--------------------------------------------------------------------------

// Generate all nucleons in a nucleus.

vector<Nucleon> WoodsSaxonModel::generate() const {
  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;

  if ( A() == 0 ) {
    nucleons.push_back(Nucleon(id(), 0, Vec4()));
    return nucleons;
  }
  if ( A() == 1 ) {
    if ( Z() == 1 ) nucleons.push_back(Nucleon(pid, 0, Vec4()));
    else  nucleons.push_back(Nucleon(nid, 0, Vec4()));
    return nucleons;
  }

  Vec4 cms;
  vector<Vec4> positions;
  while ( int(positions.size()) < A() ) {
    while ( true ) {
      Vec4 pos = generateNucleon();
      bool overlap = false;
      if (useHardCore) {
        for (int i = 0, N = positions.size(); i < N && !overlap; ++i ) {
          if ((positions[i] - pos).pAbs() < rSample() ) overlap = true;
        }
      }
      if ( overlap ) continue;
      positions.push_back(pos);
      cms += pos;
      break;
    }
  }

  cms /= A();
  nucleons.resize(A());
  int Np = Z();
  int Nn = A() - Z();
  for ( int i = 0, N= positions.size(); i < N; ++i ) {
    Vec4 pos(positions[i].px() - cms.px(),
                 positions[i].py() - cms.py());
    if ( int(rndmPtr->flat()*(Np + Nn)) >= Np ) {
      --Nn;
      nucleons[i] = Nucleon(nid, i, pos);
    } else {
      --Np;
      nucleons[i] = Nucleon(pid, i, pos);
    }
  }
  return nucleons;
}

//==========================================================================

// GLISSANDOModel is a special case of the WoodsSaxon model with specific
// parameters. The nuclear radius is calculated from A().

//--------------------------------------------------------------------------

// Initialize parameters.

bool GLISSANDOModel::init() {
  if ( A() == 0 ) return true;

  // Initialize hard core.
  initHardCore();

  // There are no parameters to be read.
  if (useHardCore) {
    RSave = (1.1*pow(double(A()),1.0/3.0) -
             0.656*pow(double(A()),-1.0/3.0))*femtometer;
    aSave = 0.459*femtometer;
  } else {
    RSave = (1.12*pow(double(A()),1.0/3.0) -
             0.86*pow(double(A()),-1.0/3.0))*femtometer;
    aSave = 0.54*femtometer;
  }

  // Calculate overestimates.
  overestimates();
  return NucleusModel::init();

}

//==========================================================================

// HOShellModel is derived from of NucleusModel and implements a
// harmonic oscillator shell model with option for hard cores. Suitable for
// nuclei with 4 <= A <= 16.

//--------------------------------------------------------------------------

// Initialize parameters.

bool HOShellModel::init() {
  if ( A() == 0 ) return true;

  // Initialize hard core parameters in base class.
  initHardCore();

  // Proton charge radius
  protonChR = isProj ? settingsPtr->parm("HeavyIonA:HOProtonChargeRadius")
                     : settingsPtr->parm("HeavyIonB:HOProtonChargeRadius");
  // Possible custom nuclear charge radius
  nucleusChR = isProj ? settingsPtr->parm("HeavyIonA:HONuclearChargeRadius")
                      : settingsPtr->parm("HeavyIonB:HONuclearChargeRadius");

  // If not changed, use defaults for a set of nuclei.
  if (nucleusChR == 0.) {
    // Helium-4
    if (A() == 4 && Z() == 2) nucleusChR = (useHardCore ? 2.45 : 2.81);
    // Lithium-6
    else if (A() == 6 && Z() == 3) nucleusChR = (useHardCore ? 6.4 : 6.7);
    // Berylium-7
    else if (A() == 7 && Z() == 4) nucleusChR = (useHardCore ? 6.69 : 7.00);
    // Lithium-8
    else if (A() == 8 && Z() == 3) nucleusChR = (useHardCore ? 5.1 : 5.47);
    // Berylium-9
    else if (A() == 9 && Z() == 4) nucleusChR = (useHardCore ? 6.0 : 6.35);
    // Boron-10
    else if (A() == 10 && Z() == 5) nucleusChR = (useHardCore ? 5.5 : 5.89);
    // Boron-11
    else if (A() == 11 && Z() == 5) nucleusChR = (useHardCore ? 5.36 : 5.79);
    // Carbon-12
    else if (A() == 12 && Z() == 6) nucleusChR = (useHardCore ? 5.66 : 6.10);
    // Carbon-13
    else if (A() == 13 && Z() == 6) nucleusChR = (useHardCore ? 5.6 : 6.06);
    // Nitrogen-14
    else if (A() == 14 && Z() == 7) nucleusChR = (useHardCore ? 6.08 : 6.54);
    // Nitrogen-15
    else if (A() == 15 && Z() == 7) nucleusChR = (useHardCore ? 6.32 : 6.79);
    // Oxygen-16
    else if (A() == 16 && Z() == 8) nucleusChR = (useHardCore ? 6.81 : 7.29);
    else {
      loggerPtr->ERROR_MSG(
        "default parameters are not defined for this nucleus",
        "(with id=" + to_string(idSave) + ")");
      return false;
    }
  }
  // Calculate C2 prefactor.
  C2 = 1./(5./2. - 4./A()) * (nucleusChR - protonChR);
  rhoMax = A() < 10 ? rho(0) : rho( (sqrt((A() - 10)*sqrt(C2))/sqrt(A() - 4)));
  NucleusModel::init();
  return true;
}

//--------------------------------------------------------------------------

// Generate the position of a single nucleon.

Vec4 HOShellModel::generateNucleon() const {
  double r = -1;
  do {
    r = -C2 * log(rndmPtr->flat());
  } while (rndmPtr->flat() * 14./8. * rhoMax * exp(-r/C2) > rho(r) );

  double costhe = 2.0*rndmPtr->flat() - 1.0;
  double sinthe = sqrt(max(1.0 - costhe*costhe, 0.0));
  double phi = 2.0*M_PI*rndmPtr->flat();

  return Vec4(r*sinthe*cos(phi), r*sinthe*sin(phi), r*costhe);
}

//--------------------------------------------------------------------------

// Generate all the nucleons.

vector<Nucleon> HOShellModel::generate() const {
  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;

  if ( A() == 0 ) {
    nucleons.push_back(Nucleon(id(), 0, Vec4()));
    return nucleons;
  }
  if ( A() == 1 ) {
    if ( Z() == 1 ) nucleons.push_back(Nucleon(pid, 0, Vec4()));
    else  nucleons.push_back(Nucleon(nid, 0, Vec4()));
    return nucleons;
  }

  Vec4 cms;
  vector<Vec4> positions;
  while ( int(positions.size()) < A() ) {
    while ( true ) {
      Vec4 pos = generateNucleon();
      bool overlap = false;
      if (useHardCore) {
        for ( int i = 0, N = positions.size(); i < N && !overlap; ++i ) {
          if ( (positions[i] - pos).pAbs() < rSample() ) overlap = true;
              }
      }
      if ( overlap ) continue;
      positions.push_back(pos);
      cms += pos;
      break;
    }
  }

  cms /= A();
  nucleons.resize(A());
  int Np = Z();
  int Nn = A() - Z();
  for ( int i = 0, N= positions.size(); i < N; ++i ) {
    Vec4 pos(positions[i].px() - cms.px(),
                 positions[i].py() - cms.py());
    if ( int(rndmPtr->flat()*(Np + Nn)) >= Np ) {
      --Nn;
      nucleons[i] = Nucleon(nid, i, pos);
    } else {
      --Np;
      nucleons[i] = Nucleon(pid, i, pos);
    }
  }

  return nucleons;
}

//==========================================================================

// The Hulthen model for deuterons.

//--------------------------------------------------------------------------

// Initialize parameters.

bool HulthenModel::init() {
  if (! (A() == 2 && Z() == 1)) {
    loggerPtr->ABORT_MSG(
      "the Hulthen distribution is only valid for deuterons");
    return false;
  }
  // Note: Hulthen model has no hard core option.

  hA = isProj ? settingsPtr->parm("HeavyIonA:HulthenA")
              : settingsPtr->parm("HeavyIonB:HulthenA");

  hB = isProj ? settingsPtr->parm("HeavyIonA:HulthenB")
              : settingsPtr->parm("HeavyIonB:HulthenB");

  // Test that b > a.
  if (hB < hA) {
    loggerPtr->ABORT_MSG(
      "you must have HeavyIonX:HulthenB > HeavyIonX:HulthenA");
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------

// Generate all the nucleons.

vector<Nucleon> HulthenModel::generate() const {
  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;

  Vec4 cms;

  // Put one at (0,0,0).
  Vec4 posA(0, 0, 0);

  // Find the distance between the nucleons.
  double r;
  do {
    r = -hB * log(1. - rndmPtr->flat())/2./hA;
  } while (rndmPtr->flat() * exp(-2.*hA*r/hB) > rho(r));
  // Add the other one on a sphere around the first one.
  double costhe = 2.0*rndmPtr->flat() - 1.0;
  double sinthe = sqrt(max(1.0 - costhe*costhe, 0.0));
  double phi = 2.0*M_PI*rndmPtr->flat();

  Vec4 posB(r*sinthe*cos(phi), r*sinthe*sin(phi), r*costhe);
  cms += posB;
  cms /= A();
  nucleons.resize(A());

  // Add them to the vector.
  bool nFirst = (rndmPtr->flat() < 0.5);
  nucleons[0] = Nucleon((nFirst ? pid : nid), 0, Vec4(posA.px() - cms.px(),
      posA.py() - cms.py()));
    nucleons[1] = Nucleon((nFirst ? nid : pid), 0, Vec4(posB.px() - cms.px(),
      posB.py() - cms.py()));
  return nucleons;
}

//==========================================================================

// A Gaussian distribution for light nuclei.

//--------------------------------------------------------------------------

// Initialize parameters.

bool GaussianModel::init() {

  if ( A() == 0 ) return true;

  // Initialize hard core.
  initHardCore();

  // A single parameter.
  nucleusChR = isProj ? settingsPtr->parm("HeavyIonA:GaussianChargeRadius")
                      : settingsPtr->parm("HeavyIonB:GaussianChargeRadius");

  return true;
}

//--------------------------------------------------------------------------

// Generate the position of a single nucleon.

Vec4 GaussianModel::generateNucleon() const {
  double r;
  do {
    r = nucleusChR * rndmPtr->gauss();
  }
  while(r > 4 * nucleusChR);
  double costhe = 2.0*rndmPtr->flat() - 1.0;
  double sinthe = sqrt(max(1.0 - costhe*costhe, 0.0));
  double phi = 2.0*M_PI*rndmPtr->flat();

  return Vec4(r*sinthe*cos(phi), r*sinthe*sin(phi), r*costhe);
}

//--------------------------------------------------------------------------

// Generate all the nucleons.

vector<Nucleon> GaussianModel::generate() const {
  int sign = id() > 0? 1: -1;
  int pid = sign*2212;
  int nid = sign*2112;
  vector<Nucleon> nucleons;

  if ( A() == 0 ) {
    nucleons.push_back(Nucleon(id(), 0, Vec4()));
    return nucleons;
  }
  if ( A() == 1 ) {
    if ( Z() == 1 ) nucleons.push_back(Nucleon(pid, 0, Vec4()));
    else  nucleons.push_back(Nucleon(nid, 0, Vec4()));
    return nucleons;
  }

  Vec4 cms;
  vector<Vec4> positions;
  while ( int(positions.size()) < A() ) {
    while ( true ) {
      Vec4 pos = generateNucleon();
      bool overlap = false;
      if (useHardCore) {
        for ( int i = 0, N = positions.size(); i < N && !overlap; ++i )
          if ( (positions[i] - pos).pAbs() < abs(rndmPtr->gauss())
             * hardCoreRadius )
            overlap = true;
      }
      if ( overlap ) continue;
      positions.push_back(pos);
      cms += pos;
      break;
    }
  }

  cms /= A();
  nucleons.resize(A());
  int Np = Z();
  int Nn = A() - Z();
  for ( int i = 0, N = positions.size(); i < N; ++i ) {
    Vec4 pos(positions[i].px() - cms.px(), positions[i].py() - cms.py());
    if ( int(rndmPtr->flat()*(Np + Nn)) >= Np ) {
      --Nn;
      nucleons[i] = Nucleon(nid, i, pos);
    } else {
      --Np;
      nucleons[i] = Nucleon(pid, i, pos);
    }
  }

  return nucleons;
}

//==========================================================================

// ClusterModel generates nucleons clustered in smaller nucleons.

//--------------------------------------------------------------------------

// Initialize. Check if this nucleon is implemented

bool ClusterModel::init() {

  // Initialize hard core.
  initHardCore();
  // So far only use this model for 4He
  vector<int> implemented = {1000020040};
  bool isImplemented = false;
  for (int i = 0, N = implemented.size(); i < N; ++i) {
    if (id() == implemented[i]) {
      isImplemented = true;
      break;
    }
  }
  if (!isImplemented) {
    loggerPtr->ABORT_MSG("nucleus has no valid cluster model",
      "(for id=" + to_string(id()) + ")");
    return false;
  }

  // Set up the internal nucleus model for clusters.
  // Hulthen
  nModelPtr = unique_ptr<HulthenModel>();
  nModelPtr->initPtr(1000010020, isProj, *infoPtr);
  nModelPtr->init();

  return true;
}

//--------------------------------------------------------------------------

// Generate the nucleons.

vector<Nucleon> ClusterModel::generate() const {
  vector<Nucleon> ret;
  auto h1 = nModelPtr->generate();
  auto h2 = nModelPtr->generate();
  ret.insert(ret.end(), h1.begin(), h1.end());
  ret.insert(ret.end(), h2.begin(), h2.end());

  return ret;
}

//==========================================================================

} // end namespace Pythia8
