// VinciaMG5MEs.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Main interface methods for MadGraph 5 C++ matrix elements
// (MEs generated with the MG5 PY8Kernels Plugin)

#include "Pythia8/VinciaMG5MEs.h"

namespace Pythia8 {

//==========================================================================

// The VinciaMG5MEs matrix-element interface

//--------------------------------------------------------------------------

// Set pointers to required PYTHIA 8 objects.

void VinciaMG5MEs::initPtr(Info* infoPtrIn,
  SusyLesHouches* slhaPtrIn, VinciaCommon* vinComPtrIn) {
  infoPtr      = infoPtrIn;
  coupSMPtr       = infoPtr->coupSMPtr;
  particleDataPtr = infoPtr->particleDataPtr;
  settingsPtr     = infoPtr->settingsPtr;
  rndmPtr         = infoPtr->rndmPtr;
  slhaPtr         = slhaPtrIn;
  vinComPtr       = vinComPtrIn;
  isInitPtr       = true;
}

//--------------------------------------------------------------------------

// Main initialization routine.

bool VinciaMG5MEs::init() {

  // Check if pointers initialized.
  verbose = settingsPtr->mode("Vincia:verbose");
  if (verbose > normal) printOut("VinciaMG5MEs::init","...");
  if (!isInitPtr) {
    printOut("VinciaMG5MEs::init",
      "Cannot initialize, pointers not set.");
    return false;
  }
  isInit = true;

  // Set colour depth (TODO: via "Vincia:matchingFullColour").
  colourDepth = 1;

  // Initialise Parameters_sm (only if not using dummy).
#ifdef MG5MES

  // Create new model instance.
  if (verbose >= quiteloud) printOut("VinciaMG5MEs::init",
    "   setting MG C++ masses, widths, couplings...");
  Parameters_sm* modelPtr = new Parameters_sm();
  modelPtr->setIndependentParameters(particleDataPtr,coupSMPtr,slhaPtr);
  modelPtr->setIndependentCouplings();
  if (verbose >= quiteloud) {
    modelPtr->printIndependentParameters();
    modelPtr->printIndependentCouplings();
  }

  // In the VINCIA context, alphaS_MGME = 1/4pi (- > gS = 1; we
  // control the running separately). Thus, even the MG5 "dependent"
  // parameters only need to be set once.

  // Alternatively, we could evaluate the QCD coupling at MZ but should
  // then use a coupling definition from a Vincia parameter list rather
  // than PYTHIA's couplings.
  //    double muDefME  = 91.2;
  //    double alpS = coupSMPtr->alphaS(pow2(muDefME));

  // The following is equivalent to
  // PY8MEs::updateModelDependentCouplingsWithPY8(...)
  double alpS = 1.0 / ( 4 * M_PI );
  modelPtr->setDependentParameters(particleDataPtr, coupSMPtr, slhaPtr,
    alpS);
  modelPtr->setDependentCouplings();

  // Construct MG5 process library.
  if (verbose >= superdebug) printOut("VinciaMG5MEs::init()",
      "   attempting to construct mg5lib");
  if (mg5libPtr != nullptr && mg5libPtr != 0) delete mg5libPtr;
  mg5libPtr = new PY8MEs(modelPtr);
  return true;
#else
  mg5libPtr = nullptr;
  isInit = false;
  return false;
#endif

}

//--------------------------------------------------------------------------

// Get pointer to matrix element, e.g. to check if process exists in
// library.x

#ifdef MG5MES
PY8ME* VinciaMG5MEs::getProcess(vector<int> idIn, vector<int> idOut,
  set<int> sChan) {
    if (verbose >= superdebug) {
      cout << " VinciaMG5interface::getProcess(): checking for process";
      for (int i = 0; i < (int)idIn.size(); ++i) cout << " " << idIn[i];
      cout << " > ";
      for (int i = 0; i < (int)idOut.size(); ++i) cout << " " << idOut[i];
      cout << endl;
    }
    if (mg5libPtr != nullptr && mg5libPtr != 0)
      return mg5libPtr->getProcess(idIn, idOut, sChan);
    cout << "      returning NULL" << endl;
    return nullptr;
}
#endif

//--------------------------------------------------------------------------

// Get the matrix element squared for a particle state.

#ifdef MG5MES
double VinciaMG5MEs::ME2(vector<Particle> state, int nIn) {

  // Prepare vector of incoming ID codes.
  if (nIn <= 0) return -1;
  else if (state.size() - nIn < 1) return -1;
  vector<int> idIn, helOrg, col, acol;
  vector<Vec4> momenta;
  idIn.push_back(state[0].id());
  momenta.push_back(state[0].p());
  helOrg.push_back(state[0].pol());
  col.push_back(state[0].col());
  acol.push_back(state[0].acol());
  if (nIn == 2) {
    idIn.push_back(state[1].id());
    momenta.push_back(state[1].p());
    helOrg.push_back(state[1].pol());
    col.push_back(state[1].col());
    acol.push_back(state[1].acol());
  }
  // Prepare vector of outgoing ID codes.
  vector<int> idOut;
  for (int i=nIn; i<(int)state.size(); ++i) {
    idOut.push_back(state[i].id());
    momenta.push_back(state[i].p());
    helOrg.push_back(state[i].pol());
    col.push_back(state[i].col());
    acol.push_back(state[i].acol());
  }
  // Currently not using the option to request specific s-channels.
  set<int> sChannels;

  // Access the process.
  process_specifier proc_spec = mg5libPtr->getProcessSpecifier(idIn, idOut,
    sChannels);
  process_accessor proc_handle = mg5libPtr->getProcess(proc_spec);

  // Return right away if unavailable.
  if (proc_handle.second.second < 0) return -1;

  // Convert momenta and colours to MG5 format (energy first entry).
  vector< vector<double> > momentaMG5;
  vector< int > colacolMG5;
  for (int i = 0; i < (int)momenta.size(); ++i) {
    vector<double> pNow;
    pNow.push_back(momenta[i].e());
    pNow.push_back(momenta[i].px());
    pNow.push_back(momenta[i].py());
    pNow.push_back(momenta[i].pz());
    momentaMG5.push_back(pNow);
    colacolMG5.push_back(col[i]);
    colacolMG5.push_back(acol[i]);
  }

  vector<int> i9;
  // Check if we are doing a (partial) helicity sum.
  for (int i = 0; i < (int)helOrg.size(); ++i) {
    // Save indices of unpolarised partons.
    if (helOrg[i] == 9) i9.push_back(i);
  }

  // Explicitly sum over any hel = 9 helicities.
  vector< vector<int> > helConf;
  helConf.push_back(helOrg);
  while (i9.size() > 0) {
    int i  = i9.back();
    int id = (i < nIn) ? idIn[i] : idOut[i-nIn];
    // How many spin states.
    int nS = particleDataPtr->spinType(id);
    // Massless particles max have max 2 (physical) spin states.
    if (particleDataPtr->m0(id) == 0.0) nS=min(nS,2);
    // Create nS copies of helConf, one for each spin state.
    int helConfSizeNow = helConf.size();
    for (int iCopy = 1; iCopy <= nS; ++iCopy) {
      // Set hel for this particle in this copy.
      // Start from -1, then 1, then 0 (if 3 states).
      int h = -1;
      if (iCopy == 2) h = 1;
      else if (iCopy == 3) h = 0;
      else if (iCopy == 4) h = -2;
      else if (iCopy == 5) h = 2;
      for (int iHelConf=0; iHelConf<helConfSizeNow; ++iHelConf) {
        vector<int> helNow = helConf[iHelConf];
        helNow[i] = h;
        // First copy: use existing.
        if (iCopy == 1) helConf[iHelConf] = helNow;
        // Subsequent copies: create new.
        else helConf.push_back(helNow);
      }
    }
    // Remove the particle whose helicities have been summed over.
    i9.pop_back();
  }
  if (verbose >= superdebug) {
    cout << " in = ";
    for (int i = 0; i < (int)idIn.size(); ++i) cout << idIn[i] << " ";
    cout << "   out = ";
    for (int i = 0; i < (int)idOut.size(); ++i) cout << idOut[i] << " ";
    cout << endl;
    cout << " number of helicity configurations = " << helConf.size() << endl;
    for (int i = 0; i < (int)helConf.size(); ++i) {
      cout << "   helConf " << i;
      for (int j = 0; j < (int)helConf[i].size(); ++j)
        cout << " " << helConf[i][j];
      cout << endl;
    }
  }

  // Set properties and return ME2 value.
  PY8ME* proc_ptr = proc_handle.first;
  vec_int perms = proc_handle.second.first;
  int proc_ID = proc_handle.second.second;
  proc_ptr->setMomenta(momentaMG5);
  proc_ptr->setProcID(proc_ID);
  proc_ptr->setPermutation(perms);
  proc_ptr->setColors(colacolMG5);

  // Compute helicity sum (and save helicity components if needed later).
  double me2 = 0.0;
  me2hel.clear();
  for (int iHC=0; iHC<(int)helConf.size(); ++iHC) {

    // Note. could check here if the helConf is physical (consistent
    // with spins of particles and conserving angular momentum).
    proc_ptr->setHelicities(helConf[iHC]);
    double me2now = proc_ptr->sigmaKin();
    // MG may produce inf/nan for unphysical hel combinations.
    if ( !isnan(me2now) && !isinf(me2now) ) {
      // Save helicity matrix element for possible later use.
      me2hel[helConf[iHC]] = me2now;
      // Add this helicity ME to me2
      me2 += me2now;
    }
  }
  me2 *= double(proc_ptr->getSymmetryFactor());
  return me2;

}
#else
// Pure dummy implementation if MG5MES plugin not linked.
double VinciaMG5MEs::ME2(vector<Particle>, int) {return -1.;}
#endif

//--------------------------------------------------------------------------

// Set helicities for a particle state.

bool VinciaMG5MEs::selectHelicities(vector<Particle>& state, int nIn){

  // Get the matrix element (automatically sums over any h=9 particles).
  double me2sum = ME2(state, nIn);
  if (verbose >= superdebug)
    cout << " VinciaMG5MEs::selectHelicities(): "
         << scientific<<me2sum << endl;

  // Did we find the ME, (ME2() returns -1 if no ME found).
  if (me2sum <= 0.) return false;

  // Check how many helicity configurations we summed over.
  int nHelConf = me2hel.size();
  if (nHelConf <= 0) return false;

  // Random number between zero and me2sum (trivial if only one helConf).
  double ranHelConf = 0.0;
  vector<int> hSelected;
  if (nHelConf >= 2) ranHelConf = rndmPtr->flat() * me2sum;
  for (map< vector<int>, double>::iterator it=me2hel.begin();
       it != me2hel.end(); ++it) {
    // Progressively subtract each me2hel and check when we cross zero.
    ranHelConf -= it->second;
    if (ranHelConf <= 0.0) {
      hSelected = it->first;
      break;
    }
  }
  if (ranHelConf > 0.) return false;

  // Set helicity configuration.
  for (int i = 0; i < (int)state.size(); ++i) state[i].pol(hSelected[i]);
  if (verbose >= superdebug)
    cout << " VinciaMG5MEs::selectHelicities(): selected " <<
      makeLabel(hSelected, nIn, false) << endl;
  return true;

}

//--------------------------------------------------------------------------

// Convert process id codes (or helicity values) to string.

string VinciaMG5MEs::makeLabel(vector<int>& id, int nIn,
  bool useNames) {
  string label = "{";
  for (int i = 0; i < (int)id.size(); ++i) {
    string idName;
    if (useNames && id[i] != 0) idName = particleDataPtr->name(id[i]);
    else idName = num2str(id[i]);
    if (i == nIn-1) idName += " ->";
    label += idName+" ";
  }
  label += "}";
  return label;
}

//==========================================================================

} // end namespace Pythia8
