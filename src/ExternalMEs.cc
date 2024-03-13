// ExternalMEs.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Peter Skands, Stefan Prestel, Philip Ilten, Torbjorn
// Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Virtual interface for external matrix element plugins for parton showers.

#include "Pythia8/ExternalMEs.h"
#include "Pythia8/SusyCouplings.h"

namespace Pythia8 {

//==========================================================================

// The parton shower matrix-element interface.

//--------------------------------------------------------------------------

// Set pointers to required PYTHIA 8 objects.

void ExternalMEs::initPtrs(Info* infoPtrIn) {
  infoPtr         = infoPtrIn;
  loggerPtr       = infoPtr->loggerPtr;
  coupSMPtr       = infoPtr->coupSMPtr;
  particleDataPtr = infoPtr->particleDataPtr;
  settingsPtr     = infoPtr->settingsPtr;
  rndmPtr         = infoPtr->rndmPtr;
  slhaPtr         = infoPtr->coupSUSYPtr->slhaPtr;
  isInitPtr       = true;
}

//--------------------------------------------------------------------------

// Fill a vector of IDs.

void ExternalMEs::fillIds(const Event& event, vector<int>& in,
  vector<int>& out, int iBeg) const {
  // Start at i = iBeg, and go to end of event.
  for (int i = iBeg; i < event.size(); ++i) {
    if ( !event[i].isFinal() ) in.push_back(event[i].id());
    else out.push_back( event[i].id() );
  }
}

//--------------------------------------------------------------------------

// Fill a vector of momenta.

void ExternalMEs::fillMoms(const Event& event, vector<Vec4>& p, int iBeg)
  const {
  // Start at i = iBeg, and go to end of event.
  for (int i = iBeg; i < event.size(); ++i) p.push_back(event[i].p());
}

//--------------------------------------------------------------------------

// Fill a vector of colors.

void ExternalMEs::fillCols(const Event& event, vector<int>& colors, int iBeg)
  const {
  // Start at i = iBeg, and go to end of event.
  for (int i = iBeg; i < event.size(); ++i) {
    colors.push_back(event[i].col());
    colors.push_back(event[i].acol());
  }
}


//--------------------------------------------------------------------------

// Return the momenta.

vector<vector<double> > ExternalMEs::fillMoms(const Event& event, int iBeg)
  const {
  vector<Vec4> p;
  fillMoms(event, p, iBeg);
  vector< vector<double> > ret;
  for (int i = 0; i < int(p.size()); i++ ) {
    vector<double> p_tmp(4, 0.);
    p_tmp[0] = isnan(p[i].e())  ? 0.0 : p[i].e() ;
    p_tmp[1] = isnan(p[i].px()) ? 0.0 : p[i].px();
    p_tmp[2] = isnan(p[i].py()) ? 0.0 : p[i].py();
    p_tmp[3] = isnan(p[i].pz()) ? 0.0 : p[i].pz();
    ret.push_back(p_tmp);
  }
  return ret;
}

//==========================================================================

// Helicity sampler which uses external matrix elements.

//--------------------------------------------------------------------------

// Set helicities for a particle state.

bool HelicitySampler::selectHelicities(vector<Particle>& state, bool force) {
  // Check we have access to a matrix element generator.
  if (!isInitPtr) return false;

  // If force = true, erase any existing helicities.
  if (force)
    for (int i=0; i<(int)state.size(); ++i) state[i].pol(9);

  // Save current settings of ME generator.
  int helModeSave = mePluginPtr->helicityMode();
  int colModeSave = mePluginPtr->colourMode();
  bool inclSymFacSave         = mePluginPtr->includeSymmetryFac();
  bool inclHelicityAvgFacSave = mePluginPtr->includeHelicityAvgFac();
  bool inclColourAvgFacSave   = mePluginPtr->includeColourAvgFac();

  // Set that we want an explicit helicity sum.
  mePluginPtr->setHelicityMode(0);
  // LC amplitudes only.
  mePluginPtr->setColourMode(1);
  // Include all averaging and symmetry factors.
  mePluginPtr->setIncludeSymmetryFac(true);
  mePluginPtr->setIncludeHelicityAvgFac(true);
  mePluginPtr->setIncludeColourAvgFac(true);

  // Calculate the matrix element and fetch helicity amplitudes.
  double me2sum = mePluginPtr->calcME2(state);
  if (me2sum <= 0.) return false;
  map<vector<int>, double> me2hel = mePluginPtr->getHelicityAmplitudes();

  // Restore ME generator settings.
  mePluginPtr->setHelicityMode(helModeSave);
  mePluginPtr->setColourMode(colModeSave);
  mePluginPtr->setIncludeSymmetryFac(inclSymFacSave);
  mePluginPtr->setIncludeHelicityAvgFac(inclHelicityAvgFacSave);
  mePluginPtr->setIncludeColourAvgFac(inclColourAvgFacSave);

  // Check how many helicity configurations we summed over.
  int nHelConf = me2hel.size();
  if (nHelConf <= 0) return false;

  // Add up all helicity matrix elements.
  double me2helsum = 0.;
  for (auto it = me2hel.begin(); it!=me2hel.end(); ++it)
    me2helsum += it->second;

  // Random number between zero and me2sum (trivial if only one helConf).
  double ranHelConf = 0.0;
  vector<int> hSelected;
  if (nHelConf >= 2) ranHelConf = rndmPtr->flat() * me2helsum;
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
  return true;

}

} // end namespace Pythia8
