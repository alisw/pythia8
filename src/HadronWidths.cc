// HadronWidths.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the HadronWidths class.

#include "Pythia8/HadronWidths.h"

namespace Pythia8 {

//--------------------------------------------------------------------------

// Static methods for reading xml files.

static string attributeValue(string line, string attribute) {
  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);
}

static int intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;
}

static double doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;
}

static void completeTag(istream& stream, string& line) {
  while (line.find(">") == string::npos) {
    string addLine;
    if (!getline(stream, addLine)) break;
    line += " " + addLine;
  }
}

//==========================================================================

// The HadronWidths class.

//--------------------------------------------------------------------------

// Initialize.

bool HadronWidths::init(string path) {

  ifstream stream(path);
  if (!stream.is_open()) {
    loggerPtr->ERROR_MSG("unable to open file");
    return false;
  }

  return init(stream);
}

//--------------------------------------------------------------------------

// Initialize.

bool HadronWidths::init(istream& stream) {

  string line;

  while (getline(stream, line)) {

    string word1;
    if (!(istringstream(line) >> word1))
      continue;

    if (word1 == "<width") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");
      auto entryIter = entries.find(id);
      if (entryIter != entries.end() && entryIter->second.isUserDefined) {
        loggerPtr->ERROR_MSG("resonance is defined more than once",
          std::to_string(id));
        continue;
      }

      double left  = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      // Insert resonance in entries
      LinearInterpolator widths(left, right, data);
      entries.emplace(id, HadronWidthEntry{ widths, {}, false });

      // Insert resonance in signature map
      int signature = getSignature(particleDataPtr->isBaryon(id),
                                   particleDataPtr->chargeType(id));

      auto iter = signatureToParticles.find(signature);
      if (iter == signatureToParticles.end())
        // If signature has not been used yet, insert a new vector into the map
        signatureToParticles.emplace(signature, vector<int> { id });
      else
        // If signature has been used already, add id to the existing vector
        iter->second.push_back(id);
    }
    else if (word1 == "<partialWidth") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");

      auto entryIter = entries.find(id);
      if (entryIter == entries.end()) {
        loggerPtr->ERROR_MSG(
          "got partial width for a particle with undefined total width",
          std::to_string(id));
        continue;
      }

      int lType = intAttributeValue(line, "lType");

      istringstream productStr(attributeValue(line, "products"));
      int prod1, prod2;
      productStr >> prod1;
      productStr >> prod2;

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      HadronWidthEntry& entry = entryIter->second;
      LinearInterpolator widths(entry.width.left(), entry.width.right(), data);

      // Generate key to ensure canonical ordering of decay products.
      pair<int, int> key = getKey(id, prod1, prod2);
      double mThreshold = particleDataPtr->mMin(key.first)
                        + particleDataPtr->mMin(key.second);

      // Insert new decay channel.
      entry.decayChannels.emplace(key, ResonanceDecayChannel {
        widths, key.first, key.second, lType, mThreshold });
    }
  }

  // Search for newly added particles.
  for (auto& pdtEntry : *particleDataPtr) {
    int id = pdtEntry.first;
    ParticleDataEntryPtr pdt = pdtEntry.second;
    auto iter = entries.find(id);
    // Parameterize particle if it is not yet listed and has variable width.
    if (iter == entries.end() && pdt->varWidth()) {
      // Check that has a resonance decay channel, i.e. into two hadrons.
      for (int iChannel = 0; iChannel < pdt->sizeChannels(); ++iChannel) {
        DecayChannel& channel = pdt->channel(iChannel);
        if (channel.multiplicity() == 2
         && particleDataPtr->isHadron(channel.product(0))
         && particleDataPtr->isHadron(channel.product(1)) ) {
          loggerPtr->INFO_MSG("parameterizing new resonance",
            to_string(id), true);
          parameterize(id);
          break;
        }
      }
    }
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Check whether input data is valid and matches particle data.

bool HadronWidths::check() {

  // Check that all resonance entries make sense.
  for (auto& entryPair : entries) {
    int id = entryPair.first;
    HadronWidthEntry& entry = entryPair.second;

    // Check that entry id actually corresponds to a particle.
    if (!particleDataPtr->isParticle(id)) {
      loggerPtr->ERROR_MSG("resonance is not a particle", to_string(id));
      return false;
    }

    // Check that entry id is positive (antiparticles are handled by symmetry).
    if (id < 0) {
      loggerPtr->ERROR_MSG("resonance is an anti-particle", to_string(id));
      return false;
    }

    // Check that entry id is hadron.
    if (!particleDataPtr->isHadron(id)) {
      loggerPtr->ERROR_MSG("resonance is not a hadron", to_string(id));
      return false;
    }

    // Check that mass boundaries are same in particle entry and widths entry.
    if (particleDataPtr->mMin(id) < entry.width.left()) {
      loggerPtr->WARNING_MSG("inconsistent lower mass bound",
        to_string(id));
    }
    if (particleDataPtr->mMax(id) > entry.width.right()) {
      loggerPtr->WARNING_MSG("inconsistent upper mass bound",
        to_string(id));
    }

    // Check that all decay channels make sense.
    for (auto channelPair : entry.decayChannels) {
      ResonanceDecayChannel& channel = channelPair.second;
      int idA = channel.prodA, idB = channel.prodB;
      string channelStr = to_string(id) + " --> "
          + to_string(idA) + " + " + to_string(idB);

      // Check that decay product ids actually correspond to particles.
      for (int idProd : { idA, idB })
      if (!particleDataPtr->isParticle(idProd)) {
        loggerPtr->ERROR_MSG("decay product is not a particle",
          to_string(idProd));
        return false;
      }
      // Check that lType makes sense.
      if (channel.lType <= 0) {
        loggerPtr->ERROR_MSG("decay channel does not specify a valid lType",
          channelStr);
        return false;
      }

      // Check that decay conserves charge.
      if (particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB)
        != particleDataPtr->chargeType(id)) {
        loggerPtr->ERROR_MSG("decay does not conserve charge", channelStr);
        return false;
      }
    }
  }

  for (auto& entry : *particleDataPtr) {
    if (entry.second->varWidth() && !hasData(entry.first)) {
      loggerPtr->WARNING_MSG(
        "particle uses mass dependent width, but width is not defined",
        to_string(entry.first));
    }
  }

  return true;
}

//--------------------------------------------------------------------------

// Gets key for the decay and flips idR if necessary

pair<int, int> HadronWidths::getKey(int& idR, int idA, int idB) const {

  if (idR < 0) {
    idR = -idR;
    idA = particleDataPtr->antiId(idA);
    idB = particleDataPtr->antiId(idB);
  }

  if (abs(idA) < abs(idB) || (idA == -idB && idA < 0)) return {idB, idA};
  else return {idA, idB};
}

//--------------------------------------------------------------------------

// Get signature of system based on total baryon number and electric charge.

int HadronWidths::getSignature(int baryonNumber, int charge) const {
  return 100 * baryonNumber + 10 * abs(charge);
}

//--------------------------------------------------------------------------

// Get whether the specified incoming particles can form a resonance.

bool HadronWidths::hasResonances(int idA, int idB) const {

  ParticleDataEntryPtr entryA = particleDataPtr->findParticle(idA);
  ParticleDataEntryPtr entryB = particleDataPtr->findParticle(idB);
  if (!entryA || !entryB) {
    loggerPtr->ERROR_MSG("invalid input particle ids");
    return false;
  }

  // Get signature for system and look only for resonances that matches it.
  int baryonNumber = entryA->isBaryon() + entryB->isBaryon();
  int charge = entryA->chargeType(idA) + entryB->chargeType(idB);
  int signature = getSignature(baryonNumber, charge);
  auto iter = signatureToParticles.find(signature);
  if (iter == signatureToParticles.end())
    return false;

  // For resonances that matches signature, check that decay channel exists.
  for (int idR : iter->second) {
    if ( canDecay(idR, idA, idB)
      || (particleDataPtr->hasAnti(idR) && canDecay(-idR, idA, idB)) )
      return true;
  }

  // No resonances found.
  return false;
}

//--------------------------------------------------------------------------

// Get all implemented resonances.

set<int> HadronWidths::getResonances() const {
  set<int> resonances;
  for (auto& p : entries) resonances.insert(p.first);
  return resonances;
}

//--------------------------------------------------------------------------

// Get resonances that can be formed by the specified incoming particles.

set<int> HadronWidths::getResonances(int idA, int idB) const {

  ParticleDataEntryPtr entryA = particleDataPtr->findParticle(idA);
  ParticleDataEntryPtr entryB = particleDataPtr->findParticle(idB);
  if (!entryA || !entryB) {
    loggerPtr->ERROR_MSG("invalid input particle ids");
    return set<int>();
  }

  // Get signature for system and look only for resonances that matches it.
  int baryonNumber = entryA->isBaryon() + entryB->isBaryon();
  int charge = entryA->chargeType(idA) + entryB->chargeType(idB);
  int signature = getSignature(baryonNumber, charge);
  auto iter = signatureToParticles.find(signature);
  if (iter == signatureToParticles.end())
    return set<int>();

  // For resonances that matches signature, check that decay channel exists.
  set<int> resonances;
  for (int idR : iter->second) {
    if (canDecay(idR, idA, idB))
      resonances.insert(idR);
    if (particleDataPtr->hasAnti(idR) && canDecay(-idR, idA, idB))
      resonances.insert(-idR);
  }

  // For pi0pi0 and pi+pi-, add f0(500) explicitly.
  if ( (idA == 111 && idB == 111)
    || (abs(idA) == 211 && abs(idB) == 211 && idA * idB < 0) )
    resonances.insert(9000221);

  // Done.
  return resonances;
}

//--------------------------------------------------------------------------

// Get whether the resonance can decay into the specified products.

bool HadronWidths::canDecay(int idR, int idA, int idB) const {

  // Get key and flip idR if necessary.
  pair<int, int> key = getKey(idR, idA, idB);

  auto entryIter = entries.find(idR);
  if (entryIter == entries.end())
    return false;

  auto channelIter = entryIter->second.decayChannels.find(key);
  return channelIter != entryIter->second.decayChannels.end();
}

//--------------------------------------------------------------------------

// Get the total width of the specified particle at the specified mass.

double HadronWidths::width(int id, double m) const {

  // Find particle data entry.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(id);
  if (!entry) {
    loggerPtr->ERROR_MSG("particle does not exist", to_string(id));
    return 0.;
  }

  // Check that mass is within range.
  if (m < entry->mMin() || m > entry->mMax())
    return 0.;

  // If particle is not resonance, use Breit-Wigner width.
  if (!entry->varWidth())
    return entry->mWidth();

  // For resonances, get width from parameterization.
  auto iter = entries.find(abs(id));
  if (iter == entries.end()) {
    loggerPtr->WARNING_MSG("particle is resonance, but is not parameterized",
      to_string(id));
    return entry->mWidth();
  }
  return iter->second.width(m);
}

//--------------------------------------------------------------------------

// Get the partial width for the specified decay channel of the particle.

double HadronWidths::partialWidth(int idR, int idA, int idB, double m) const {

  // Get key and flip idR if necessary.
  pair<int, int> key = getKey(idR, idA, idB);

  // Find particle data entry.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(idR);
  if (!entry) {
    loggerPtr->ERROR_MSG("particle does not exist", to_string(idR));
    return 0.;
  }

  // Check that mass is within range.
  if (m < entry->mMin() || m > entry->mMax())
    return 0.;

  // If particle is not resonance, use Breit-Wigner width.
  if (!particleDataPtr->varWidth(idR))
    return particleDataPtr->mWidth(idR) * br(idR, idA, idB, m);

  // For resonances, get width from parameterization.
  auto entryIter = entries.find(idR);
  if (entryIter == entries.end()) {
    loggerPtr->WARNING_MSG("particle is resonance, but is not parameterized",
      to_string(idR));
    return 0.;
  }

  auto channelIter = entryIter->second.decayChannels.find(key);
  if (channelIter == entryIter->second.decayChannels.end())
    return 0.;

  return (m <= channelIter->second.mThreshold) ? 0.
        : channelIter->second.partialWidth(m);
}

//--------------------------------------------------------------------------

// Get the branching ratio for the specified decay channel of the particle.

double HadronWidths::br(int idR, int idA, int idB, double m) const {

  // Get key and flip idR if necessary.
  pair<int, int> key = getKey(idR, idA, idB);

  // Find particle data entry.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(idR);
  if (!entry) {
    loggerPtr->ERROR_MSG("particle does not exist", to_string(idR));
    return 0.;
  }

  // Check that mass is within range.
  if (m < entry->mMin() || m > entry->mMax())
    return 0.;

  // If particle is not resonance, get BR from particle database.
  if (!entry->varWidth()) {
    for (int iChannel = 0; iChannel < entry->sizeChannels(); ++iChannel) {
      DecayChannel& channel = entry->channel(iChannel);
      if (channel.multiplicity() != 2)
        continue;
      if ( (channel.product(0) == idA && channel.product(1) == idB)
        || (channel.product(0) == idB && channel.product(1) == idA) )
        return channel.bRatio() * entry->mWidth();
    }

    // If resonance cannot decay into A B, return 0.
    return 0.;
  }

  // For resonances, get branching ratio from parameterization.
  auto entryIter = entries.find(idR);
  if (entryIter == entries.end()) {
    loggerPtr->WARNING_MSG("particle is resonance, but is not parameterized",
      to_string(idR));
    return 0.;
  }

  auto channelIter = entryIter->second.decayChannels.find(key);
  if (channelIter == entryIter->second.decayChannels.end())
    return 0.;

  double widthTotal = entryIter->second.width(m);
  if (widthTotal == 0.)
    return 0.;
  else
    return (m <= channelIter->second.mThreshold) ? 0.
          : channelIter->second.partialWidth(m) / widthTotal;
}

//--------------------------------------------------------------------------

// Get the mass distribution density for the particle at the specified mass,
// using the Breit-Wigner formula with a mass-dependent width.

double HadronWidths::mDistr(int id, double m) const  {
  double w = width(id, m);
  if (w == 0) return 0.;
  double m0 = particleDataPtr->m0(id);
  return 0.5 / M_PI * w / (pow2(m - m0) + 0.25 * w * w);
}

//--------------------------------------------------------------------------

// Pick a decay channel for the specified particle, together with phase
// space configuration. If successful, the results are written to the
// output arguments.

bool HadronWidths::pickDecay(int idDec, double m, int& idAOut, int& idBOut,
    double& mAOut, double& mBOut) {

  // Find particle data entry for decaying particle.
  ParticleDataEntryPtr pdEntry = particleDataPtr->findParticle(idDec);
  if (!pdEntry) {
    loggerPtr->ERROR_MSG("particle not found", to_string(idDec));
    return false;
  }

  // If antiparticle, flip the resonance and then later flip decay products.
  bool isAnti = (idDec < 0);
  if (isAnti) idDec = -idDec;

  // Find table entry for decaying particle.
  auto entriesIter = entries.find(idDec);
  if (entriesIter == entries.end()) {
    loggerPtr->ERROR_MSG("particle is not parameterized", to_string(idDec));
    return false;
  }
  HadronWidthEntry& entry = entriesIter->second;

  // Get list of channels that are currently on.
  vector<ResonanceDecayChannel*> channelsList;
  vector<double> branchingRates;
  bool gotAny = false;
  for (auto& channel : entry.decayChannels) {

    // Check that channel is open at this mass.
    if (m <= channel.second.mThreshold)
      continue;

    // Check that channel is on in the particle database.
    int prodA = channel.first.first, prodB = channel.first.second;
    bool isOn = false;
    for (int iChannel = 0; iChannel < pdEntry->sizeChannels(); ++iChannel) {
      DecayChannel& pdChannel = pdEntry->channel(iChannel);
      if (pdChannel.multiplicity() != 2) continue;
      int pdProdA = pdChannel.product(0), pdProdB = pdChannel.product(1);
      if ( (pdProdA == prodA && pdProdB == prodB)
        || (pdProdB == prodA && pdProdA == prodB) ) {
        isOn = (pdChannel.onMode() == 1)
            || (pdChannel.onMode() == 2 && idDec == pdEntry->id())
            || (pdChannel.onMode() == 3 && idDec == pdEntry->antiId());
        break;
      }
    }
    if (!isOn)
      continue;

    // If channel is on, add to list if width is positive.
    double widthNow = channel.second.partialWidth(m);
    if (widthNow > 0.) {
      gotAny = true;
      channelsList.push_back(&channel.second);
      branchingRates.push_back(widthNow);
    }
  }

  // If no channels are on, the decay fails.
  if (!gotAny) {
     loggerPtr->ERROR_MSG("no channels have positive widths",
       "for " + to_string(idDec) + " @ " + to_string(m) + " GeV");
    return false;
  }

  auto* channelNow = channelsList[rndmPtr->pick(branchingRates)];
  int prodA = channelNow->prodA;
  int prodB = channelNow->prodB;

  // Pick masses of decay products.
  double mA, mB;
  if (!pickMasses(prodA, prodB, m, mA, mB, channelNow->lType)) {
    loggerPtr->ERROR_MSG("failed to pick masses",
      "for " + to_string(idDec) + " --> " + to_string(prodA)
      + " + " + to_string(prodB) + " @ " + to_string(m));
    return false;
  }

  // Write output values and done.
  idAOut = isAnti ? particleDataPtr->antiId(prodA) : prodA;
  idBOut = isAnti ? particleDataPtr->antiId(prodB) : prodB;
  mAOut = mA;
  mBOut = mB;
  return true;

}

//--------------------------------------------------------------------------

// Constants used when sampling masses.
static constexpr int    MAXLOOP        = 100;
static constexpr double MINWIDTH       = 0.001;
static constexpr double MAXWIDTHGROWTH = 2.;

//--------------------------------------------------------------------------

// Pick a pair of masses given pair invariant mass and angular momentum.

bool HadronWidths::pickMasses(int idA, int idB, double eCM,
  double& mAOut, double& mBOut, int lType) {

  // Minimal masses must be a possible choice.
  double mAMin = particleDataPtr->mMin(idA);
  double mBMin = particleDataPtr->mMin(idB);
  if (mAMin + mBMin >=  eCM) {
    loggerPtr->ERROR_MSG("energy is smaller than minimum masses");
    return false;
  }

  if (lType <= 0) {
    loggerPtr->ERROR_MSG("invalid angular momentum",
      "2l+1 = " + to_string(lType));
    return false;
  }

  // Done if none of the daughters have a width.
  double mAFix      = particleDataPtr->m0(idA);
  double gammaAFix  = particleDataPtr->mWidth(idA);
  bool hasFixWidthA = (gammaAFix > MINWIDTH);
  double mBFix      = particleDataPtr->m0(idB);
  double gammaBFix  = particleDataPtr->mWidth(idB);
  bool hasFixWidthB = (gammaBFix > MINWIDTH);
  mAOut             = mAFix;
  mBOut             = mBFix;
  if (!hasFixWidthA && !hasFixWidthB) return true;

  // Get width entries for particles with mass-depedent widths.
  bool hasVarWidthA = hasData(idA) && particleDataPtr->varWidth(idA);
  bool hasWidthA    = hasFixWidthA || hasVarWidthA;
  HadronWidthEntry* entryA = nullptr;
  if (hasVarWidthA) {
    auto iterA = entries.find( abs(idA) );
    if (iterA == entries.end()) {
      loggerPtr->ERROR_MSG("mass distribution for particle is not defined",
        to_string(idA));
      return false;
    }
    entryA = &iterA->second;
  }
  bool hasVarWidthB = hasData(idB) && particleDataPtr->varWidth(idB);
  bool hasWidthB    = hasFixWidthB || hasVarWidthB;
  HadronWidthEntry* entryB = nullptr;
  if (hasVarWidthB) {
    auto iterB = entries.find( abs(idB) );
    if (iterB == entries.end()) {
      loggerPtr->ERROR_MSG("mass distribution for particle is not defined",
        to_string(idB));
      return false;
    }
    entryB = &iterB->second;
  }

  // Parameters for mass selection.
  double mAMax  = min(particleDataPtr->mMax(idA), eCM - mBMin);
  if (hasVarWidthA) gammaAFix = entryA->width(mAFix);
  double bwAMin = (hasWidthA) ? atan(2. * (mAMin - mAFix) / gammaAFix) : 0.;
  double bwAMax = (hasWidthA) ? atan(2. * (mAMax - mAFix) / gammaAFix) : 0.;
  double mBMax  = min(particleDataPtr->mMax(idB), eCM - mAMin);
  if (hasVarWidthB) gammaBFix = entryB->width(mBFix);
  double bwBMin = (hasWidthB) ? atan(2. * (mBMin - mBFix) / gammaBFix) : 0.;
  double bwBMax = (hasWidthB) ? atan(2. * (mBMax - mBFix) / gammaBFix) : 0.;
  double p2Max  = (eCM*eCM - pow2(mAMin + mBMin))
                * (eCM*eCM - pow2(mAMin - mBMin));

  // Loop over attempts to pick the two masses simultaneously.
  double wtTot, gammaAVar, gammaBVar, bwAFix, bwBFix, bwAVar, bwBVar;
  for (int i = 0; i < MAXLOOP; ++i) {
    wtTot = 1.;

    // Simplify handling if full procedure does not seem to work.
    if (2 * i > MAXLOOP) {
      hasVarWidthA = false;
      hasVarWidthB = false;
    }
    if (4 * i > 3 * MAXLOOP) lType = 0;

    // Initially pick according to simple Breit-Wigner.
    if (hasWidthA) mAOut = mAFix + 0.5 * gammaAFix * tan(bwAMin
      + rndmPtr->flat() * (bwAMax - bwAMin));
    if (hasWidthB) mBOut = mBFix + 0.5 * gammaBFix * tan(bwBMin
      + rndmPtr->flat() * (bwBMax - bwBMin));

    // Correction given by BW(Gamma_now)/BW(Gamma_fix) for variable width.
    // Note: width not allowed to explode at large masses.
    if (hasVarWidthA) {
      gammaAVar = min(entryA->width(mAOut), MAXWIDTHGROWTH * gammaAFix);
      bwAVar    = gammaAVar / (pow2( mAOut - mAFix) + 0.25 * pow2(gammaAVar));
      bwAFix    = gammaAFix / (pow2( mAOut - mAFix) + 0.25 * pow2(gammaAFix));
      wtTot    *= bwAVar / (bwAFix * MAXWIDTHGROWTH);
    }
    if (hasVarWidthB) {
      gammaBVar = min(entryB->width(mBOut), MAXWIDTHGROWTH * gammaBFix);
      bwBVar    = gammaBVar / (pow2( mBOut - mBFix) + 0.25 * pow2(gammaBVar));
      bwBFix    = gammaBFix / (pow2( mBOut - mBFix) + 0.25 * pow2(gammaBFix));
      wtTot    *= bwBVar / (bwBFix * MAXWIDTHGROWTH);
    }

    // Weight by (p/p_max)^lType.
    if (mAOut + mBOut >= eCM) continue;
    double p2Ratio = (eCM*eCM - pow2(mAOut + mBOut))
                   * (eCM*eCM - pow2(mAOut - mBOut)) / p2Max;
    if (lType > 0) wtTot *= pow(p2Ratio, 0.5 * lType);
    if (wtTot > rndmPtr->flat()) {
      // Give warning message only for more severe cases
      if (4 * i > 3 * MAXLOOP)
      loggerPtr->WARNING_MSG("angular momentum and running widths not used");
      return true;
    }
  }

  // Last resort: pick masses within limits known to work, without weight.
  loggerPtr->WARNING_MSG("using last-resort simplified description");
  double mSpanNorm = (eCM - mAMin - mBMin) / (gammaAFix + gammaBFix);
  mAOut = mAMin + rndmPtr->flat() * mSpanNorm * gammaAFix;
  mBOut = mBMin + rndmPtr->flat() * mSpanNorm * gammaBFix;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Calculate the total width of the particle without using interpolation.

double HadronWidths::widthCalc(int id, double m) const {

  // Find particle data entry.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(id);
  if (entry == nullptr) {
    loggerPtr->ERROR_MSG("particle not found", to_string(id));
    return 0.;
  }
  if (m < entry->mMin() || m > entry->mMax())
    return 0.;

  // Sum contributions from all channels.
  double w = 0.;
  for (int iChannel = 0; iChannel < entry->sizeChannels(); ++iChannel)
    w += widthCalc(id, entry->channel(iChannel), m);
  return w;
}

//--------------------------------------------------------------------------

// Calculate partial width of the particle without using interpolation.

double HadronWidths::widthCalc(int id, int prodA, int prodB, double m) const {

  // Find particle data entry.
  pair<int, int> key = getKey(id, prodA, prodB);
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(id);
  if (entry == nullptr)
    return 0.;

  // Search for the matching decay channel.
  for (int iChannel = 0; iChannel < entry->sizeChannels(); ++iChannel) {
    DecayChannel& channel = entry->channel(iChannel);
    if (channel.multiplicity() > 2)
      continue;
    if ( (channel.product(0) == key.first && channel.product(1) == key.second)
      || (channel.product(1) == key.first && channel.product(0) == key.second))
      return widthCalc(id, channel, m);
  }

  // Decay channel not found.
  loggerPtr->ERROR_MSG("decay channel not found",
    to_string(id) + " --> " + to_string(prodA) + " " + to_string(prodB));
  return 0.;
}

//--------------------------------------------------------------------------

// Calculate partial width of the particle without using interpolation.

double HadronWidths::widthCalc(int id, Pythia8::DecayChannel& channel,
  double m) const {

  // Find particle data entry.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(id);
  if (entry == nullptr) {
    loggerPtr->ERROR_MSG("particle not found", to_string(id));
    return 0.;
  }
  if (m < entry->mMin() || m > entry->mMax())
    return 0.;

  // Multibody channels cannot have mass-dependent widths.
  if (channel.multiplicity() != 2)
    return entry->mWidth() * channel.bRatio();

  // Find particle data entries for decay products.
  auto prodA = particleDataPtr->findParticle(channel.product(0));
  auto prodB = particleDataPtr->findParticle(channel.product(1));
  if (m < prodA->mMin() + prodB->mMin())
    return 0.;

  // Get angular momentum for the outgoing two-body system.
  int lType;
  if (channel.meMode() >= 3 && channel.meMode() <= 7)
    lType = 2 * (channel.meMode() - 3) + 1;
  else if (channel.meMode() == 2)
    lType = 3;
  else
    lType = 1;

  // Calculate phase space at the specified mass.
  double pM = psSize(m, prodA, prodB, lType);
  if (pM == 0.)
    return 0.;
  double pMS = psSize(m, prodA, prodB, lType - 1);
  if (pMS == 0.)
    return 0.;

  // Calculate phase space at on-shell mass.
  double m0 = entry->m0();
  double pM0  = psSize(m0, prodA, prodB, lType);
  double pM0S = psSize(m0, prodA, prodB, lType - 1);
  if (pM0 <= 0 || pM0S <= 0) {
    loggerPtr->ERROR_MSG("on-shell decay is not possible",
      to_string(id) + " --> " + to_string(prodA->id())
       + " " + to_string(prodB->id()));
    return numeric_limits<double>::quiet_NaN();
  }

  // Return mass-dependent partial width, using UrQMD approach.
  double gamma0 = channel.bRatio() * entry->mWidth();
  return gamma0 * (m0 / m) * (pM / pM0) * 1.2 / (1. + 0.2 * pMS / pM0S);
}

//--------------------------------------------------------------------------

// Regenerate parameterization for particle and its decay products if needed.

bool HadronWidths::parameterize(int id, int precision) {

  // Get particle entry and validate input.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(id);

  if (entry == nullptr) {
    loggerPtr->ERROR_MSG("particle does not exist", to_string(id));
    return false;
  }
  if (precision <= 1) {
    loggerPtr->ERROR_MSG("precision must be at least 2");
    return false;
  }
  if (entry->mMin() >= entry->mMax()) {
    loggerPtr->ERROR_MSG("particle has fixed mass", to_string(id));
    return false;
  }
  if (!entry->varWidth())
    loggerPtr->WARNING_MSG("particle does not have mass-dependent width",
      to_string(id));

  // Do recursive parameterization (no validation checks in recursive calls)
  return parameterizeRecursive(id, precision);
}

//--------------------------------------------------------------------------

// Generate parameterization for particle and its decay products if needed.

bool HadronWidths::parameterizeRecursive(int id, int precision) {

  // End recursion if data has already been generated.
  if (hasData(id))
    return true;

  // Get particle entry.
  ParticleDataEntryPtr entry = particleDataPtr->findParticle(id);

  // Check whether any decay products must be parameterized first.
  for (int iChannel = 0; iChannel < entry->sizeChannels(); ++iChannel) {
    DecayChannel& channel = entry->channel(iChannel);
    if (channel.multiplicity() == 2) {
      auto prodA = particleDataPtr->findParticle(channel.product(0));
      auto prodB = particleDataPtr->findParticle(channel.product(1));

      // Recursive call to parameterize decay product widths if necessary.
      if (prodA->varWidth() && !hasData(prodA->id()))
        if (!parameterizeRecursive(prodA->id(), precision)) return false;
      if (prodB->varWidth() && !hasData(prodB->id()))
        if (!parameterizeRecursive(prodB->id(), precision)) return false;
    }
  }

  // Perform the actual parameterization of this particle.
  map<pair<int, int>, ResonanceDecayChannel> partialWidths;
  vector<double> totalWidthData(precision);

  double mMin = entry->mMin(), mMax = entry->mMax();
  double dm = (mMax - mMin) / (precision - 1);
  entry->rescaleBR(1.0);

  // Parameterize all channels.
  for (int iChannel = 0; iChannel < entry->sizeChannels(); ++iChannel) {
    DecayChannel& channel = entry->channel(iChannel);

    // For multibody channels, mass-dependent widths are not defined,
    // but they should still count towards the total width.
    if (channel.multiplicity() != 2) {
      for (double& dataPoint : totalWidthData)
        dataPoint += entry->mWidth() * channel.bRatio();
      continue;
    }

    // Create key to put decay products in canonical order.
    pair<int, int> key = getKey(id, channel.product(0), channel.product(1));
    int prodA = key.first, prodB = key.second;

    // Calculate widths at regular mass intervals.
    vector<double> widthData(precision);
    for (int j = 0; j < precision; ++j) {
      double m = mMin + j * dm;
      widthData[j] = widthCalc(entry->id(), channel, m);
      totalWidthData[j] += widthData[j];
    }

    // Get two-body angular momentum.
    int lType;
    if (channel.meMode() >= 3 && channel.meMode() <= 7)
      lType = 2 * (channel.meMode() - 3) + 1;
    else if (channel.meMode() == 2)
      lType = 3;
    else
      lType = 1;

    // Add new ResonanceDecayChannel to map.
    partialWidths.emplace(make_pair(prodA, prodB), ResonanceDecayChannel {
      LinearInterpolator(mMin, mMax, widthData),
      prodA, prodB, lType,
      max(mMin, particleDataPtr->mMin(prodA) + particleDataPtr->mMin(prodB))
    });
  }

  // Create new or update existing HadronWidthEntry.
  HadronWidthEntry newEntry {
    LinearInterpolator(mMin, mMax, totalWidthData),
    partialWidths,
    true
  };
  auto iter = entries.find(id);
  if (iter == entries.end()) {
    entries.emplace(id, newEntry);

    // If particle is new, insert it in signature index.
    int signature = getSignature(entry->isBaryon(), entry->chargeType(id));
    auto signatureIter = signatureToParticles.find(signature);
    if (signatureIter == signatureToParticles.end())
      signatureToParticles.emplace(signature, vector<int> { id });
    else
      signatureIter->second.push_back(id);
  }
  else
    entries[id] = newEntry;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Helper function to calculate CM momentum.

static double pCMS(double eCM, double mA, double mB) {
  if (eCM <= mA + mB) return 0;
  double sCM = eCM * eCM;
  return sqrt((sCM - pow2(mA + mB)) * (sCM - pow2(mA - mB))) / (2. * eCM);
}

//--------------------------------------------------------------------------

// Get total available phase space, integrating over mass-dependent widths.

double HadronWidths::psSize(double eCM, ParticleDataEntryPtr prodA,
  ParticleDataEntryPtr prodB, double lType) const {

  // Store some important values.
  int idA      = prodA->id(), idB = prodB->id();
  double m0A   = prodA->m0(), m0B = prodB->m0();
  double mMinA = prodA->mMin(), mMinB = prodB->mMin();
  double mMaxA = prodA->mMax(), mMaxB = prodB->mMax();
  bool varA = mMaxA > mMinA, varB = mMaxB > mMinB;

  if (eCM < mMinA + mMinB)
    return 0.;

  double result;
  bool success = true;

  // No resonances.
  if (!varA && !varB)
    return pow(pCMS(eCM, m0A, m0B), lType);

  // A is resonance.
  else if (varA && !varB) {
    if (eCM <= mMinA + m0B)
      return 0.;

    // Integrate mass of A.
    auto f = [=](double mA) {
      return pow(pCMS(eCM, mA, m0B), lType) * mDistr(idA, mA); };
    if (!integrateGauss(result, f, mMinA, min(mMaxA, eCM - m0B)))
      success = false;
  }

  // B is resonance.
  else if (!varA && varB) {
    if (eCM <= m0A + mMinB)
      return 0.;

    // Integrate mass of B.
    auto f = [=](double mB) {
      return pow(pCMS(eCM, m0A, mB), lType) * mDistr(idB, mB); };
    if (!integrateGauss(result, f, mMinB, min(mMaxB, eCM - m0A)))
      success = false;
  }

  // Both are resonances.
  else {
    if (eCM <= mMinA + mMinB)
      return 0.;

    // Define integrand of outer integral.
    auto I = [=, &success](double mA) {

      // Define integrand of inner integral.
      auto f = [=](double mB) {
        return pow(pCMS(eCM, mA, mB), lType)
              * mDistr(idA, mA) * mDistr(idB, mB); };
      double res;

      // Integrate mass of B.
      if (!integrateGauss(res, f, mMinB, min(mMaxB, eCM - mA)))
        success = false;

      return res;
    };

    // Integrate mass of A.
    if (!integrateGauss(result, I, mMinA, min(mMaxA, eCM - mMinB)))
      success = false;
  }

  // Return result if successful.
  if (success)
    return result;
  else {
    loggerPtr->ERROR_MSG("unable to integrate");
    return numeric_limits<double>::quiet_NaN();
  }
}

//--------------------------------------------------------------------------

// Regenerate parameterization for all particles.

void HadronWidths::parameterizeAll(int precision) {

  // Get all particles with varWidth from particle database.
  vector<ParticleDataEntryPtr> variableWidthEntries;
  for (auto& mapEntry : *particleDataPtr) {
    ParticleDataEntryPtr entry = mapEntry.second;
    if (entry->varWidth())
      variableWidthEntries.push_back(entry);
  }

  // Clear existing data and parameterize new data.
  entries.clear();

  for (ParticleDataEntryPtr entry : variableWidthEntries) {
    loggerPtr->INFO_MSG("parameterizing resonance",
      to_string(entry->id()), true);
    if (!parameterizeRecursive(entry->id(), precision)) {
      loggerPtr->ABORT_MSG("parameterization failed");
      return;
    }
  }
}

//--------------------------------------------------------------------------

// Write all widths data to an xml file.

bool HadronWidths::save(ostream& stream) const {

  if (!stream.good())
    return false;

  stream << "\n";

  for (auto& mapEntry : entries) {
    int id = mapEntry.first;
    const HadronWidthEntry& entry = mapEntry.second;

    // Counter for number of entries on current line, maximum 8 per line.
    int c = 0;

    // Write total width.
    stream << "<width id=\"" << id << "\" "
           << "left=\"" << entry.width.left() << "\" "
           << "right=\"" << entry.width.right() << "\" "
           << "data=\" \n";
    for (double dataPoint : entry.width.data()) {
      stream << " " << dataPoint;
      if (++c >= 7) {
        c = 0;
        stream << " \n";
      }
    }
    stream << "\"/> \n \n";

    // Write partial widths.
    for (auto& channelEntry : entry.decayChannels) {
      const ResonanceDecayChannel& channel = channelEntry.second;
      stream << "<partialWidth id=\"" << id << "\" "
        << "products=\"" << channel.prodA << " " << channel.prodB << "\" "
        << "lType=\"" << channel.lType << "\" data=\" \n";
      c = 0;
      for (double dataPoint : channel.partialWidth.data()){
        stream << " " << dataPoint;
        if (++c >= 7) {
          c = 0;
          stream << " \n";
        }
      }
      stream << "\"/> \n \n";
    }

    stream << " \n \n";
  }

  // Done.
  return true;
}

//==========================================================================

}
