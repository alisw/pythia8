// main184.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.m.utheim@jyu.fi>
//          Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: cosmic ray cascade; switch beam; switch collision energy

// This example is based on work from Eur. Phys. J. C82 (2022) 21 and
// arXiv:2108.03481 [hep-ph]. main184 should be equivalent with
// main183, but with collisions and decays now handled by the
// PythiaCascade class, which streamlines the main program.  The
// stepwise addition to the event record is also more transparent.
// Reminder: for vertex Vec4 the components are labelled (px, py, pz,
// e), but actually represent (x, y, z, t) values.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PythiaCascade.h"

using namespace Pythia8;

//==========================================================================

// Note: when a hadron interacts with a medium particle, the latter is added
// to the event record. This program uses these additional status codes for
// target particles in the medium:
//  181: the first (or only) target nucleon in a collision.
//  182: nucleons from further subcollisions with the same target nucleus.

// Conversion factor between mm and km.
constexpr double KM2MM = 1.e6;

// Conversion factor for kg/m^3 to GeV/(c^2 mm mb)
constexpr double kg_m3_to_GeV_c2mm1mb1
  =  5.60958865e26 // kg to GeV/c2
    * 0.001         // m^-1 to mm^-1
    * 1e-31;        // m^-2 to mb^-1

// Medium parameters
constexpr double mAir = 0.9315;
constexpr double mediumDensity = 1.225 * kg_m3_to_GeV_c2mm1mb1 / mAir;
constexpr double rhoAir = 1.225e-4; // g cm^-2 mm^-1
constexpr double mediumHeight = 100 * KM2MM;
constexpr double H = 10.4 * KM2MM;

//==========================================================================

// Configuration of the atmospheric model.

struct Configuration {
  // Description of the configuration, to be used when plotting.
  string legend;
  // Whether the atmosphere should be exponentially attenuated, or uniform.
  bool doExponential;
  // Whether the medium particles are clustered into ions or a p/n gas.
  bool doHeavyIon;
  // The zenith angle of the primary particle.
  double zenithAngle;
  // Hadrons below the kinetic energy threshold cannot interact with medium.
  double eKinMin;

  // Get the height above ground of the primary particle when cascade begins.
  double initialHeight() {
    return doExponential ? mediumHeight : H;
  }

  // Get the medium depth of a particle at the specified height.
  double getDepth(double h) {
    return 1. / cos(zenithAngle)
        * ( doExponential ? rhoAir * H * exp(-h / H)
                          : rhoAir * (initialHeight() - h) );
  }
};

//==========================================================================

int main() {

  // Set up four different "atmospheres" to compare.
  constexpr int nCases = 4;
  vector<Configuration> configurations = {
    { "Uniform p/n atmosphere",                false, false, 0.,       0.2 },
    { "Uniform nitrogen",                      false, true,  0.,       0.2 },
    { "Exponential nitrogen",                  true,  true,  0.,       0.2 },
    { "Exponential nitrogen at 45$^{\\circ}$", true,  true,  M_PI / 4, 0.2 },
  };
  vector<string> col = { "r", "b", "k", "g" };

  // Energy of primary incident proton (in GeV).
  double pPri = 1e6;

  // Number of events per case. Only do a few since each shower is so big.
  int nEvent = 5;

  // Set maximum size on the event record, to limit runaway code.
  int maxSize = 2000000;

  // Pythia wrapper for managing the cascade evolution and particle decays.
  // See the PythiaCascade.h file for documentation, such as init() arguments,
  // notably whether a previous *.mpi file can and should be reused.
  PythiaCascade pythiaCascade;
  Rndm& rndm = pythiaCascade.rndm();
  double mp = pythiaCascade.particleData().m0(2212);
  pythiaCascade.init(pPri + mp);

  // Event record for full cascade evolution.
  Event eventMain;
  eventMain.init( "Full cascade event", &pythiaCascade.particleData());

  // Book histograms.
  double depthMax = 1.5 * rhoAir * H;
  Hist nInt[nCases], diffHad[nCases], diffMuon[nCases], prodnue[nCases],
       prodnumu[nCases];
  for (int iCase = 0; iCase < nCases; ++iCase) {
    nInt[iCase].book("depth of interactions",            100, 0., depthMax);
    diffHad[iCase].book("hadron production-decay depth", 100, 0., depthMax);
    diffMuon[iCase].book("muon production-decay depth",  100, 0., depthMax);
    prodnue[iCase].book("nu_e production depth",         100, 0., depthMax);
    prodnumu[iCase].book("nu_mu production depth",       100, 0., depthMax);
  }

  // Begin loops over cases and events.
  for (int iCase = 0; iCase < nCases; ++iCase)
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Four-momentum of incoming initiator.
    Configuration& config = configurations[iCase];
    double pxPri = 0.;
    double pyPri = pPri * sin(config.zenithAngle);
    double pzPri = pPri * cos(config.zenithAngle);
    Vec4 p0(pxPri, pyPri, -pzPri, sqrt(pow2(pPri) + pow2(mp)));

    // Insert primary particle in cleared main event record.
    eventMain.clear();
    eventMain.append(90,  -11, 0, 0, 1, 1, 0, 0, p0, mp);
    eventMain.append(2212, 12, 0, 0, 0, 0, 0, 0, p0, mp);

    // Set initial location of initiator, where z is distance above ground.
    double heightNow = config.initialHeight();
    eventMain[0].yProd(-heightNow * tan(config.zenithAngle));
    eventMain[0].zProd(heightNow);
    eventMain[1].yProd(-heightNow * tan(config.zenithAngle));
    eventMain[1].zProd(heightNow);

    // Loop over particles (usually hadrons) in the main event record.
    for (int iHad = 1; iHad < eventMain.size(); ++iHad) {
      Particle& hadNow = eventMain[iHad];
      int idNow        = hadNow.id();
      Vec4 pNow        = hadNow.p();
      double mNow      = hadNow.m();
      double eNow      = hadNow.e();
      bool mustDecayNow    = false;

      // Skip already fragmented/decayed or upwards-moving particles.
      if (!hadNow.isFinal() || hadNow.pz() > 0.) continue;

      // Find decay vertex for unstable hadrons. (Below ground if no decay.)
      Vec4 vDec = hadNow.canDecay() ? hadNow.vDec() : Vec4( 0., 0., -1., 0.);
      bool canDecayNow = hadNow.canDecay() && vDec.pz() > 0.;

      // Low energy hadrons should not interact with medium.
      // Decay non-hadrons or low-energy ones if decay happens above ground.
      if (!hadNow.isHadron() || eNow - mNow < config.eKinMin) {
        if (canDecayNow) mustDecayNow = true;
        else continue;
      }

      // Get hadron-nucleon cross section for hydrogen or nitrogen atmosphere.
      int  Znow       = (config.doHeavyIon) ?  7 : 1;
      int  Anow       = (config.doHeavyIon) ? 14 : 1;
      double sigmaNow = 0.;
      if (!mustDecayNow) {
        // Set up hadron-nucleon cross sections for the current hadron.
        if (!pythiaCascade.sigmaSetuphN( idNow, pNow, mNow)) {
          if (canDecayNow) mustDecayNow = true;
          else continue;
        }
        else {
          // Get hadron-ion cross section. For media consisting of several
          // different atom types, this should be called once for each type.
          sigmaNow = pythiaCascade.sigmahA(Anow);

          // In this run density is nucleons per volume, not nuclei.
          sigmaNow /= Anow;

          // If the cross section vanishes, decay is the only option.
          if (sigmaNow <= 0.) {
            if (canDecayNow) mustDecayNow = true;
            else continue;
          }
        }
      }

      // Calculate potential interaction vertex, depending on medium.
      Vec4 vInt( 0., 0., -1., 0.);
      if (!mustDecayNow) {
        Vec4 dirNow = pNow / pNow.pAbs();
        // Exponential atmosphere.
        if (config.doExponential) {
          double zNow  = hadNow.zProd();
          double dzds  = hadNow.pz() / hadNow.pAbs();
          double logR  = log(rndm.flat());
          double zNext = -H * log( exp(-zNow / H)
                       + dzds / (H * sigmaNow * mediumDensity) * logR );
          vInt = hadNow.vProd() + (zNext - zNow) * dirNow / dzds;
        // Homogeneous atmosphere.
        } else {
          double freePath = rndm.exp() / (mediumDensity * sigmaNow);
          vInt = hadNow.vProd() + freePath * dirNow;
        }

        // Done if hadron reaches surface before both interaction and decay.
        if (vInt.pz() < 0. && !canDecayNow) continue;

        // Do decay if it happens first.
        if (vDec.pz() > vInt.pz()) mustDecayNow = true;
      }

      // Do decay or collision. Empty event means failure, so skip,
      // which could lead to unintended punch-through of a particle.
      int sizeSave    = eventMain.size();
      Event& eventTmp = (mustDecayNow)
                      ? pythiaCascade.nextDecay(idNow, pNow, mNow, vDec)
                      : pythiaCascade.nextColl( Znow, Anow, vInt);
      if (eventTmp.size() == 0) continue;

      // Update properties of incoming hadron.
      hadNow.statusNeg();
      hadNow.daughters( sizeSave, sizeSave);
      double dTau = ( (mustDecayNow ? vDec.e() : vInt.e())
        - hadNow.tProd() ) * mNow / eNow;
      hadNow.tau( dTau);

      // Append new event to existing and link history.
      eventMain += eventTmp;
      eventMain[sizeSave].mothers( iHad, iHad);

      // Stop generation if event record extremely long.
      if (eventMain.size() > maxSize) {
        cout << " Error: maximum event size exceeded for iCase = "
             << iCase << " and iEvent = " << iEvent << endl;
        break;
      }

    // End of loop over interactions + decays inside a single cascade.
    }

    // Begin analysis. Loop over all particles to find interaction depths.
    for (Particle& h : eventMain) {
      if (h.status() == -12) continue;
      double depthProd = config.getDepth(h.zProd());
      double depthDec  = config.getDepth(h.isFinal() ? 0. : h.zDec());

      // If particle came from the medium, record the interaction depth.
      if (h.status() == -181 || h.status() == -182) {
        nInt[iCase].fill(depthProd);
      }
      // Otherwise, track depths where particles are created/destroyed.
      else if (h.e() - h.m() > config.eKinMin) {
        if (h.isHadron()) {
          if (h.isFinal())
            diffHad[iCase].fill(depthProd, 1.);
          else {
            diffHad[iCase].fill(min(depthProd, depthDec), 1.);
            diffHad[iCase].fill(max(depthProd, depthDec), -1.);
          }
        }
        else if (h.idAbs() == 13) {
          if (h.isFinal())
            diffMuon[iCase].fill(depthProd, 1.);
          else {
            diffMuon[iCase].fill(min(depthProd, depthDec), 1.);
            diffMuon[iCase].fill(max(depthProd, depthDec), -1.);
          }
        }
        else if (h.idAbs() == 12)
          prodnue[iCase].fill(depthProd, 1.);
        else if (h.idAbs() == 14)
          prodnumu[iCase].fill(depthProd, 1.);
      }
    }

  // End loops over events and cases.
  }

  // Print statistics, mainly for errors.
  pythiaCascade.stat();

  // Book histograms.
  Hist nHad[nCases], nMuon[nCases], nnue[nCases], nnumu[nCases];
  for (int iCase = 0; iCase < nCases; ++iCase) {
    nHad[iCase] .book("", 100, 0., depthMax);
    nMuon[iCase].book("", 100, 0., depthMax);
    nnue[iCase] .book("", 100, 0., depthMax);
    nnumu[iCase].book("", 100, 0., depthMax);

    // Integrate production minus depletion to find particle number by depth.
    double nHadSum = 0., nMuonSum = 0., nnueSum = 0., nnumuSum = 0.;
    for (int i = 1; i <= 100; ++i) {
      double depthNow = depthMax * (i - 0.5) / 100.;
      if (depthNow > configurations[iCase].getDepth(0.)) break;
      nHadSum  += diffHad[iCase] .getBinContent(i);
      nMuonSum += diffMuon[iCase].getBinContent(i);
      nnueSum  += prodnue[iCase] .getBinContent(i);
      nnumuSum += prodnumu[iCase].getBinContent(i);
      nHad[iCase] .fill(depthNow, nHadSum );
      nMuon[iCase].fill(depthNow, nMuonSum);
      nnue[iCase] .fill(depthNow, nnueSum );
      nnumu[iCase].fill(depthNow, nnumuSum);
    }

    // Normalize histograms.
    nInt[iCase].normalizeSpectrum(nEvent);
    nHad[iCase]  /= nEvent;
    nMuon[iCase] /= nEvent;
    nnue[iCase]  /= nEvent;
    nnumu[iCase] /= nEvent;
  }

  // Normalize and plot histograms.
  HistPlot plt("plot184");
  plt.frame("fig184", "Atmospheric depth of p/n interactions",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) dn_{int}/dX$", 6.4, 4.8);
  for (int iCase = 0; iCase < nCases; ++iCase)
    plt.add(nInt[iCase], "-,"+col[iCase], configurations[iCase].legend);
  plt.plot(true);
  plt.frame("", "Number of hadrons at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{had}$", 6.4, 4.8);
  for (int iCase = 0; iCase < nCases; ++iCase)
    plt.add(nHad[iCase], "-,"+col[iCase], configurations[iCase].legend);
  plt.plot(true);
  plt.frame("", "Number of muons at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{\\mu}$", 6.4, 4.8);
  for (int iCase = 0; iCase < nCases; ++iCase)
    plt.add(nMuon[iCase], "-,"+col[iCase], configurations[iCase].legend);
  plt.plot(true);
  plt.frame("", "Number of e neutrinos at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{\\nu_e}$", 6.4, 4.8);
  for (int iCase = 0; iCase < nCases; ++iCase)
    plt.add(nnue[iCase], "-,"+col[iCase], configurations[iCase].legend);
  plt.plot(true);
  plt.frame("", "Number of mu neutrinos at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{\\nu_\\mu}$", 6.4, 4.8);
  for (int iCase = 0; iCase < nCases; ++iCase)
    plt.add(nnumu[iCase], "-,"+col[iCase], configurations[iCase].legend);
  plt.plot(true);

}
