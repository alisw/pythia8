// main487.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>
// largely based on code by Marius Utheim

// Keywords: cosmic ray cascade; switch beam; switch collision energy

// Cosmic ray evolution in the atmosphere, with a direct comparison
// between Angantyr and PythiaCascade. The atmosphere is a mixture of
// nitrogen, oxygen and argon, with exponential density fall-off.
// Hadron-nucleus cross sections are taken from PythiaCascade only,
// since Angantyr ones are not available on the fly.
// Also, since Angantyr does not generate elastic cross sections, an
// arbitrary 0.9 reduction factor is applied to its total cross section.
// By default hadron interactions below 10.5 GeV CM energy are not
// considered, since Angantyr cannot go below 10 GeV. (PythiaCascade
// could go lower, opening for a possible switchover at low energies.)

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
constexpr double kg_m3_to_GeV_c2mm1mb1 =  5.60958865e26 // kg to GeV/c2
  * 0.001  // m^-1 to mm^-1
  * 1e-31; // m^-2 to mb^-1

// Medium parameters. Assume exponential atmosphere.
constexpr double mp   = 0.93827;
constexpr double mAir = 0.9315;
constexpr double mediumDensity = 1.225 * kg_m3_to_GeV_c2mm1mb1 / mAir;
constexpr double rhoAir = 1.225e-4; // g cm^-2 mm^-1
constexpr double mediumHeight = 100 * KM2MM;
constexpr double H = 10.4 * KM2MM;

// Atmosphere with 78% 14N, 21% 16O and 1% 40Ar by number of nuclei (mol).
constexpr double fracN  = 0.78;
constexpr double fracO  = 0.21;
constexpr double fracAr = 0.01;

// Model parameters.
// The zenith angle of the primary particle (0 = straight down).
constexpr double zenithAngle = 0.;

// Get the medium depth of a particle at the specified height.
double getDepth(double h) {
  return rhoAir * H * exp(-h / H) /  cos(zenithAngle); }

//==========================================================================

int main() {

  // Energy of primary incident proton (in GeV).
  double pPri = 1e6;

  // Number of events per case. Only do a few since each shower is so big.
  int nEvent = 100;

  // Minimal collision CM-frame energies allowed in the shower evolution,
  // restricted by Angantyr for these comnpåarisons.
  double eCMMin = 10.5;

  // Minimal kinetic energy for histogramming of particles.
  double eKinMin = 10.;

  // Set maximum size on the event record, to limit runaway code.
  int maxSize = 2000000;

  // Book histograms.
  const int nCases = 2;
  double depthMax = rhoAir * H / cos(zenithAngle);
  Hist nInt[nCases], diffHad[nCases], diffMuon[nCases], prodnue[nCases],
       prodnumu[nCases];
  for (int iC = 0; iC < nCases; ++iC) {
    nInt[iC].book("depth of interactions",            100, 0., depthMax);
    diffHad[iC].book("hadron production-decay depth", 100, 0., depthMax);
    diffMuon[iC].book("muon production-decay depth",  100, 0., depthMax);
    prodnue[iC].book("nu_e production depth",         100, 0., depthMax);
    prodnumu[iC].book("nu_mu production depth",       100, 0., depthMax);
  }

  // Four-momentum of incoming initiator.
  double pxPri = 0.;
  double pyPri = pPri * sin(zenithAngle);
  double pzPri = pPri * cos(zenithAngle);
  Vec4 p0(pxPri, pyPri, -pzPri, sqrt(pow2(pPri) + pow2(mp)));

  // Commonly reused variables.
  bool   canDecayNow, mustDecayNow;
  int    idNow, idTarg, Ztarg, Atarg, angCount[3] = {0,0,0};
  double mNow, eNow, sigmaN, sigmaO, sigmaAr, sigmaPick, sigmaNow,
         zNow, dzds, logR, zNext, eCMNow;
  Vec4   pNow, vDec, vInt, dirNow;
  clock_t timeBeg, timeEnd;
  double timeSec[2];

  // Pythia object for storing and parsing the full atmospheric cascade.
  // Also for handling Angantyr particle decays.
  Pythia pythiaMain;
  Event& eventMain = pythiaMain.event;
  // Prepare to do decays but no hard processes.
  pythiaMain.readString("ProcessLevel:all = off");
  pythiaMain.readString("211:mayDecay = on");
  pythiaMain.readString("13:mayDecay  = on");
  pythiaMain.readString("321:mayDecay = on");
  pythiaMain.readString("130:mayDecay = on");

  // Return error if pythiaMain initialization fails.
  if (!pythiaMain.init()) return 1;

  // Pythia object for performing individual Angantyr collisions.
  Pythia pythiaAng;
  Event& eventAng = pythiaAng.event;
  // Enable Angantyr. Reuse MPI and cross sections initialization files.
  // If they don't exist, you can generate them by running main424.
  pythiaAng.readFile("main424.cmnd");
  pythiaAng.readString("HeavyIon:mode = 2");
  // Set up for fixed-target collisions.
  pythiaAng.readString("Beams:frameType = 3");
  pythiaAng.settings.parm("Beams:pzA", -pPri);
  pythiaAng.readString("Beams:pzB = 0.");
  // Unit weight needed in cascade.
  pythiaAng.readString("HeavyIon:forceUnitWeight = on");
  // Decays to be done by pythiaMain.
  pythiaAng.readString("HadronLevel:Decay = off");
  // Reduce printout and relax energy-momentum conservation.
  pythiaAng.readString("Print:quiet = on");
  pythiaAng.readString("Check:epTolErr = 0.1");

  // Return error if pythiaAng initialization fails.
  if (!pythiaAng.init()) return 1;

  // PythiaCas object for performing individual PythiaCascade collisions.
  PythiaCascade pythiaCas;

  // Return error if PythiaCas initialization fails.
  // See PythiaCascade.h for explanation of reuseMPI options. In this run
  // it is sensible but not necessary to use same data file as for Angantyr.
  bool listFinal   = false;
  bool rapidDecays = false;
  int reuseMPI     = -1;
  string reuseFile = (reuseMPI < 0) ? "main424.cmnd" : "pythiaCascade.mpi";
  if ( !pythiaCas.init( pPri + mp, listFinal, rapidDecays, 0.,
    reuseMPI, reuseFile) ) return 1;

  // Begin loop over Angantyr (iCase = 0) and PythiaCas (iCase = 1).
  for (int iCase = 0; iCase < 2; ++iCase) {
    timeBeg = clock();

    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Insert primary particle in cleared main event record.
      eventMain.clear();
      eventMain.append(90,  -11, 0, 0, 1, 1, 0, 0, p0, mp);
      eventMain.append(2212, 12, 0, 0, 0, 0, 0, 0, p0, mp);

      // Set initial location of initiator, where z is distance above ground.
      double heightNow = mediumHeight;
      eventMain[0].yProd(-heightNow * tan(zenithAngle));
      eventMain[0].zProd(heightNow);
      eventMain[1].yProd(-heightNow * tan(zenithAngle));
      eventMain[1].zProd(heightNow);

      // Loop over particles (usually hadrons) in the main event record.
      for (int iHad = 1; iHad < eventMain.size(); ++iHad) {
        Particle& hadNow = eventMain[iHad];

        // Skip already fragmented/decayed or upwards-moving particles.
        if (!hadNow.isFinal() || hadNow.pz() > 0.) continue;

        //  Projectile properties. Invariant mass  with a p/n nucleon at rest.
        idNow        = hadNow.id();
        pNow         = hadNow.p();
        mNow         = hadNow.m();
        eNow         = hadNow.e();
        mustDecayNow = false;
        eCMNow       = (pNow + Vec4(0, 0, 0, mp)).mCalc();

        // Find decay vertex for unstable hadrons. (Below ground if no decay.)
        vDec = hadNow.canDecay() ? hadNow.vDec() : Vec4( 0., 0., -1., 0.);
        canDecayNow = hadNow.canDecay() && vDec.pz() > 0.;

        // Low energy hadrons should not interact with medium.
        // Decay non-hadrons or low-energy ones if decay happens above ground.
        if (!hadNow.isHadron() || eCMNow < eCMMin) {
          if (canDecayNow) mustDecayNow = true;
          else continue;
        }

        // Get hadron-nucleus cross sections and pick interaction vertex.
        // This part is temporary until Angantyr cross sections are reliable.
        if (!mustDecayNow) {

          // Set up hadron-nucleon cross section. Fail if vanishing.
          if (!pythiaCas.sigmaSetuphN( idNow, pNow, mNow)) {
            if (canDecayNow) mustDecayNow = true;
            else continue;

          } else {
            // Get hadron-ion cross sections, corrected by fractions.
            sigmaN  = fracN  * pythiaCas.sigmahA(14);
            sigmaO  = fracO  * pythiaCas.sigmahA(16);
            sigmaAr = fracAr * pythiaCas.sigmahA(40);

            // Pick nucleus by relative cross sections.
            sigmaPick = (sigmaN + sigmaO + sigmaAr) * pythiaCas.rndm().flat();
            if (sigmaPick < sigmaN)               idTarg = 1000070140;
            else if (sigmaPick < sigmaN + sigmaO) idTarg = 1000080160;
            else                                  idTarg = 1000180400;
            Ztarg = (idTarg / 10000) % 100;
            Atarg = (idTarg / 10) % 1000;

            // Atmosphere medium density is nucleons per volume (not nuclei),
            // so need cross section per nucleon.
            sigmaNow = (sigmaN + sigmaO + sigmaAr)
                     / (fracN * 14. + fracO * 16. + fracAr * 40.);

            // Ad hoc Angantyr cross section reduction from excluding elastic
            // collisions handing in Angantyr, while included in PythiaCascade.
            if (iCase == 0) sigmaNow *= 0.9;

            // Calculate potential interaction vertex.
            dirNow = pNow / pNow.pAbs();
            zNow  = hadNow.zProd();
            dzds  = hadNow.pz() / hadNow.pAbs();
            logR  = log(pythiaCas.rndm().flat());
            zNext = -H * log( exp(-zNow / H)
                  + dzds / (H * sigmaNow * mediumDensity) * logR );
            vInt  = hadNow.vProd() + (zNext - zNow) * dirNow / dzds;

            // Done if hadron reaches surface before both interaction
            // and decay.
            if (vInt.pz() < 0. && !canDecayNow) continue;

            // Do decay if it happens first.
            if (vDec.pz() > vInt.pz()) mustDecayNow = true;
          }
        }

        // Handling of decay or interaction with Angantyr.
        if (iCase == 0) {

          // Perform the decay if it happens first.
          if (mustDecayNow) {
            // Decay incoming particle. Skip particle on a failure, which
            // could lead to unintended punch-through of a particle.
            if (!pythiaMain.moreDecays(iHad)) continue;

          // Perform interaction if it happens first.
          } else {
            if ( pythiaAng.setBeamIDs(idNow, idTarg)
            && pythiaAng.setKinematics(pNow, Vec4(0., 0., 0., 0.)) ) {
            } else {
              ++angCount[0];
              continue;
            }

            // Simulate interaction. Skip particle if failure, which
            // could lead to unintended punch-through of a particle.
            if (!pythiaAng.next()) {
              ++angCount[1];
              continue;
            }
            ++angCount[2];

            // Append target.
            int iTarg = eventMain.append(idTarg, -181, iHad, iHad, 0, 0, 0, 0,
                0., 0., 0., mp, mp);
            eventMain[iTarg].vProdAdd(vInt);

            // Copy final-state particles.
            int sizeOld = eventMain.size();
            for (int iSub = 3; iSub < eventAng.size(); ++iSub) {
              if (eventAng[iSub].isFinal()) {
                int iNew = eventMain.append(eventAng[iSub]);
                eventMain[iNew].mothers(iHad, iTarg);
                eventMain[iNew].vProdAdd(vInt);
              }
            }

            // Update daughters of the interacting particles.
            int sizeNew = eventMain.size();
            eventMain[iTarg].daughters(sizeOld, sizeNew - 1);
            eventMain[iHad].daughters(sizeOld, sizeNew - 1);
            eventMain[iHad].statusNeg();
            eventMain[iHad].tau(( vInt.e() - eventMain[iHad].tProd() ) * mNow
              / eNow);
          }

        // Handling of decay or interaction with PythiaCascade.
        } else {

          // Do decay or collision. Empty event means failure, so skip,
          // which could lead to unintended punch-through of a particle.
          int sizeSave    = eventMain.size();
          Event& eventTmp = (mustDecayNow)
                          ? pythiaCas.nextDecay(idNow, pNow, mNow, vDec)
                          : pythiaCas.nextColl( Ztarg, Atarg, vInt);
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
        }

        // Stop generation if event record extremely long.
        if (eventMain.size() > maxSize) break;

      // End of loop over interactions + decays inside a single cascade.
      }

      // Begin analysis. Loop over all particles to find interaction depths.
      for (Particle& h : eventMain) {
        if (h.status() == -12) continue;
        double depthProd = getDepth(h.zProd());
        double depthDec  = getDepth(h.isFinal() ? 0. : h.zDec());

        // If particle came from the medium, record the interaction depth.
        if (h.status() == -181) {
          nInt[iCase].fill(depthProd);
        }
        // Otherwise, track depths where particles are created/destroyed.
        else if (h.e() - h.m() > eKinMin) {
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

    // End loops over events and over Angantyr/PythiaCascade cases.
    }
    timeEnd = clock();
    timeSec[iCase] = double(timeEnd - timeBeg) / double(CLOCKS_PER_SEC);
  }

  // Print time.
  cout << " \n Runtime for Angantyr events is " << fixed
       << setprecision(3) << timeSec[0] << " s" << endl
       << " Runtime for PythiaCascade events is " << timeSec[1] << " s"
       << endl;

  // Print Angantyr errors.
  cout << "\n Angantyr failed setup collision = " << angCount[0] << endl
       << " Angantyr failed performing collision = " << angCount[1] << endl
       << " Angantyr successful collision = " << angCount[2] << endl;

  // Print statistics, mainly for errors.
  cout << "\n Statistics from PythiaMain: " << endl;
  pythiaMain.stat();
  cout << "\n Statistics from PythiaAng: " << endl;
  pythiaAng.stat();
  cout << "\n Statistics from PythiaCas: " << endl;
  pythiaCas.stat();

  // Book histograms.
  Hist nHad[nCases], nMuon[nCases], nnue[nCases], nnumu[nCases];
  for (int iC = 0; iC < nCases; ++iC) {
    nHad[iC] .book("", 100, 0., depthMax);
    nMuon[iC].book("", 100, 0., depthMax);
    nnue[iC] .book("", 100, 0., depthMax);
    nnumu[iC].book("", 100, 0., depthMax);

    // Integrate production minus depletion to find particle number by depth.
    double nHadSum = 0., nMuonSum = 0., nnueSum = 0., nnumuSum = 0.;
    for (int i = 1; i <= 100; ++i) {
      double depthNow = depthMax * (i - 0.5) / 100.;
      if (depthNow > getDepth(0.)) break;
      nHadSum  += diffHad[iC] .getBinContent(i);
      nMuonSum += diffMuon[iC].getBinContent(i);
      nnueSum  += prodnue[iC] .getBinContent(i);
      nnumuSum += prodnumu[iC].getBinContent(i);
      nHad[iC] .fill(depthNow, nHadSum );
      nMuon[iC].fill(depthNow, nMuonSum);
      nnue[iC] .fill(depthNow, nnueSum );
      nnumu[iC].fill(depthNow, nnumuSum);
    }

    // Normalize histograms.
    nInt[iC].normalizeSpectrum(nEvent);
    nHad[iC]  /= nEvent;
    nMuon[iC] /= nEvent;
    nnue[iC]  /= nEvent;
    nnumu[iC] /= nEvent;
  }

  // Normalize and plot histograms.
  HistPlot plt("plot487");
  vector<string> col = { "r", "b" };
  vector<string> legend = {  "Angantyr", "PythiaCascade"};
  plt.frame("fig487", "Atmospheric depth of nuclear interactions",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) dn_{int}/dX$", 6.4, 4.8);
  for (int iC = 0; iC < nCases; ++iC)
    plt.add(nInt[iC], "-,"+col[iC], legend[iC]);
  plt.plot( 0., depthMax, 0.02, 20., true);
  plt.frame("", "Number of hadrons at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{had}$", 6.4, 4.8);
  for (int iC = 0; iC < nCases; ++iC)
    plt.add(nHad[iC], "-,"+col[iC], legend[iC]);
  plt.plot( 0., depthMax, 5., 5000., true);
  plt.frame("", "Number of muons at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{\\mu}$", 6.4, 4.8);
  for (int iC = 0; iC < nCases; ++iC)
    plt.add(nMuon[iC], "-,"+col[iC], legend[iC]);
  plt.plot( 0., depthMax, 5., 5000., true);
  plt.frame("", "Number of e neutrinos at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{\\nu_e}$", 6.4, 4.8);
  for (int iC = 0; iC < nCases; ++iC)
    plt.add(nnue[iC], "-,"+col[iC], legend[iC]);
  plt.plot( 0., depthMax, 0.2, 200., true);
  plt.frame("", "Number of mu neutrinos at depth",
    "$X$ (g/cm$^2$)", "$(1/n_{ev}) \\int_0^{X} dn_{\\nu_\\mu}$", 6.4, 4.8);
  for (int iC = 0; iC < nCases; ++iC)
    plt.add(nnumu[iC], "-,"+col[iC], legend[iC]);
  plt.plot( 0., depthMax, 2., 2000., true);

}
