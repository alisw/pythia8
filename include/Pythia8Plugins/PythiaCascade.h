// PythiaCascade.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Torbjorn Sjostrand.

#ifndef Pythia8_PythiaCascade_H
#define Pythia8_PythiaCascade_H

#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// Wrapper class for Pythia usage in cosmic ray or detector cascades.
// Here it is assumed that the user wants to control the full
// evolution.  The complete event record is then kept somewhere
// separate from PYTHIA, and specific input on production and decay
// generates the next subevent.

// Intended flow:

// - init sets up all generation, given a maximum energy. This may be
//   time-consuming, less so if MPI initialization data can be reused.

// - sigmaSetuphN prepares a possible collision for a given hadron by
//   calculating the hadron-nucleon cross section (if possible).
//   Moderately time-consuming.

// - sigmahA uses the hN cross section calculated by sigaSetuphN and
//   returns hadron-ion cross section for a collision with a specified
//   nucleus.  It can be called several times to cover a mix of target
//   nuclei, with minimal time usage.

// - nextColl performs the hadron-nucleus collision, as a sequences of
//   hadron-nucleon ones. Can be quite time-consuming.

// - nextDecay can be used anytime to decay a particle. Each
//   individual decay is rather fast, but there may be many of them.

// - stat can be used at end of run to give a summary of error
//   messages.  Negligible time usage.

// - references to particleData() and rndm() can be used in the main
// - program.

// Note: when a hadron interacts with a medium particle, the latter is
// added to the event record. This program uses these additional
// status codes for target particles in the medium:
//  181: the first (or only) target nucleon in a collision.
//  182: nucleons from further subcollisions with the same target nucleus.

// Warning: the large boosts involved at the higher cosmic ray
// energies are not perfectly managed, which leads to non-negligible
// energy-momentum non-conservation. The worst (sub)events are caught
// and regenerated.

//--------------------------------------------------------------------------

class PythiaCascade {

public:

  // Default constructor, all setup is done in init().
  PythiaCascade() = default;

  //--------------------------------------------------------------------------

  // Initialize PythiaCascade for a given maximal incoming energy.

  // Set eMax to the maximal incoming energy (in GeV) on fixed target.

  // Keep listFinal = false if you want to get back the full event
  // record, else only the "final" particles of the collision or decay
  // are returned.

  // Keep rapidDecays = false if you want to do every decay yourself,
  // else all particle types with a tau0 below smallTau0 will be
  // decayed immediately.  The tau0 units are mm/c, so default gives c
  // * tau0 = 1e-10 mm = 100 fm.  Note that time dilation effects are
  // not included, and they can be large.

  // Set reuseMPI = false in case you cannot reuse the
  // pythiaCascade.mpi file, or any other relevant stored file. The
  // saved initialization must have been done with the same or a
  // larger eMax than the one now being used.

  void init(double eMaxIn = 1e9, bool listFinalIn = false,
    bool rapidDecaysIn = false, double smallTau0In = 1e-10,
    bool reuseMPI = true, string initFile = "pythiaCascade.mpi") {

    // Store input for future usage.
    eMax        = eMaxIn;
    listFinal   = listFinalIn;
    rapidDecays = rapidDecaysIn;
    smallTau0   = smallTau0In;

    // Proton mass.
    mp      = pythiaMain.particleData.m0(2212);

    // Main Pythia object for managing the cascade evolution in a
    // nucleus. Can also do decays, but no hard processes.
    pythiaMain.readString("ProcessLevel:all = off");
    pythiaMain.readString("13:mayDecay  = on");
    pythiaMain.readString("211:mayDecay = on");
    pythiaMain.readString("321:mayDecay = on");
    pythiaMain.readString("130:mayDecay = on");
    pythiaMain.settings.flag("ParticleDecays:limitTau0", rapidDecays);
    pythiaMain.settings.parm("ParticleDecays:tau0Max", smallTau0);

    // Redure statistics printout to relevant ones.
    pythiaMain.readString("Stat:showProcessLevel = off");
    pythiaMain.readString("Stat:showPartonLevel = off");

    // Initialize.
    pythiaMain.init();

    // Secondary Pythia object for performing individual collisions,
    // or decays. Variable incoming beam type and energy.
    pythiaColl.readString("Beams:allowVariableEnergy = on");
    pythiaColl.readString("Beams:allowIDAswitch = on");

    // Set up for fixed-target collisions.
    pythiaColl.readString("Beams:frameType = 3");
    pythiaColl.settings.parm("Beams:pzA", eMax);
    pythiaColl.readString("Beams:pzB = 0.");

    // Must use the soft and low-energy QCD processes.
    pythiaColl.readString("SoftQCD:all = on");
    pythiaColl.readString("LowEnergyQCD:all = on");

    // Primary (single) decay to be done by pythiaColl, to circumvent
    // limitTau0.
    pythiaColl.readString("13:mayDecay  = on");
    pythiaColl.readString("211:mayDecay = on");
    pythiaColl.readString("321:mayDecay = on");
    pythiaColl.readString("130:mayDecay = on");

    // Secondary decays to be done by pythiaMain, respecting limitTau0.
    pythiaColl.readString("HadronLevel:Decay = off");

    // Reduce printout and relax energy-momentum conservation.
    // (Unusually large errors unfortunate consequence of large boosts.)
    pythiaColl.readString("Print:quiet = on");
    pythiaColl.readString("Check:epTolErr = 0.01");
    pythiaColl.readString("Check:epTolWarn = 0.0001");
    pythiaColl.readString("Check:mTolErr = 0.01");

    // Redure statistics printout to relevant ones.
    pythiaColl.readString("Stat:showProcessLevel = off");
    pythiaColl.readString("Stat:showPartonLevel = off");

    // Reuse MPI initialization file if it exists; else create a new one.
    if (reuseMPI)
         pythiaColl.readString("MultipartonInteractions:reuseInit = 3");
    else pythiaColl.readString("MultipartonInteractions:reuseInit = 1");
    pythiaColl.settings.word("MultipartonInteractions:initFile", initFile);

    // Initialize.
    pythiaColl.init();

  }

  //--------------------------------------------------------------------------

  // Calculate the average number of collisions. Average number of
  // hadron-nucleon collisions in a hadron-nucleus one. If not in
  // table then interpolate, knowing that <n> - 1 propto A^{2/3}.

  double nCollAvg(int A) {
    for (int i = 0; i < nA; ++i) {
      if (A == tabA[i]) {
        return (sigmaNow < tabBorder) ? 1. + tabSlopeLo[i] * sigmaNow
          : 1. + tabOffset[i] + tabSlope[i] * sigmaNow;
      } else if (A < tabA[i]) {
        double nColl1 = (sigmaNow < tabBorder) ? tabSlopeLo[i - 1] * sigmaNow
          : tabOffset[i - 1] + tabSlope[i - 1] * sigmaNow;
        double nColl2 = (sigmaNow < tabBorder) ? tabSlopeLo[i] * sigmaNow
          : tabOffset[i] + tabSlope[i] * sigmaNow;
        double wt1 = (tabA[i] - A) / (tabA[i] - tabA[i - 1]);
        return 1. + wt1 * pow( A / tabA[i - 1], 2./3.) * nColl1
          + (1. - wt1) * pow( A / tabA[i], 2./3.) * nColl2;
      }
    }

    return numeric_limits<double>::quiet_NaN();
  }

  //--------------------------------------------------------------------------

  // Calculate the hadron-proton collision cross section.
  // Return false if not possible to find.

  bool sigmaSetuphN(int idNowIn, Vec4 pNowIn, double mNowIn) {

    // Cannot handle low-energy hadrons.
    if (pNowIn.e() - mNowIn < eKinMin) return false;

    // Cannot handle hadrons above maximum energy set at initialization.
    if (pNowIn.e() > eMax) {
      logger.ERROR_MSG("too high energy");
      return false;
    }

    // Save incoming quantities for reuse in later methods.
    idNow = idNowIn;
    pNow  = pNowIn;
    mNow  = mNowIn;

    // Calculate hadron-nucleon cross section. Check if cross section
    // vanishes.
    eCMNow = (pNow + Vec4(0, 0, 0, mp)).mCalc();
    sigmaNow = pythiaColl.getSigmaTotal(idNow, 2212, eCMNow, mNow, mp);
    if (sigmaNow <= 0.) {
      if (eCMNow - mNow - mp > eKinMin)
        logger.ERROR_MSG("vanishing cross section");
       return false;
    }

    // Done.
    return true;

  }

  //--------------------------------------------------------------------------

  // Calculate the hadron-nucleus cross section for a given nucleon
  // number A, using hN cross section from sigmaSetuphN. Interpolate
  // where not (yet) available.

  double sigmahA(int A) {

    // Restrict to allowed range 1 <= A <= 208.
    if (A < 1 || A > 208) {
      logger.ERROR_MSG("A is outside of valid range (1 <= A <= 208)");
      return 0.;
    }

    // Correction factor for number of h-nucleon collisions per
    // h-nucleus one.
    double sigmahA = A * sigmaNow / nCollAvg(A);

    // Done.
    return sigmahA;

  }

  //--------------------------------------------------------------------------

  // Generate a collision, and return the event record.
  // Input (Z, A) of nucleus, and optionally collision vertex.

  Event& nextColl(int Znow, int Anow, Vec4 vNow = Vec4() ) {

    // References to the two event records. Clear main event record.
    Event& eventMain = pythiaMain.event;
    Event& eventColl = pythiaColl.event;
    eventMain.clear();

    // Restrict to allowed range 1 <= A <= 208.
    if (Anow < 1 || Anow > 208) {
      logger.ERROR_MSG("A is outside of valid range (1 <= A <= 208)");
      return eventMain;
    }

    // Insert incoming particle in cleared main event record.
    eventMain.append(90,   -11, 0, 0, 1, 1, 0, 0, pNow, mNow);
    int iHad = eventMain.append(idNow, 12, 0, 0, 0, 0, 0, 0, pNow, mNow);
    eventMain[iHad].vProd(vNow);

    // Set up for collisions on a nucleus.
    int np      = Znow;
    int nn      = Anow - Znow;
    int sizeOld = 0;
    int sizeNew = 0;
    Vec4 dirNow = pNow / pNow.pAbs();
    Rndm& rndm  = pythiaMain.rndm;

    // Drop rate of geometric series. (Deuterium has corrected nCollAvg.)
    double probMore = 1. - 1. / nCollAvg(Anow);

    // Loop over varying number of hit nucleons in target nucleus.
    for (int iColl = 1; iColl <= Anow; ++iColl) {
      if (iColl > 1 && rndm.flat() > probMore) break;

      // Pick incoming projectile: trivial for first subcollision, else ...
      int iProj    = iHad;
      int procType = 0;

      // ... find highest-pLongitudinal particle from latest subcollision.
      if (iColl > 1) {
        iProj = 0;
        double pMax = 0.;
        for (int i = sizeOld; i < sizeNew; ++i)
        if ( eventMain[i].isFinal() && eventMain[i].isHadron()) {
          double pp = dot3(dirNow, eventMain[i].p());
          if (pp > pMax) {
            iProj = i;
            pMax  = pp;
          }
        }

        // No further subcollision if no particle with enough energy.
        if ( iProj == 0 || eventMain[iProj].e() - eventMain[iProj].m()
          < eKinMin) break;

        // Choose process; only SD or ND at perturbative energies.
        double eCMSub = (eventMain[iProj].p() + Vec4(0, 0, 0, mp)).mCalc();
        if (eCMSub > 10.) procType = (rndm.flat() < probSD) ? 4 : 1;
      }

      // Pick one p or n from target.
      int idProj = eventMain[iProj].id();
      bool doProton = rndm.flat() < (np / double(np + nn));
      if (doProton) np -= 1;
      else          nn -= 1;
      int idNuc = (doProton) ? 2212 : 2112;

      // Do a projectile-nucleon subcollision. Return empty event if failure.
      pythiaColl.setBeamIDs(idProj, idNuc);
      pythiaColl.setKinematics(eventMain[iProj].p(), Vec4());
      if (!pythiaColl.next(procType)) {
        eventMain.clear();
        return eventMain;
      }

      // Insert target nucleon. Mothers are (0,iProj) to mark who it
      // interacted with. Always use proton mass for simplicity.
      int statusNuc = (iColl == 1) ? -181 : -182;
      int iNuc = eventMain.append( idNuc, statusNuc, 0, iProj, 0, 0, 0, 0,
        0., 0., 0., mp, mp);
      eventMain[iNuc].vProdAdd(vNow);

      // Update full energy of the event with the proton mass.
      eventMain[0].e( eventMain[0].e() + mp);
      eventMain[0].m( eventMain[0].p().mCalc() );

      // Insert secondary produced particles (but skip intermediate partons)
      // into main event record and shift to correct production vertex.
      sizeOld = eventMain.size();
      for (int iSub = 3; iSub < eventColl.size(); ++iSub) {
        if (!eventColl[iSub].isFinal()) continue;
        int iNew = eventMain.append(eventColl[iSub]);
        eventMain[iNew].mothers(iNuc, iProj);
        eventMain[iNew].vProdAdd(vNow);
      }
      sizeNew = eventMain.size();

      // Update daughters of colliding hadrons and other history.
      eventMain[iProj].daughters(sizeOld, sizeNew - 1);
      eventMain[iNuc].daughters(sizeOld, sizeNew - 1);
      eventMain[iProj].statusNeg();
      eventMain[iProj].tau(0.);

    // End of loop over interactions in a nucleus.
    }

    // Optionally do decays of short-lived particles.
    if (rapidDecays) pythiaMain.moreDecays();

    // Optionally compress event record.
    if (listFinal) compress();

    // Return generated collision.
    return eventMain;

  }

  //--------------------------------------------------------------------------

  // Generate a particle decay, and return the event record.
  // You can allow sequential decays, if they occur rapidly enough.

  Event& nextDecay(int idNowIn, Vec4 pNowIn, double mNowIn,
    Vec4 vNow = Vec4() ) {

    // Save incoming quantities. (Not needed, but by analogy with
    // collisions.)
    idNow = idNowIn;
    pNow  = pNowIn;
    mNow  = mNowIn;

    // References to the event records. Clear them.
    Event& eventMain = pythiaMain.event;
    Event& eventColl = pythiaColl.event;
    eventMain.clear();
    eventColl.clear();

    // Insert incoming particle in cleared collision event record.
    eventColl.append(90,   -11, 0, 0, 1, 1, 0, 0, pNow, mNow);
    int iHad = eventColl.append(idNow, 12, 0, 0, 0, 0, 0, 0, pNow, mNow);
    eventColl[iHad].vProd(vNow);

    // Decay incoming particle. Return empty event if fail. Copy event record.
    if (!pythiaColl.moreDecays(iHad)) return eventMain;
    eventMain = eventColl;

    // Optionally do secondary decays of short-lived particles.
    if (rapidDecays) pythiaMain.moreDecays();

    // Optionally compress event record.
    if (listFinal) compress();

    // Return generated collision.
    return eventMain;

  }

  //--------------------------------------------------------------------------

  // Compress the event record by removing initial and intermediate
  // particles.  Keep line 0, since the += operator for event records
  // only adds from 1 on.

  void compress() {

    // Reference to the main event record. Original and new size.
    Event& eventMain = pythiaMain.event;
    int sizeOld = eventMain.size();
    int sizeNew = 1;

    // Loop through all particles and move up the final ones to the top.
    // Remove history information. Update event record size.
    for (int i = 1; i < sizeOld; ++i) if (eventMain[i].isFinal()) {
      eventMain[sizeNew] = eventMain[i];
      eventMain[sizeNew].mothers( 0, 0);
      eventMain[sizeNew].daughters( 0, 0);
      ++sizeNew;
    }

    // Shrink event record to new size.
    eventMain.popBack( sizeOld - sizeNew);

  }

  //--------------------------------------------------------------------------

  // Summary of aborts, errors and warnings.

  void stat() {
    pythiaMain.stat();
    pythiaColl.stat();
    logger.errorStatistics();
  }

  //--------------------------------------------------------------------------

  // Possibility to access particle data and random numbers from
  // pythiaMain.

  ParticleData& particleData() {return pythiaMain.particleData;}

  Rndm& rndm() {return pythiaMain.rndm;}

//--------------------------------------------------------------------------

private:

  // The Pythia instances used for decays and for collisions. Could
  // be made public, but cleaner to allow only limited access as
  // above.
  Pythia pythiaMain, pythiaColl;

  // Logger instance for errors in this class.
  Logger logger;

  // Save quantities.
  bool   listFinal, rapidDecays;
  int    idNow;
  double eMax, smallTau0, mp, mNow, eCMNow, sigmaNow;
  Vec4   pNow;

  // Fraction of single diffractive events beyond first collision in nucleus.
  static constexpr double probSD = 0.3;

  // Hadrons below the kinetic energy threshold cannot interact with medium.
  // (At least not using Pythia.) Used both in lab and in CM frame.
  static constexpr double eKinMin = 0.2;

  // Studied nuclei by A number, with offset and slope of <nColl>(sigma):
  // 1H, 2H, 4He, 9Be, 12C, 14N, 16O, 27Al, 40Ar, 56Fe, 63Cu, 84Kr,
  // 107Ag, 129Xe, 197Au, 208Pb.
  static constexpr int nA = 16;
  static constexpr double tabBorder = 20.;
  static const double tabA[nA];
  static const double tabOffset[nA];
  static const double tabSlope[nA];
  static const double tabSlopeLo[nA];

};

// Constant definitions.
const double PythiaCascade::tabA[] = {1, 2, 4, 9, 12, 14, 16, 27, 40, 56,
  63, 84, 107, 129, 197, 208};
const double PythiaCascade::tabOffset[] = {0., 0.03, 0.08, 0.15, 0.20, 0.20,
  0.20, 0.26, 0.30, 0.34, 0.40, 0.40, 0.40, 0.50, 0.50, 0.60};
const double PythiaCascade::tabSlope[] = {0., 0.0016, 0.0033, 0.0075, 0.0092,
  0.0105, 0.012, 0.017, 0.022, 0.027, 0.028, 0.034, 0.040, 0.044, 0.055,
  0.055};
const double PythiaCascade::tabSlopeLo[] = {0., 0.0031, 0.0073, 0.015, 0.0192,
  0.0205, 0.022, 0.03, 0.037, 0.044, 0.048, 0.054, 0.06, 0.069, 0.08, 0.085};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_PythiaCascade_H
