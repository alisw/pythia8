// main367.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: top; R-hadrons;

// Authors: Torbjörn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Fictitious scenario where the top is long-lived enough to hadronize.
// This is studied using the R-hadron machinery, in e+e- or pp collisions.
// Kindly disregard warning that wzDecayShowers find undecayed tops.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Key parameters: e+e- or pp, # events, t mass.
  bool doee    = true;
  int nEvent   = 1000;
  double mTop  = 171.;

  // Generator. Shorthand for the event and particle data.
  Pythia pythia;
  Event& event = pythia.event;
  ParticleData& pData = pythia.particleData;

  // Set up beams: e+e- or pp. Switch off ISR for former.
  if (doee) {
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM = 500.");
    pythia.readString("PDF:lepton = off");
  } else {
    pythia.readString("Beams:eCM = 13000.");
  }

  // Set up top pair production.
  if (doee) {
    pythia.readString("Top:ffbar2ttbar(s:gmZ) = on");
  } else {
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.readString("PhaseSpace:pTHatMin = 100.");
  }

  // Set top mass and width.
  pData.m0( 6, mTop);
  pData.mWidth( 6, 0.1);
  pData.doForceWidth( 6, true);

  // Allow R-hadron formation with top marked as stop.
  pythia.readString("Rhadrons:allow = on");
  pythia.readString("RHadrons:idStop = 6");

  // Set R-hadrons stable.
  //pythia.readString("RHadrons:allowDecay = off");

  // Switch off key components. Useful for first checks, but not full run.
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("HadronLevel:Hadronize = off");

  // Possibility to switch off particle data and event listings.
  //pythia.readString("Init:showChangedSettings = off");
  //pythia.readString("Init:showChangedParticleData = off");
  //pythia.readString("Next:numberShowInfo = 0");
  //pythia.readString("Next:numberShowProcess = 0");
  //pythia.readString("Next:numberShowEvent = 0");

  // Debug options.
  pythia.readString("Check:nErrList = 5");

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Histograms.
  Hist nChargedH("charged multiplicity", 100, -0.5, 799.5);
  Hist dndyChargedH("dn/dy charged", 100, -10., 10.);
  Hist dndyRH("dn/dy R-hadrons", 100, -5., 5.);
  Hist pTRH("pT R-hadrons", 100, 0., 1000.);
  Hist xRH("p_RHadron / p_sparticle", 100, 0.9, 1.1);
  Hist mDiff("m(Rhadron) - m(sparticle)", 100, 0., 5.);

  // R-hadron flavour composition.
  map<int, int> flavours;

  // Begin event loop. Generate events. Skip. if failure.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Loop over final charged particles in the event.
    // The R-hadrons may not yet have decayed here.
    int nCharged = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal()) {
        pSum += event[i].p();
        if (event[i].isCharged()) {
          ++nCharged;
          dndyChargedH.fill( event[i].y() );
        }
      }
    }
    nChargedH.fill( nCharged );

    // Loop over final R-hadrons in the event: kinematic distribution
    for (int i = 0; i < event.size(); ++i) {
      int idAbs = event[i].idAbs();
      if (idAbs > 1000100 && idAbs < 2000000 && idAbs != 1009002) {
        ++flavours[ event[i].id() ];
        dndyRH.fill( event[i].y() );
        pTRH.fill( event[i].pT() );
        // Trace back to mother; compare momenta and masses.
        int iMother = i;
        while( event[iMother].statusAbs() > 100)
          iMother = event[iMother].mother1();
        double xFrac = event[i].pAbs() / event[iMother].pAbs();
        xRH.fill( xFrac);
        double mShift = event[i].m() - event[iMother].m();
        mDiff.fill( mShift );

      // End of loop over final R-hadrons.
      }
    }

  // End of event loop.
  }

  // Final statistics, flavour composition and histogram output.
  pythia.stat();
  cout << "\n Composition of produced R-hadrons \n    code            "
       << "name   times " << endl;
  for (map<int, int>::iterator flavNow = flavours.begin();
    flavNow != flavours.end(); ++flavNow)  cout << setw(8)
    << flavNow->first << setw(16) << pythia.particleData.name(flavNow->first)
    << setw(8) << flavNow->second << endl;
  cout << nChargedH << dndyChargedH << dndyRH << pTRH << xRH << mDiff;

  // Done.
  return 0;
}
