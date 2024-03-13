// main185.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.m.utheim@jyu.fi>;
//          Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: cosmic ray cascade; switch beam; switch collision energy

// This example is based on work from Eur. Phys. J. C82 (2022) 21 and
// arXiv:2108.03481 [hep-ph]. This program illustrates how the
// PythiaCascade class can be called for an interaction or a decay,
// with properties as provided by the user.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PythiaCascade.h"

using namespace Pythia8;

//==========================================================================

// The main program tests both the collision and decay possibilities,
// using four slightly different setups to illustrate variability.

int main() {

  // Basic parameter choices. Includes beam, target and decaying particle.
  double eMax      = 1e8;
  double smallTau0 = 1e-10;
  int    nEvent    = 1000;
  int    nList     = 1;
  int    idHad     = 321;
  int    Ztarg     = 18;
  int    Atarg     = 40;
  int    idDec     = 4122;

  // Histograms for number of subcollisions or final particles.
  Hist nhA( "number of subcollisions", 40, 0.5, 40.5);
  Hist nFin("number of final decay particles", 20, 0.5, 20.5);

  // Main features of four subruns, illustrating variability.
  for (int subrun = 0; subrun < 4; ++subrun) {
    bool listFinal   = (subrun > 1);
    bool rapidDecays = (subrun%2 == 1);

    // Set up Pythia wrapper.
    // Reuse of MPI data is here implicitly accepted.
    PythiaCascade pythiaCascade;
    pythiaCascade.init( eMax, listFinal, rapidDecays, smallTau0);

    // Event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Test 1 : collision with fixed-target nucleus.
      // Set up incoming particle.
      double mHad  = pythiaCascade.particleData().m0(idHad);
      double pzHad = eMax * pythiaCascade.rndm().flat();
      double eHad  = sqrt( mHad * mHad + pzHad * pzHad);
      Vec4   pHad( 0., 0., pzHad, eHad);

      // Cross section for various target nucleons. Print for first event.
      pythiaCascade.sigmaSetuphN( idHad, pHad, mHad);
      if (subrun == 0 && iEvent < nList) {
        cout << " Cross sections for incoming K+ on various targets for E = "
             << scientific << setprecision(3) << eHad << fixed << endl;
        cout << " p : " << setw(8) << pythiaCascade.sigmahA( 1) << endl;
        cout << " N : " << setw(8) << pythiaCascade.sigmahA(14) << endl;
        cout << " O : " << setw(8) << pythiaCascade.sigmahA(16) << endl;
        cout << " Ar: " << setw(8) << pythiaCascade.sigmahA(40) << endl;
      }

      // Generate event. Skip empty = aborted event. List first event.
      Event& eventColl = pythiaCascade.nextColl( Ztarg, Atarg);
      if (eventColl.size() == 0) continue;
      if (iEvent < nList) eventColl.list();

      // Count number of subcollisions. Only meaningful if history remains.
      if (!listFinal) {
        int nSub = 0;
        for (int i = 1; i < eventColl.size(); ++i)
          if (eventColl[i].status() == -181 || eventColl[i].status() == -182)
            ++nSub;
        nhA.fill( nSub);
      }

      // Test 2: decay of given particle.
      // Set up decaying particle. Random momentum and decay vertex.
      double mDec  = pythiaCascade.particleData().m0(idDec);
      double pxDec = 10. * pythiaCascade.rndm().flat() - 5.;
      double pyDec = 10. * pythiaCascade.rndm().flat() - 5.;
      double pzDec = 1000. * pythiaCascade.rndm().flat();
      double eDec  = sqrt(pxDec*pxDec + pyDec*pyDec + pzDec*pzDec + mDec*mDec);
      Vec4   pDec( pxDec, pyDec, pzDec, eDec);
      Vec4   vDec( pythiaCascade.rndm().flat(), pythiaCascade.rndm().flat(),
                   pythiaCascade.rndm().flat(), pythiaCascade.rndm().flat());

      // Generate decay. List first event, including production vertices.
      Event& eventDec = pythiaCascade.nextDecay( idDec, pDec, mDec, vDec);
      if (eventDec.size() == 0) continue;
      if (iEvent < nList) eventDec.list( true);

      // Count final multiplicity. Only done when secondary decays included.
      if (rapidDecays) {
        int nFinal = 0;
        for (int i = 1; i < eventDec.size(); ++i)
          if (eventDec[i].isFinal()) ++nFinal;
        nFin.fill( nFinal);
      }

    // End event loop. Error statistics. End subrun loop.
    }
    pythiaCascade.stat();
  }

  // Histogram printout. Done.
  cout << nhA << nFin;
}
