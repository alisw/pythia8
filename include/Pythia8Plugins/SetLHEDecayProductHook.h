// SetLHEDecayProductHook.h is part of the PYTHIA event generator.
// Copyright (C) 2024 Stephen Mrenna, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Stephen Mrenna, December 2022.

// This class is used to modify the resonance decay products from an
// input LHE file.

// The following settings are available:
// flag SetLHEDecayProduct:filter = false
//   Activate filter if true.

#ifndef Pythia8_SetLHEDecayProductHooks_H
#define Pythia8_SetLHEDecayProductHooks_H

// Includes.
#include "Pythia8/Pythia.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/Event.h"

namespace Pythia8 {

class SetLHEDecayProductHook : public UserHooks {
public:

  // Constructor.
  SetLHEDecayProductHook(Settings &settings, const ParticleData* pdtPtrIn);

  // Override base class methods.
  bool canVetoProcessLevel() override {return true;}
  bool doVetoProcessLevel(Event& process) override {
    return checkVetoProcessLevel(process);}
  bool initAfterBeams() override;

  // Class specific.
  bool checkVetoProcessLevel(Event& process);
  unsigned long int returnCounter() {return counter;};

private:

  // Return the custom particle mass.
  double m0powheg(const int idIn) {
    double mReturn(-1.0);
    switch( idIn ) {
    case 3: mReturn=0.2; break;
    case 4: mReturn=1.5; break;
    case 1: mReturn=0.1; break;
    case 2: mReturn=0.1; break;
    default: mReturn = pdtPtr->m0(idIn); break;
    }
    return mReturn;
  }

  // Data members.
  bool filter;
  const ParticleData* pdtPtr;
  unsigned long int counter;

};

//--------------------------------------------------------------------------

// Constructor.

SetLHEDecayProductHook::SetLHEDecayProductHook(Settings &settings,
  const ParticleData* pdtPtrIn) :
  pdtPtr(pdtPtrIn), counter(0) {
  settings.addFlag("SetLHEDecayProduct:filter", false);
}

//--------------------------------------------------------------------------

// Intialize the user hook after the beams.

bool SetLHEDecayProductHook::initAfterBeams() {
  filter = settingsPtr->flag("SetLHEDecayProduct:filter");
  return true;
}

//--------------------------------------------------------------------------

// Return true if the resonance decays are vetoed.

bool SetLHEDecayProductHook::checkVetoProcessLevel(Event& process) {

  if (!filter) return false;
  counter++;
  // Determine which W decays hadronically
  int hadronW = ( rndmPtr->flat() < 0.5 ) ? 1 : 2;
  int idQuark = ( rndmPtr->flat() < 0.5 ) ? 1 : 3;
  int newCol  = process.nextColTag();

  int countW(0);
  for (int i = 0; i < process.size(); ++i) {
    if (process[i].idAbs() == 24 && process[i].status() == -22) {
      countW++;
      if( countW > 2 ) break;
      int idFermion = ( hadronW == countW ) ? idQuark : 11 ;
      // Identify W+- daughters, properly ordered.
      int idV = process[i].id();
      int i1  = process[i].daughter1();
      int i2  = process[i].daughter2();
      int id1  = process[i1].idAbs();
      int id2  = process[i2].idAbs();
      if( id1 == idFermion || id2 == idFermion ) continue;

      double mV = process[i].m();
      Vec4 p1 = process[i1].p();
      p1.bstback( process[i].p() );
      double pMag = p1.pAbs();

      // The current distributions are for the particle with the
      // same charge sign as the mother, i.e. W- -> e-.
      if (process[i1].id() * idV > 0) swap( i1, i2);
      process[i1].id(idFermion * process[i1].id()/abs(process[i1].id()) );
      process[i2].id((idFermion+1) * process[i2].id()/abs(process[i2].id()) );
      if( idFermion != 11 ) {
        if( process[i1].id() > 0 ) {
          process[i1].col(newCol);
          process[i2].acol(newCol);

        } else {
          process[i2].col(newCol);
          process[i1].acol(newCol);
        }
      }

      double m1 = m0powheg(process[i1].idAbs());
      double m2 = m0powheg(process[i2].idAbs());

      // Energy and absolute momentum of first decay product in W rest frame.
      double e1 = 0.5* (pow2(mV) + pow2(m1) - pow2(m2))/mV;
      double pA = sqrt(pow2(e1) - pow2(m1));

      p1.rescale3( pA/ pMag );
      p1.e( e1 );
      Vec4 p2   = Vec4(0,0,0,mV) - p1;


      p1.bst( process[i].p() );
      p2.bst( process[i].p() );
      process[i1].p( p1 );
      process[i2].p( p2 );
      process[i1].m( m1 );
      process[i2].m( m2 );
    }
    // End of loop over W's. Do not veto any events.

  }
  return false;
}

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_SetLHEDecayProductHooks_H
