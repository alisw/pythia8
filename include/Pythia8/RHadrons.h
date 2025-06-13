// RHadrons.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains a class for the production and decay
// of long-lived heavy coloured particles, for now the gluino.

#ifndef Pythia8_RHadrons_H
#define Pythia8_RHadrons_H

#include "Pythia8/FragmentationModel.h"

namespace Pythia8 {

//==========================================================================

// The RHadrons class contains the routines for the production and decay
// of long-lived heavy coloured particles.

class RHadrons : public FragmentationModel {

public:

  // Constructor.
  RHadrons() : FragmentationModel(), allowRH(), allowRSb(), allowRSt(),
    allowRGo(), allowSomeR(), setMassesRH(), idRSb(), idRSt(), idRGo(),
    maxWidthRH(), probGluinoballRH(), mOffsetCloudRH(), mCollapseRH(),
    diquarkSpin1RH(), m0Sb(), m0St(), m0Go(), nRHad(0), iRHad(), iBef(),
    iSys(), systemPtr() {}

  // Initialize and save pointers.
  bool init(StringFlav* flavSelPtrIn = nullptr, StringPT* pTSelPtrIn = nullptr,
    StringZ* zSelPtrIn = nullptr, FragModPtr fragModPtrIn = nullptr) override;

  // Fragment the event.
  bool fragment(int iSub, ColConfig& colConfig, Event& event,
    bool isDiff = false, bool systemRecoil = true) override;

  // Decay R-hadrons.
  bool decay( Event& event);

  // Tell whether a given particle is supposed to form R-hadrons.
  bool givesRHadron(int id);

  // Tell whether any R-hadrons have been formed.
  bool exist() {return (nRHad > 0);}

  // Tell whether a R-hadron production+decay happened, and trace down.
  int trace(int i) {
    for (int iR = 0; iR < nRHad; ++iR)
      if (iBefRHad[iR] == i || iCreRHad[iR] == i) return iAftRHad[iR];
    return 0;}

private:

  // Constants: could only be changed in the code itself.
  static const int    IDRHADSB[14], IDRHADST[14], IDRHADGO[38], NTRYMAX;
  static const double MSAFETY, EGBORROWMAX;

  // Initialization data, mainly read from Settings.
  bool   allowRH, allowRSb, allowRSt, allowRGo, allowSomeR, setMassesRH;
  int    idRSb, idRSt, idRGo;
  double maxWidthRH, probGluinoballRH, mOffsetCloudRH, mCollapseRH,
         diquarkSpin1RH, m0Sb, m0St, m0Go;

  // Current event properties.
  vector<int>  iBefRHad, iCreRHad, iRHadron, iAftRHad;
  vector<bool> isTriplet;
  int          nRHad, iRHad, iBef, iSys;
  ColSinglet*  systemPtr;

  // Split a system that contains both a sparticle and a junction.
  bool splitOffJunction( ColConfig& colConfig, Event& event);

  // Open up a closed gluon/gluino loop.
  bool openClosedLoop( ColConfig& colConfig, Event& event);

  // Split a single colour singlet that contains two sparticles.
  bool splitSystem( ColConfig& colConfig, Event& event);

  // Produce a R-hadron from a squark.
  bool produceSquark( ColConfig& colConfig, Event& event);

  // Produce a R-hadron from a gluino.
  bool produceGluino( ColConfig& colConfig, Event& event);

  // Construct R-hadron code from squark and (di)quark codes.
  int toIdWithSquark( int id1, int id2);

  // Construct squark and (di)quark codes from R-hadron code.
  pair<int,int> fromIdWithSquark( int idRHad);

  // Construct R-hadron code from endpoints and a gluino.
  int toIdWithGluino( int id1, int id2);

  // Construct endpoint codes from R-hadron code with a gluino.
  pair<int,int> fromIdWithGluino( int idRHad);

  // Construct modified four-vectors to match modified masses.
  bool newKin( Vec4 pOld1, Vec4 pOld2, double mNew1, double mNew2,
    Vec4& pNew1, Vec4& pNew2, bool checkMargin = true);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_RHadrons_H
