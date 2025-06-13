// StringFragmentation.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the classes for string fragmentation.
// StringEnd: keeps track of the fragmentation step.
// StringFragmentation: is the top-level class.

#ifndef Pythia8_StringFragmentation_H
#define Pythia8_StringFragmentation_H

#include "Pythia8/FragmentationModel.h"

namespace Pythia8 {

//==========================================================================

// The StringEnd class contains the information related to
// one of the current endpoints of the string system.
// Only to be used inside StringFragmentation, so no private members.

class StringEnd {

public:

  // Constructor.
  StringEnd() : particleDataPtr(), flavSelPtr(), pTSelPtr(), zSelPtr(),
    fromPos(), thermalModel(), mT2suppression(), iEnd(), iMax(), idHad(),
    iPosOld(), iNegOld(), iPosNew(), iNegNew(), hadSoFar(), colOld(), colNew(),
    pxOld(), pyOld(), pxNew(), pyNew(), pxHad(), pyHad(), mHad(), mT2Had(),
    zHad(), GammaOld(), GammaNew(), xPosOld(), xPosNew(), xPosHad(), xNegOld(),
    xNegNew(), xNegHad(), aLund(), bLund(), iPosOldPrev(), iNegOldPrev(),
    colOldPrev(), pxOldPrev(), pyOldPrev(), GammaOldPrev(), xPosOldPrev(),
    xNegOldPrev(), mVecRatio(1.), tinyEq(), pT2tiny() {}

  // Save pointers.
  void init( ParticleData* particleDataPtrIn, StringFlav* flavSelPtrIn,
    StringPT* pTSelPtrIn, StringZ* zSelPtrIn, Settings& settings) {
    particleDataPtr = particleDataPtrIn; flavSelPtr = flavSelPtrIn;
    flavSelNow = *flavSelPtr;
    pTSelPtr = pTSelPtrIn; zSelPtr = zSelPtrIn;
    bLund = zSelPtr->bAreaLund(); aLund = zSelPtr->aAreaLund();
    thermalModel   = settings.flag("StringPT:thermalModel");
    mT2suppression = settings.flag("StringPT:mT2suppression");
    closePacking = settings.flag("ClosePacking:doClosePacking"); }

  // Set up initial endpoint values from input.
  void setUp(bool fromPosIn, int iEndIn, int idOldIn, int iMaxIn,
    double pxIn, double pyIn, double GammaIn, double xPosIn,
    double xNegIn, int colIn, double mVecRatioIn);

  // Fragment off one hadron from the string system, in flavour and pT.
  void newHadron(double kappaModifier, bool forbidPopcornNow = false,
    double strangeJunc = 0., double probQQmod = 1.);

  // Fragment off one hadron from the string system, in momentum space,
  // by taking steps either from positive or from negative end.
  Vec4 kinematicsHadron(StringSystem& system, StringVertex& newVertex,
    double zHadIn);

  // Generate momentum for some possible next hadron, based on mean values
  // to get an estimate for rapidity and pT.
  Vec4 kinematicsHadronTmp(StringSystem system, Vec4 pRem, double phi,
    double mult);

  // Update string end information after a hadron has been removed.
  void update();
  void storePrev();
  void updateToPrev();

  // Constants: could only be changed in the code itself.
  static const double TINY, PT2SAME, MEANMMIN, MEANM, MEANPT;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointers to classes for flavour, pT and z generation.
  StringFlav*   flavSelPtr;
  StringPT*     pTSelPtr;
  StringZ*      zSelPtr;

  // Local copy of flavSelPtr for modified flavour selection.
  StringFlav    flavSelNow;

  // Data members.
  bool   fromPos, thermalModel, mT2suppression, closePacking;
  int    iEnd, iMax, idHad, iPosOld, iNegOld, iPosNew, iNegNew, hadSoFar,
         colOld, colNew;
  double pxOld, pyOld, pxNew, pyNew, pxHad, pyHad, mHad, mT2Had, zHad,
         GammaOld, GammaNew, xPosOld, xPosNew, xPosHad, xNegOld, xNegNew,
         xNegHad, aLund, bLund;
  int    iPosOldPrev, iNegOldPrev, colOldPrev;
  double pxOldPrev, pyOldPrev, GammaOldPrev, xPosOldPrev, xNegOldPrev,
         mVecRatio, tinyEq, pT2tiny;
  FlavContainer flavOld, flavNew, flavOldPrev;
  Vec4   pHad, pSoFar;

};

//==========================================================================

// The StringFragmentation class contains the top-level routines
// to fragment a colour singlet partonic system.

class StringFragmentation : public FragmentationModel {

public:

  // Constructor.
  StringFragmentation() :
    FragmentationModel(), flavRopePtr(), closePacking(),
    setVertices(), constantTau(), smearOn(), traceColours(false),
    hadronVertex(), stopMass(), stopNewFlav(), stopSmear(),
    pNormJunction(), pMaxJunction(), eBothLeftJunction(),
    eMaxLeftJunction(), eMinLeftJunction(), mJoin(), bLund(),
    closePackingFluxRatio(1.), closePackingPT20(1.), pT20(),
    xySmear(), maxSmear(), maxTau(), kappaVtx(), mc(), mb(),
    hasJunction(), isClosed(), iPos(), iNeg(), nExtraJoin(),
    w2Rem(), stopMassNow(), mVecRatio(1.), closedM2max(),
    idDiquark(), legMin(), legMid() {}

  // Initialize and save pointers.
  bool init(StringFlav* flavSelPtrIn = nullptr, StringPT* pTSelPtrIn = nullptr,
    StringZ* zSelPtrIn = nullptr, FragModPtr fragModPtrIn = nullptr) override;

  // Do the fragmentation: driver routine.
  bool fragment(int iSub, ColConfig& colConfig, Event& event,
    bool isDiff = false, bool systemRecoil = true) override;

  // Local copy of flavSelPtr for modified flavour selection.
  StringFlav flavSelNow;

  // Find the boost matrix to the rest frame of a junction.
  Vec4 junctionRestFrame(const Vec4& p0, const Vec4& p1, const Vec4& p2,
    const bool angleCheck = true) const;

  // Set the vector mass ratio.
  void setMVecRatio(double mVecRatioIn) {mVecRatio = mVecRatioIn;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRYFLAV, NTRYJOIN, NSTOPMASS,
                      NTRYJNMATCH, NTRYJRFEQ, NTRYSMEAR, MAXVETOFINTWO;
  static const double FACSTOPMASS, CLOSEDM2MAX, CLOSEDM2FRAC, EXPMAX,
                      MATCHPOSNEG, M2MINJRF, EMINJRF, EEXTRAJNMATCH,
                      MDIQUARKMIN, CONVJRFEQ, CHECKPOS;

  // Pointer to flavour-composition-changing ropes.
  FragModPtr  flavRopePtr;

  // Initialization data, read from Settings.
  bool   closePacking, setVertices, constantTau, smearOn,
         traceColours, hardRemn, doStrangeJunc;
  int    hadronVertex;
  double stopMass, stopNewFlav, stopSmear, pNormJunction, pMaxJunction,
         eBothLeftJunction, eMaxLeftJunction, eMinLeftJunction,
         mJoin, bLund, closePackingFluxRatio, closePackingPT20,
         qqSupPar, qqSupAnti, pT20, xySmear, maxSmear, maxTau,
         kappaVtx, mc, mb, dampPopcorn, aRemn, bRemn, strangeJuncParm;

  // Data members.
  bool   hasJunction, isClosed;
  int    iPos, iNeg, nExtraJoin;
  double w2Rem, stopMassNow, kappaModifier, probQQmod, mVecRatio, closedM2max;
  Vec4   pSum, pRem, pJunctionHadrons;

  // UserHooks flags.
  bool doChangeFragPar = false, doVetoFrag = false;

  // List of partons in string system.
  vector<int> iParton, iPartonMinLeg, iPartonMidLeg, iPartonMax;

  // Vertex information from the fragmentation process.
  vector<StringVertex> stringVertices, legMinVertices, legMidVertices;
  StringVertex newVertex;

  // Variables used for JRF finding.
  vector<Vec4>   listJRF;
  vector<double> weightJRF;
  int    iLeg[3], idLeg[3], legEnd[3];
  double weightSum, pSumJRF, m2Leg[3];
  Vec4   pLeg[3];
  bool   lastJRF, endpoint[3];

  // Boost from/to rest frame of a junction to original frame.
  RotBstMatrix MfromJRF, MtoJRF;

  // Information on diquark created at the junction.
  int    idDiquark;

  // Fictitious opposing partons in JRF: string ends for vertex location.
  Vec4 pMinEnd, pMidEnd;

  // Temporary event record for the produced particles.
  Event hadrons;

  // Information on the system of string regions.
  StringSystem system, systemMin, systemMid;

  // Information on the two current endpoints of the fragmenting
  // system, with backup copies if there is a veto.
  StringEnd posEnd, negEnd, posEndSave, negEndSave;

  // Find region where to put first string break for closed gluon loop.
  vector<int> findFirstRegion(int iSub, const ColConfig& colConfig,
    const Event& event) const;

  // Set flavours and momentum position for initial string endpoints.
  void setStartEnds(int idPos, int idNeg, const StringSystem& systemNow,
    int legNow = 3);

  // Check remaining energy-momentum whether it is OK to continue.
  bool energyUsedUp(bool fromPos);

  // Produce the final two partons to complete the system.
  bool finalTwo(bool fromPos, const Event& event, bool usedPosJun,
    bool usedNegJun);

  // Final region information.
  Vec4 pPosFinalReg, pNegFinalReg, eXFinalReg, eYFinalReg;

  // Set hadron production points in space-time picture.
  bool setHadronVertices(Event& event);

  // Construct a special joining region for the final two hadrons.
  StringRegion finalRegion();

  // Store the hadrons in the normal event record, ordered from one end.
  void store(Event& event);

  // Fragment off two of the string legs in to a junction.
  bool fragmentToJunction(Event& event,
    vector< vector< pair<double,double> > >& rapPairs);

  // Initially considered legs from the junction.
  int legMin, legMid;

  // Functions used in JRF iterative procedure.
  bool   collinearPair(Event& event);
  bool   perturbedJRF(Event& event);
  int    updateLegs(Event& event, Vec4 vJunIn, bool juncCoM = false);
  double updateWeights(double pSmall, Vec4 vJunIn);
  void   nextParton(Event& event, int leg);

  // Join extra nearby partons when stuck.
  int extraJoin(double facExtra, Event& event);

  // Get the number of nearby strings given the energies.
  void kappaEffModifier(StringSystem& systemNow,
    StringEnd end, bool fromPos, vector<int> partonList,
    vector< vector< pair<double,double> > >& rapPairs,
    double mRem, Event& event);

  double yMax(Particle pIn, double mTiny) {
    double temp = log( ( pIn.e() + abs(pIn.pz()) ) / max( mTiny, pIn.mT()) );
    return (pIn.pz() > 0) ? temp : -temp; }

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_StringFragmentation_H
