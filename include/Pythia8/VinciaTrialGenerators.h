// VinciaTrialGenerators.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef Pythia8_VinciaTrial_H
#define Pythia8_VinciaTrial_H

// Pythia headers.
#include "Pythia8/Basics.h"
#include "Pythia8/PythiaStdlib.h"

// Vincia headers.
#include "Pythia8/VinciaCommon.h"

namespace Pythia8 {

// Helpful enums.
enum class TrialGenType { Void = 0, FF = 1, RF = 2, IF = 3, II = 4 };
// (Default is used for soft, global, or splittings as appropriate.)
enum class BranchType { Void = -1, Emit = 0, SplitF = 1, SplitI = 2,
  Conv = 3 };
enum class Sector { Void = -99, ColI = -1, Default = 0, ColK = 1 };

// Forward declarations.
class ZetaGenerator;
class ZetaGeneratorSet;

//==========================================================================

// Helper struct for passing trial-alphaS information.

struct EvolutionWindow {

  int runMode{};
  double alphaSmax{}, b0{}, kMu2{}, lambda2{}, qMin{};
  map<int, double> mass;

};

//==========================================================================

// Base class for trial generators.

class TrialGenerator {

 public:

  // Main constructor.
  TrialGenerator(bool isSectorIn, TrialGenType trialGenTypeIn,
    BranchType branchTypeIn, ZetaGeneratorSet* zetaGenSet)
    : isSector(isSectorIn), trialGenTypeSav(trialGenTypeIn),
        branchType(branchTypeIn) { setupZetaGens(zetaGenSet); }

  // Destructor.
  virtual ~TrialGenerator() = default;

  // Set pointers to zetaGenerators.
  void setupZetaGens(ZetaGeneratorSet* zetaGenSet);

  // Re-calculate the current zeta limits and integrals.
  virtual void reset(double Q2min, double s, const vector<double> & masses,
    enum AntFunType antFunType, double xA = 1., double xB = 1.);

  // Generate the next trial scale.
  virtual double genQ2(double Q2MaxNow, Rndm* rndmPtr,
    const EvolutionWindow* evWindowPtrIn, double colFac,
    double wtIn, Logger* loggerPtr, int verboseIn);

  // Get the invariants.
  virtual bool genInvariants(double sAnt, const vector<double>& masses,
    vector<double>& invariants, Rndm* rndmPtr, Logger* loggerPtr,
    int verboseIn);

  // Calculate the trial based on invariants and saved quantities.
  virtual double aTrial(vector<double>& invariants,
    const vector<double>& masses, int verboseIn);

  // Calculate the colour and coupling stripped antenna function.
  virtual double aTrialStrip(vector<double>& invariants,
    const vector<double>& masses, int verboseIn);

  // Delete the current trial.
  virtual void resetTrial();

  // Mark trial as used.
  virtual void needsNewTrial();

  // Return the sector.
  int getSector() {return (int)sectorSav;}

 protected:

  // Calculate the Kallen factor.
  virtual void calcKallenFac(double, const vector<double>&) {
    kallenFacSav = 1.0;}

  // Calculate the PDF ratio.
  virtual void calcRpdf(const vector<double>&) {Rpdf = 1.0;}

  void addGenerator(ZetaGeneratorSet* zetaGenSet,
    Sector sector = Sector::Default);

  // True when init succeeds.
  bool isInit{false};

  // Information set at construction.
  const bool isSector;
  const TrialGenType trialGenTypeSav;
  const BranchType branchType;

  // Common prefactors to the trial integral.
  double kallenFacSav{1.};
  double Rpdf{1.};

  // Information about the antenna.
  double sAntSav{};
  vector<double> massesSav;

  // Information about the trial.
  bool hasTrial{false};
  double q2Sav{}, colFacSav{};
  const EvolutionWindow* evWindowSav{};
  Sector sectorSav;

  // Map from sector to the correct zeta generator.
  // (note these live inside a ZetaGeneratorSet)
  map<Sector, ZetaGeneratorPtr> zetaGenPtrs;

  // Map from sector to the corresponding zeta phase-space limits.
  map<Sector, pair<double, double>> zetaLimits;

  // Save the zeta integrals.
  map<Sector, double> IzSav;

  // Save which sectors are currently active.
  map<Sector, bool> isActiveSector;

};

//==========================================================================

// Trial generator for final-final branchings.

class TrialGeneratorFF : public TrialGenerator {

 public:

  // Default constructor/destructor.
  TrialGeneratorFF(bool isSectorIn, BranchType branchTypeIn,
    ZetaGeneratorSet* zetaGenSet) : TrialGenerator(isSectorIn,
      TrialGenType::FF, branchTypeIn, zetaGenSet) {;}
  ~TrialGeneratorFF() = default;

 private:

  void calcKallenFac(double sIK, const vector<double>& masses);

};

//==========================================================================

// Trial generator for resonance-final branchings.

class TrialGeneratorRF : public TrialGenerator{

 public:

  // Default constructor/destructor.
  TrialGeneratorRF(bool isSectorIn, BranchType branchTypeIn,
    ZetaGeneratorSet* zetaGenSet) : TrialGenerator(isSectorIn,
      TrialGenType::RF, branchTypeIn, zetaGenSet) {;}
  ~TrialGeneratorRF() = default;

 private:

  void calcKallenFac(double sAK, const vector<double>& masses);

};

//==========================================================================

// Trial generator for initial-final branchings.

class TrialGeneratorIF : public TrialGenerator {

 public:

  // Default constructor/destructor.
  TrialGeneratorIF(bool isSectorIn, BranchType branchTypeIn,
    ZetaGeneratorSet* zetaGenSet) : TrialGenerator(isSectorIn,
      TrialGenType::IF, branchTypeIn, zetaGenSet) {;}
  ~TrialGeneratorIF() = default;

};

//==========================================================================

// Trial generator for initial-initial branchings.

class TrialGeneratorII : public TrialGenerator {

 public:

  // Default constructor/destructor.
  TrialGeneratorII(bool isSectorIn, BranchType branchTypeIn,
    ZetaGeneratorSet* zetaGenSet) : TrialGenerator(isSectorIn,
      TrialGenType::II, branchTypeIn, zetaGenSet) {;}
  ~TrialGeneratorII() = default;

};

//==========================================================================

// Place to store all types of zeta trial generators.
// To live in VinicaFSR, VinciaISR.

class ZetaGeneratorSet {

 public:

  // Construct all zeta generators for a given type.
  ZetaGeneratorSet(TrialGenType trialGenTypeIn);

  // Destructor.
  ~ZetaGeneratorSet() = default;

  // Get ptr to ZetaGenerator for a sector.
  ZetaGeneratorPtr getZetaGenPtr(BranchType branchType, Sector sectIn);

  TrialGenType getTrialGenType() {return trialGenTypeSav;}

 protected :

  const TrialGenType trialGenTypeSav;

  void addGenerator(ZetaGeneratorPtr zGenPtr);

  map<pair<BranchType, Sector>, ZetaGeneratorPtr> zetaGenPtrs;

};

//==========================================================================

// Base class for zeta trial generators.

class ZetaGenerator {

 public:

  // Constructor and destructor.
  ZetaGenerator(TrialGenType trialGenTypeIn, BranchType branchTypeIn,
    Sector sectorIn, double globalIn) : trialGenType(trialGenTypeIn),
    branchType(branchTypeIn), sector(sectorIn), globalFactSav(globalIn) {;}
  virtual ~ZetaGenerator() = default;

  // Get (best/physical) limits given a set of input parameters.
  virtual double getzMin(double Q2min,double sAnt,
    const vector<double>& masses, double xA = 1., double xB = 1.) = 0;
  virtual double getzMax(double Q2min,double sAnt,
    const vector<double>& masses, double xA = 1., double xB = 1.) = 0;

  // Get hull of physical phase space in zeta.
  virtual double getzMinHull(double Q2min,double sAnt,
    const vector<double>& masses, double xA = 1., double xB = 1.) {
    return getzMin(Q2min, sAnt, masses, xA, xB);}
  virtual double getzMaxHull(double Q2min,double sAnt,
    const vector<double>& masses, double xA = 1., double xB = 1.) {
    return getzMax(Q2min, sAnt, masses, xA, xB);}

  // Get constant factor for zeta integral.
  // NOTE: only used in II conversion trial.
  virtual double getConstFactor(double,
    const vector<double>&) {return 1.;}

  // Set the invariants for the current value of the evolution variables.
  virtual void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) = 0;

  // Evaluate the trial antenna given invariants and masses.
  virtual double aTrial(const vector<double>& invariants,
    const vector<double>& masses) = 0;

  // Check if this trial is active for specific AntFunType.
  virtual bool isActive(enum AntFunType) {return false;}

  // Return information about this generator.
  TrialGenType getTrialGenType() {return trialGenType;}
  Sector getSector() {return sector;}
  BranchType getBranchType() {return branchType;}

  // Return multiplier to convert to global.
  double globalMultiplier() {return globalFactSav;}

  // The zeta integral.
  // Optionally with exponent gamma for PDF overestimate.
  double getIz(double zMinIn, double zMaxIn, double gammaPDF = 1.) {
    return zetaIntSingleLim(zMaxIn, gammaPDF)
      -zetaIntSingleLim(zMinIn, gammaPDF);}

  // Generate a value of zeta.
  double genZeta(Rndm* rndmPtr, double zMinIn, double zMaxIn,
    double gammaPDF = 1.);

  // Print the trial generator.
  void print();

 protected:

  // The functional form of the zeta integral.
  // Optionally with exponent gamma for PDF overestimate.
  virtual double zetaIntSingleLim(double z, double gammaPDF = 1.) = 0;

  // The function form of the inverse of the zeta integral.
  // Optionally with exponent gamma for PDF overestimate.
  virtual double inverseZetaIntegral(double Iz, double gammaPDF = 1.) = 0;

  // Check if invariants are valid.
  bool valid(const string& method, Logger* loggerPtr, int verbose, double zIn);
  bool valid(const string& method, Logger* loggerPtr, int verbose, double zIn,
    const double& Q2In);

  // Labels to define this trial generator (set in derived constructors).
  const TrialGenType trialGenType{TrialGenType::Void};
  const BranchType branchType{BranchType::Void};
  const Sector sector{Sector::Void};

  // Multiplier to convert trial to global.
  const double globalFactSav;

};
//==========================================================================

// Final-final trial generators.

//==========================================================================

// The final-final default sector generator.

class ZGenFFEmitSoft : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenFFEmitSoft() : ZetaGenerator(TrialGenType::FF , BranchType::Emit,
    Sector::Default, 1.0) {;}
  ~ZGenFFEmitSoft() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QQEmitFF || antFunType == QGEmitFF ||
      antFunType == GQEmitFF || antFunType == GGEmitFF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The final-final ColI sector emission generator.

class ZGenFFEmitColI: public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenFFEmitColI() : ZetaGenerator(TrialGenType::FF, BranchType::Emit,
    Sector::ColI,1.0) {;}
  ~ZGenFFEmitColI() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt,const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt,const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == GQEmitFF || antFunType == GGEmitFF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The final-final ColK sector emission generator.

class ZGenFFEmitColK : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenFFEmitColK() : ZetaGenerator(TrialGenType::FF, BranchType::Emit,
    Sector::ColK, 1.0) {;}
  ~ZGenFFEmitColK() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt,const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt,const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QGEmitFF || antFunType == GGEmitFF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The final-final default sector splitting generator.

class ZGenFFSplit : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenFFSplit() : ZetaGenerator(TrialGenType::FF , BranchType::SplitF,
    Sector::Default, 0.5) {;}
  ~ZGenFFSplit() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == GXSplitFF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// Resonance-final trial generators.

//==========================================================================

// The resonance-final default sector generator.

class ZGenRFEmitSoft : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenRFEmitSoft() : ZetaGenerator(TrialGenType::RF, BranchType::Emit,
    Sector::Default, 1.0) {;}
  ~ZGenRFEmitSoft() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QQEmitRF || antFunType == QGEmitRF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;
};

//==========================================================================

// The resonance-final default sector alternate generator.

class ZGenRFEmitSoftAlt : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenRFEmitSoftAlt() : ZetaGenerator(TrialGenType::RF, BranchType::Emit,
    Sector::Default, 1.0) {;}
  ~ZGenRFEmitSoftAlt() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA=1., double xB=1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA=1., double xB=1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses ) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QQEmitRF || antFunType == QGEmitRF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The resonance-final ColK sector generator.

class ZGenRFEmitColK : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenRFEmitColK() : ZetaGenerator(TrialGenType::RF, BranchType::Emit,
    Sector::ColK, 1.0) {;}
  ~ZGenRFEmitColK() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial( const vector<double>& invariants,
    const vector<double>& masses ) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QGEmitRF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The resonance-final default sector splitting generator.

class ZGenRFSplit : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenRFSplit() : ZetaGenerator(TrialGenType::RF, BranchType::SplitF,
    Sector::Default, 0.5) {;}
  ~ZGenRFSplit() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses ) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == XGSplitRF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// Intial-final trial generators.

//==========================================================================

// The initial-final default sector generator.

class ZGenIFEmitSoft : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIFEmitSoft() :
    ZetaGenerator(TrialGenType::IF, BranchType::Emit, Sector::Default, 1.0) {;}
  ~ZGenIFEmitSoft() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {return
      antFunType == QQEmitIF || antFunType == QGEmitIF ||
      antFunType == GQEmitIF || antFunType == GGEmitIF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-final ColI sector generator.

class ZGenIFEmitColA : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIFEmitColA() :
    ZetaGenerator(TrialGenType::IF, BranchType::Emit, Sector::ColI, 1.0) {;}
  ~ZGenIFEmitColA() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType  == GQEmitIF || antFunType == GGEmitIF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-final ColK sector generator.

class ZGenIFEmitColK : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIFEmitColK() : ZetaGenerator(TrialGenType::IF, BranchType::Emit,
    Sector::ColK, 1.0) {;}
  ~ZGenIFEmitColK() = default;

  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA=1., double xB=1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA=1., double xB=1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {return
      antFunType == QGEmitIF || antFunType == GGEmitIF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-final initial antenna splitting generator.

class ZGenIFSplitA: public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIFSplitA() : ZetaGenerator(TrialGenType::IF, BranchType::SplitI,
    Sector::Default, 1.) {;}
  ~ZGenIFSplitA() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial( const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QXConvIF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-final final antenna splitting generator.

class ZGenIFSplitK : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIFSplitK() : ZetaGenerator(TrialGenType::IF, BranchType::SplitF,
    Sector::Default, .5) {;}
  ~ZGenIFSplitK() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double> & invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == XGSplitIF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-final splitting generator.

class ZGenIFConv : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIFConv() : ZetaGenerator(TrialGenType::IF, BranchType::Conv,
    Sector::Default, 1.) {;}
  ~ZGenIFConv() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == GXConvIF;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-initial trial generators.

//==========================================================================

// The initial-initial default sector generator.

class ZGenIIEmitSoft : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIIEmitSoft() : ZetaGenerator(TrialGenType::II, BranchType::Emit,
    Sector::Default, 1.0) {;}
  ~ZGenIIEmitSoft() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QQEmitII || antFunType == GQEmitII ||
      antFunType == GGEmitII;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-initial ColI sector generator.

class ZGenIIEmitCol : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIIEmitCol() : ZetaGenerator(TrialGenType::II, BranchType::Emit,
    Sector::ColI, 1.0) {;}
  ~ZGenIIEmitCol() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses ) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType  == GQEmitII || antFunType == GGEmitII;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-initial initial splitting generator.

class ZGenIISplit : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIISplit() : ZetaGenerator(TrialGenType::II, BranchType::SplitI,
    Sector::Default, 1.) {;}
  ~ZGenIISplit() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == QXConvII;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

//==========================================================================

// The initial-initial splitting generator.

class ZGenIIConv : public ZetaGenerator {

 public:

  // Constructor/destructor.
  ZGenIIConv() : ZetaGenerator(TrialGenType::II, BranchType::Conv,
    Sector::Default, 1.) {;}
  ~ZGenIIConv() = default;

  // Overridden methods.
  double getzMin(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getzMax(double Q2,double sAnt, const vector<double>& masses,
    double xA = 1., double xB = 1.) override;
  double getConstFactor(double sAnt, const vector<double>& masses) override;
  void genInvariants(double Q2In, double zIn, double sAnt,
    const vector<double>& masses, vector<double>& invariants,
    Logger* loggerPtr, int verboseIn) override;
  double aTrial(const vector<double>& invariants,
    const vector<double>& masses) override;
  bool isActive(enum AntFunType antFunType) override {
    return antFunType == GXConvII;}

 private:

  double zetaIntSingleLim(double z, double gammaPDF = 1.) override;
  double inverseZetaIntegral(double Iz, double gammaPDF = 1.) override;

};

} // end namespace Pythia8

#endif // end Pythia8_VinciaTrial_H
