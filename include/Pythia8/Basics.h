// Basics.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for basic, often-used helper classes.
// RndmEngine: base class for external random number generators.
// Rndm: random number generator.
// Vec4: simple four-vectors.
// RotBstMatrix: matrices encoding rotations and boosts of Vec4 objects.
// Hist: simple one-dimensional histograms.

#ifndef Pythia8_Basics_H
#define Pythia8_Basics_H

#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/SharedPointers.h"

namespace Pythia8 {

//==========================================================================

// Forward reference to RotBstMatrix class; needed in Vec4 class.
class RotBstMatrix;

//--------------------------------------------------------------------------

// Vec4 class.
// This class implements four-vectors, in energy-momentum space.
// (But can equally well be used to hold space-time four-vectors.)

class Vec4 {

public:

  // Constructors.
  Vec4(double xIn = 0., double yIn = 0., double zIn = 0., double tIn = 0.)
    : xx(xIn), yy(yIn), zz(zIn), tt(tIn) { }
  Vec4(const Vec4& v) : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) { }
  Vec4& operator=(const Vec4& v) { if (this != &v) { xx = v.xx; yy = v.yy;
    zz = v.zz; tt = v.tt; } return *this; }
  Vec4& operator=(double value) { xx = value; yy = value; zz = value;
    tt = value; return *this; }

  // Member functions for input.
  void reset() {xx = 0.; yy = 0.; zz = 0.; tt = 0.;}
  void p(double xIn, double yIn, double zIn, double tIn)
    {xx = xIn; yy = yIn; zz = zIn; tt = tIn;}
  void p(Vec4 pIn) {xx = pIn.xx; yy = pIn.yy; zz = pIn.zz; tt = pIn.tt;}
  void px(double xIn) {xx = xIn;}
  void py(double yIn) {yy = yIn;}
  void pz(double zIn) {zz = zIn;}
  void e(double tIn) {tt = tIn;}

  // Member functions for output.
  double px() const {return xx;}
  double py() const {return yy;}
  double pz() const {return zz;}
  double e() const {return tt;}
  double& operator[](int i) {
    switch(i) {
      case 1: return xx;
      case 2: return yy;
      case 3: return zz;
      default: return tt;
    }
  }
  double mCalc() const {double temp = tt*tt - xx*xx - yy*yy - zz*zz;
    return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
  double m2Calc() const {return tt*tt - xx*xx - yy*yy - zz*zz;}
  double pT() const {return sqrt(xx*xx + yy*yy);}
  double pT2() const {return xx*xx + yy*yy;}
  double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
  double pAbs2() const {return xx*xx + yy*yy + zz*zz;}
  double eT() const {double temp = xx*xx + yy*yy;
    return tt * sqrt( temp / (temp + zz*zz) );}
  double eT2() const {double temp = xx*xx + yy*yy;
    return tt*tt * temp / (temp + zz*zz);}
  double theta() const {return atan2(sqrt(xx*xx + yy*yy), zz);}
  double phi() const {return atan2(yy,xx);}
  double thetaXZ() const {return atan2(xx,zz);}
  double pPos() const {return tt + zz;}
  double pNeg() const {return tt - zz;}
  double rap() const {
    double txyz = (tt > 0.) ? tt : sqrt(xx*xx + yy*yy + zz*zz);
    if (zz >= txyz) return 20.;
    if (zz <= -txyz) return -20.;
    return 0.5 * log( (txyz + zz) / (txyz - zz) );}
  double eta() const {double xyz = sqrt(xx*xx + yy*yy + zz*zz);
    if (zz >= xyz) return 20.;
    if (zz <= -xyz) return -20.;
    return 0.5 * log( (xyz + zz) / (xyz - zz) );}

  // Member functions that perform operations.
  void rescale3(double fac) {xx *= fac; yy *= fac; zz *= fac;}
  void rescale4(double fac) {xx *= fac; yy *= fac; zz *= fac; tt *= fac;}
  void flip3() {xx = -xx; yy = -yy; zz = -zz;}
  void flip4() {xx = -xx; yy = -yy; zz = -zz; tt = -tt;}
  void rot(double thetaIn, double phiIn);
  void rotaxis(double phiIn, double nx, double ny, double nz);
  void rotaxis(double phiIn, const Vec4& n);
  void bst(double betaX, double betaY, double betaZ);
  void bst(double betaX, double betaY, double betaZ, double gamma);
  void bst(const Vec4& pIn);
  void bst(const Vec4& pIn, double mIn);
  void bstback(const Vec4& pIn);
  void bstback(const Vec4& pIn, double mIn);
  void rotbst(const RotBstMatrix& M);

  // Operator overloading with member functions
  inline Vec4 operator-() const {Vec4 tmp; tmp.xx = -xx; tmp.yy = -yy;
    tmp.zz = -zz; tmp.tt = -tt; return tmp;}
  inline Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz;
    tt += v.tt; return *this;}
  inline Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz;
    tt -= v.tt; return *this;}
  inline Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f;
    tt *= f; return *this;}
  inline Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f;
    tt /= f; return *this;}
  inline Vec4 operator+(const Vec4& v) const {
    Vec4 tmp = *this; return tmp += v;}
  inline Vec4 operator-(const Vec4& v) const {
    Vec4 tmp = *this; return tmp -= v;}
  inline Vec4 operator*(double f) const {
    Vec4 tmp = *this; return tmp *= f;}
  inline Vec4 operator/(double f) const {
    Vec4 tmp = *this; return tmp /= f;}
  inline double operator*(const Vec4& v) const {
    return tt*v.tt - xx*v.xx - yy*v.yy - zz*v.zz;}

  // Operator overloading with friends.
  friend Vec4 operator*(double f, const Vec4& v1);

  // Print a four-vector.
  friend ostream& operator<<(ostream&, const Vec4& v) ;

  // Check if NaN, INF, or finite.
  friend inline bool isnan(const Vec4 &v) {
    return isnan(v.tt) || isnan(v.xx) || isnan(v.yy) || isnan(v.zz);}
  friend inline bool isinf(const Vec4 &v) {
    return isinf(v.tt) || isinf(v.xx) || isinf(v.yy) || isinf(v.zz);}
  friend inline bool isfinite(const Vec4 &v) {
    return isfinite(v.tt) && isfinite(v.xx) && isfinite(v.yy)
      && isfinite(v.zz);}

  // Invariant mass and its square.
  friend double m(const Vec4& v1);
  friend double m(const Vec4& v1, const Vec4& v2);
  friend double m2(const Vec4& v1);
  friend double m2(const Vec4& v1, const Vec4& v2);
  friend double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3);
  friend double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3,
                   const Vec4& v4);

  // Scalar and cross product of 3-vector parts.
  friend double dot3(const Vec4& v1, const Vec4& v2);
  friend Vec4 cross3(const Vec4& v1, const Vec4& v2);

  // Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c).
  friend Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c);

  // theta is polar angle between v1 and v2.
  friend double theta(const Vec4& v1, const Vec4& v2);
  friend double costheta(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around z axis.
  friend double phi(const Vec4& v1, const Vec4& v2);
  friend double cosphi(const Vec4& v1, const Vec4& v2);

  // phi is azimuthal angle between v1 and v2 around n axis.
  friend double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
  friend double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

  // R is distance in cylindrical (y/eta, phi) coordinates.
  friend double RRapPhi(const Vec4& v1, const Vec4& v2);
  friend double REtaPhi(const Vec4& v1, const Vec4& v2);

  // Shift four-momenta within pair from old to new masses.
  friend bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New);

  // Create two vectors that are perpendicular to both input vectors.
  friend pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2);

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // The four-vector data members.
  double xx, yy, zz, tt;

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of Vec4 class.

// Implementation of operator overloading with friends.
inline Vec4 operator*(double f, const Vec4& v1)
  {Vec4 v = v1; return v *= f;}

// Invariant mass and its square.
double m(const Vec4& v1);
double m(const Vec4& v1, const Vec4& v2);
double m2(const Vec4& v1);
double m2(const Vec4& v1, const Vec4& v2);
double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3);
double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3,
          const Vec4& v4);

// Scalar and cross product of 3-vector parts.
double dot3(const Vec4& v1, const Vec4& v2);
Vec4 cross3(const Vec4& v1, const Vec4& v2);

// Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c).
Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c);

// theta is polar angle between v1 and v2.
double theta(const Vec4& v1, const Vec4& v2);
double costheta(const Vec4& v1, const Vec4& v2);
double costheta(double e1, double e2, double m1, double m2, double s12);

// phi is azimuthal angle between v1 and v2 around z axis.
double phi(const Vec4& v1, const Vec4& v2);
double cosphi(const Vec4& v1, const Vec4& v2);

// phi is azimuthal angle between v1 and v2 around n axis.
double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

// R is distance in cylindrical (y/eta, phi) coordinates.
double RRapPhi(const Vec4& v1, const Vec4& v2);
double REtaPhi(const Vec4& v1, const Vec4& v2);

// Print a four-vector.
ostream& operator<<(ostream&, const Vec4& v) ;

// Shift four-momenta within pair from old to new masses.
bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New);

// Create two vectors that are perpendicular to both input vectors.
pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2);

//==========================================================================

// RotBstMatrix class.
// This class implements 4 * 4 matrices that encode an arbitrary combination
// of rotations and boosts, that can be applied to Vec4 four-vectors.

class RotBstMatrix {

public:

  // Constructors.
  RotBstMatrix() : M() {for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) { M[i][j] = (i==j) ? 1. : 0.; } } }
  RotBstMatrix(const RotBstMatrix& Min) : M() {
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } }
  RotBstMatrix& operator=(const RotBstMatrix& Min) {if (this != &Min) {
    for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) {
    M[i][j] = Min.M[i][j]; } } } return *this; }

  // Member functions.
  void rot(double = 0., double = 0.);
  void rot(const Vec4& p);
  void bst(double = 0., double = 0., double = 0.);
  void bst(const Vec4&);
  void bstback(const Vec4&);
  void bst(const Vec4&, const Vec4&);
  void toCMframe(const Vec4&, const Vec4&);
  void fromCMframe(const Vec4&, const Vec4&, bool flip = false);
  void toSameVframe(const Vec4&, const Vec4&);
  void fromSameVframe(const Vec4&, const Vec4&);
  void rotbst(const RotBstMatrix&);
  void invert();
  RotBstMatrix inverse() const { RotBstMatrix tmp = *this;
    tmp.invert(); return tmp; }
  void reset();

  // Return value of matrix element.
  double value(int i, int j) { return M[i][j];}

  // Crude estimate deviation from unit matrix.
  double deviation() const;

  // Print a transformation matrix.
  friend ostream& operator<<(ostream&, const RotBstMatrix&) ;

  // Private members to be accessible from Vec4.
  friend class Vec4;

  // Multiplication.
  Vec4 operator*(Vec4 p) const { p.rotbst(*this); return p; }
  RotBstMatrix operator*(RotBstMatrix R) const { R.rotbst(*this); return R; }

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // The rotation-and-boost matrix data members.
  double M[4][4];

};

//--------------------------------------------------------------------------

// Namespace function declaration; friend of RotBstMatrix class.

// Print a transformation matrix.
ostream& operator<<(ostream&, const RotBstMatrix&) ;

// Get a RotBstMatrix to rest frame of p.
inline RotBstMatrix toCMframe(const Vec4& p) {
  RotBstMatrix tmp; tmp.bstback(p); return tmp; }

// Get a RotBstMatrix from rest frame of p.
inline RotBstMatrix fromCMframe(const Vec4& p) {
  RotBstMatrix tmp; tmp.bst(p); return tmp; }

// Get a RotBstMatrix to rest frame of p1 and p2, where p1 is along
// the z-axis.
inline RotBstMatrix toCMframe(const Vec4& p1, const Vec4& p2) {
  RotBstMatrix tmp; tmp.toCMframe(p1, p2); return tmp; }

// Get a RotBstMatrix from rest frame of p1 and p2, where p1 is
// assumed by default to be along the z-axis. The flip option
// handles the case when p1 is along the negative z-axis.
inline RotBstMatrix fromCMframe(const Vec4& p1, const Vec4& p2,
  bool flip = false) {
  RotBstMatrix tmp; tmp.fromCMframe(p1, p2, flip); return tmp; }

// Get a RotBstMatrix to rest frame of ptot where pz is along the
// z-axis and pxz is in the xz-plane with positive x.
inline RotBstMatrix toCMframe(const Vec4& ptot, const Vec4& pz,
  const Vec4 & pxz) { RotBstMatrix tmp = toCMframe(ptot);
  Vec4 pzp = tmp*pz; tmp.rot(0.0, -pzp.phi()); tmp.rot(-pzp.theta());
  tmp.rot(0.0, -(tmp*pxz).phi()); return tmp; }

// Get a RotBstMatrix from rest frame of ptot where pz is along the
// z-axis and pxz is in the xz-plane with positive x.
inline RotBstMatrix fromCMframe(const Vec4& ptot, const Vec4& pz,
  const Vec4 & pxz) { return toCMframe(ptot, pz, pxz).inverse(); }

//==========================================================================

// RndmEngine is the base class for external random number generators.
// There is only one pure virtual method, that should do the generation.

class RndmEngine {

public:

  // Destructor.
  virtual ~RndmEngine() {}

  // A virtual method, wherein the derived class method
  // generates a random number uniformly distributed between 0 and 1.
  virtual double flat() {return 1;}

};

//==========================================================================

// RndmState class.
// This class describes the configuration of a Rndm object.

struct RndmState {
  int    i97{}, j97{}, seed{0};
  long   sequence{0};
  double u[97]{}, c{}, cd{}, cm{};

  // Test whether two random states would generate the same random sequence.
  bool operator==(const RndmState& other) const;
};

//==========================================================================

// Rndm class.
// This class handles random number generation according to the
// Marsaglia-Zaman-Tsang algorithm.

class Rndm {

public:

  // Constructors.
  Rndm() : initRndm(false), stateSave(), useExternalRndm(false) { }
  Rndm(int seedIn) : initRndm(false), stateSave(), useExternalRndm(false) {
    init(seedIn);}

  // Possibility to pass in pointer for external random number generation.
  bool rndmEnginePtr( RndmEnginePtr rndmEngPtrIn);

  // Initialize, normally at construction or in first call.
  void init(int seedIn = 0) ;

  // Generate next random number uniformly between 0 and 1.
  double flat() ;

  // Generate random numbers according to exp(-x).
  double exp() ;

  // Generate random numbers according to x * exp(-x).
  double xexp() { return -log(flat() * flat()) ;}

  // Generate random numbers according to exp(-x^2/2).
  double gauss() {return sqrt(-2. * log(flat())) * cos(M_PI * flat());}

  // Generate two random numbers according to exp(-x^2/2-y^2/2).
  pair<double, double> gauss2() {double r = sqrt(-2. * log(flat()));
    double phi = 2. * M_PI * flat();
    return { r * sin(phi), r * cos(phi) };}

  // Generate a random number according to a Gamma-distribution.
  double gamma(double k0, double r0);

  // Generate two random vectors according to the phase space distribution
  pair<Vec4, Vec4> phaseSpace2(double eCM, double m1, double m2);

  // Pick one option among  vector of (positive) probabilities.
  int pick(const vector<double>& prob) ;

  // Randomly shuffle a vector, standard Fisher-Yates algorithm.
  template<typename T> void shuffle(vector<T>& vec);

  // Save or read current state to or from a binary file.
  bool dumpState(string fileName);
  bool readState(string fileName);

  // Get or set the state of the random number generator.
  RndmState getState() const {return stateSave;}
  void setState(const RndmState& state) {stateSave = state;}

  // The default seed, i.e. the Marsaglia-Zaman random number sequence.
  static constexpr int DEFAULTSEED = 19780503;

#ifdef RNGDEBUG
  // Random number methods used for debugging only.
  double flatDebug(string loc, string file, int line);
  double xexpDebug(string loc, string file, int line);
  double gaussDebug(string loc, string file, int line);
  pair<double, double> gauss2Debug(string loc, string file, int line);
  double gammaDebug(string loc, string file, int line, double k0, double r0);
  pair<Vec4, Vec4> phaseSpace2Debug(string loc, string file, int line,
    double eCM, double m1, double m2);

  // Static members for debugging to print call file location or filter.
  static bool debugNow, debugLocation, debugIndex;
  static int debugPrecision, debugCall;
  static set<string> debugStarts, debugEnds, debugContains, debugMatches;
#endif

private:

  // State of the random number generator.
  bool      initRndm;
  RndmState stateSave;

  // Pointer for external random number generation.
  bool   useExternalRndm;
  RndmEnginePtr rndmEngPtr{};

};

//==========================================================================

// Hist class.
// This class handles a single histogram at a time.

class Hist {

public:

  // Constructors, including copy constructors.
  Hist() : titleSave(""), nBin(), nFill(), nNonFinite(), xMin(),
    xMax(), linX(), doStats(), dx(), under(), inside(), over(), sumxNw()
    { }
  Hist(string titleIn, int nBinIn = 100, double xMinIn = 0.,
    double xMaxIn = 1., bool logXIn = false, bool doStatsIn = false) :
    nBin(), nFill(), nNonFinite(), xMin(), xMax(), linX(), doStats(), dx(),
      under(), inside(), over(), sumxNw()
  { book(titleIn, nBinIn, xMinIn, xMaxIn, logXIn, doStatsIn); }
  Hist(const Hist& h)
    : titleSave(h.titleSave), nBin(h.nBin), nFill(h.nFill),
      nNonFinite(h.nNonFinite), xMin(h.xMin), xMax(h.xMax), linX(h.linX),
      doStats(h.doStats), dx(h.dx), under(h.under), inside(h.inside),
      over(h.over), res(h.res), res2(h.res2), sumxNw() {
    for (int i = 0; i < nMoments; ++i) sumxNw[i] = h.sumxNw[i];
  }
  Hist(string titleIn, const Hist& h)
    : titleSave(titleIn), nBin(h.nBin), nFill(h.nFill),
      nNonFinite(h.nNonFinite), xMin(h.xMin), xMax(h.xMax), linX(h.linX),
      doStats(h.doStats), dx(h.dx), under(h.under), inside(h.inside),
      over(h.over), res(h.res), res2(h.res2), sumxNw() {
    for (int i = 0; i < nMoments; ++i) sumxNw[i] = h.sumxNw[i];
  }
  Hist& operator=(const Hist& h) { if(this != &h) {
    nBin = h.nBin; nFill = h.nFill; nNonFinite = h.nNonFinite; xMin = h.xMin;
    xMax = h.xMax; linX = h.linX; doStats = h.doStats; dx = h.dx;
    under = h.under; inside = h.inside; over = h.over;
    for (int i = 0; i < nMoments; ++i) sumxNw[i] = h.sumxNw[i];
    res = h.res; res2 = h.res2; } return *this; }

  // Create a histogram that is the plot of the given function.
  static Hist plotFunc(function<double(double)> f, string titleIn,
    int nBinIn, double xMinIn, double xMaxIn, bool logXIn = false);

  // Book a histogram.
  void book(string titleIn = "  ", int nBinIn = 100, double xMinIn = 0.,
    double xMaxIn = 1., bool logXIn = false, bool doStatsIn = false) ;

  // Set title of a histogram.
  void title(string titleIn = "  ") {titleSave = titleIn; }

  // Reset bin contents.
  void null() ;

  // Fill bin with weight.
  void fill(double x, double w = 1.) ;

  // Print a histogram with overloaded << operator.
  friend ostream& operator<<(ostream& os, const Hist& h) ;

  // Print histogram contents as a table (e.g. for Gnuplot, Rivet or Pyplot),
  // optionally with statistical errors.
  void table(ostream& os = cout, bool printOverUnder = false,
    bool xMidBin = true, bool printError = false) const ;
  void table(string fileName, bool printOverUnder = false,
    bool xMidBin = true, bool printError = false) const {
    ofstream streamName(fileName.c_str());
    table(streamName, printOverUnder, xMidBin, printError);}
  void rivetTable(ostream& os = cout, bool printError = true) const ;
  void rivetTable(string fileName, bool printError = true) const {
    ofstream streamName(fileName.c_str()); rivetTable(streamName, printError);}
  void pyplotTable(ostream& os = cout, bool isHist = true,
    bool printError = false) const ;
  void pyplotTable(string fileName, bool isHist = true,
    bool printError = false) const {ofstream streamName(fileName.c_str());
    pyplotTable(streamName, isHist, printError);}

  // Fill contents of a two-column (x,y) table, e.g. written by table() above.
  void fillTable(istream& is = cin);
  void fillTable(string fileName) { ifstream streamName(fileName.c_str());
    fillTable(streamName);}

  // Print a table out of two histograms with same x axis (no errors printed).
  friend void table(const Hist& h1, const Hist& h2, ostream& os,
    bool printOverUnder, bool xMidBin) ;
  friend void table(const Hist& h1, const Hist& h2, string fileName,
    bool printOverUnder, bool xMidBin) ;

  // Return title and size of histogram. Also if logarithmic x scale.
  string getTitle() const {return titleSave;}
  int    getBinNumber() const {return nBin;}
  int    getNonFinite() const {return nNonFinite;}
  bool   getLinX() const {return linX;}

  // Return min and max in x and y directions.
  double getXMin() const {return xMin;}
  double getXMax() const {return xMax;}
  double getYMin() const { if (nBin == 0) return 0.;
    double yMin = res[0];
    for (int ix = 1; ix < nBin; ++ix)
      if (res[ix] < yMin ) yMin = res[ix];
    return yMin;}
  double getYMax() const { if (nBin == 0) return 0.;
    double yMax = res[0];
    for (int ix = 1; ix < nBin; ++ix)
      if (res[ix] > yMax ) yMax = res[ix];
    return yMax;}
  double getYAbsMin() const { double yAbsMin = 1e20; double yAbs;
    for (int ix = 0; ix < nBin; ++ix) { yAbs = abs(res[ix]);
      if (yAbs > 1e-20 && yAbs < yAbsMin) yAbsMin = yAbs; }
    return yAbsMin;}

  // Return <X> and error on <X>, unbinned from saved weight sums (default)
  // or directly from the histogram bins (unbinned = false). In the latter
  // case, the error estimate includes the difference between the binned and
  // unbinned value summed in quadrature with the statistical error, as a
  // measure of bin granularity error.
  double getXMean(bool unbinned=true) const;
  double getXMeanErr(bool unbinned=true) const;

  // Return Median in X and its statistical error, ignoring underflow and
  // overflow (default) or including them (includeOverUnder = true). By
  // default, error includes granularity estimate obtained by comparing binned
  // vs unbinned mean value, but this can be switched off (unbinned = false).
  double getXMedian(bool includeOverUnder=false) const;
  double getXMedianErr(bool unbinned=true) const;

  // Return average <Y> value.
  double getYMean() const { return inside / nFill; }

  // Return RMS and equivalent n'th roots of n'th moments about the mean,
  // and their error estimates. Up to n = 6, both unbinned and binned moments
  // can be calculated. For n >= 7, and for all error estimates, only
  // binned values are available. Note that (unlike ROOT), the error estimates
  // do not assume normal distributions.
  // RMN(2) = RMS is the standard root-mean-square deviation from the mean.
  // RMN(3) is the cube root of the mean-cube deviation,
  //         cbrt(<(x - <x>)^3>). It is sensitive to single-sided tails,
  //         as are characteristic of many particle-physics distributions.
  // RMN(4) adds sensitivity to double-sided long tails (eg BW vs
  //         Gaussian), and further sensitivity to long single-sided ones.
  // Etc.
  double getXRMN(int n=2, bool unbinned=true) const;
  double getXRMS(bool unbinned=true) const {return getXRMN(2, unbinned);}

  double getXRMNErr(int n=2, bool unbinned=true) const;
  double getXRMSErr(bool unbinned=true) const {
    return getXRMNErr(2, unbinned);}

  // Return content of specific bin: 0 gives underflow and nBin+1 overflow.
  double getBinContent(int iBin) const;

  // Return the lower edge of the bin.
  double getBinEdge(int iBin) const;

  // Return the width of the bin.
  double getBinWidth(int iBin=1) const;

  // Return bin contents.
  vector<double> getBinContents() const;

  // Return bin edges.
  vector<double> getBinEdges() const;

  // Return number of entries.
  int getEntries(bool alsoNonFinite = true) const {
    return alsoNonFinite ? nNonFinite + nFill : nFill; }

  // Return sum of weights.
  double getWeightSum(bool alsoOverUnder = true) const {
    return alsoOverUnder ? inside + over + under : inside; }

  // Return effective entries (for weighted histograms = number
  // of equivalent unweighted events for same statistical power).
  double getNEffective() const {
    double sumw2 = 0.;
    for (int ix = 0; ix < nBin; ++ix) sumw2 += res2[ix];
    if (sumw2 <= Hist::TINY) return 0.;
    else return pow2(sumxNw[0]) / sumw2;
  }

  // Check whether another histogram has same size and limits.
  bool sameSize(const Hist& h) const ;

  // Take an arbitrary function of bin contents.
  void takeFunc(function<double(double)> func);

  // Take logarithm (base 10 or e) of bin contents.
  void takeLog(bool tenLog = true);

  // Take square root of bin contents.
  void takeSqrt();

  // Normalize sum of bin contents to value, with or without overflow bins.
  void normalize(double f = 1, bool overflow = true) ;

  // Normalize sum of bin areas to value, with or without overflow bins.
  void normalizeIntegral(double f = 1, bool overflow = true);

  // Scale each bin content by 1 / (wtSum * bin width).
  void normalizeSpectrum(double wtSum);

  // Operator overloading with member functions
  Hist& operator+=(const Hist& h) ;
  Hist& operator-=(const Hist& h) ;
  Hist& operator*=(const Hist& h) ;
  Hist& operator/=(const Hist& h) ;
  Hist& operator+=(double f) ;
  Hist& operator-=(double f) ;
  Hist& operator*=(double f) ;
  Hist& operator/=(double f) ;
  Hist operator+(double f) const;
  Hist operator+(const Hist& h2) const;
  Hist operator-(double f) const;
  Hist operator-(const Hist& h2) const;
  Hist operator*(double f) const;
  Hist operator*(const Hist& h2) const;
  Hist operator/(double f) const;
  Hist operator/(const Hist& h2) const;

  // Operator overloading with friends
  friend Hist operator+(double f, const Hist& h1);
  friend Hist operator-(double f, const Hist& h1);
  friend Hist operator*(double f, const Hist& h1);
  friend Hist operator/(double f, const Hist& h1);

private:

  // Constants: could only be changed in the code itself.
  static const int    NBINMAX, NCOLMAX, NLINES;
  static const double TOLERANCE, TINY, LARGE, SMALLFRAC, DYAC[];
  static const char   NUMBER[];

  // Properties and contents of a histogram.
  string titleSave;
  int    nBin, nFill, nNonFinite;
  double xMin, xMax;
  bool   linX, doStats;
  double dx, under, inside, over;
  vector<double> res, res2;

  // Sum x^N w, for different powers N, for calculation of unbinned moments.
  static constexpr int nMoments = 7;
  double sumxNw[nMoments];

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of Hist class.

// Print a histogram with overloaded << operator.
ostream& operator<<(ostream& os, const Hist& h) ;

// Print a table out of two histograms with same x axis.
void table(const Hist& h1, const Hist& h2, ostream& os = cout,
  bool printOverUnder = false, bool xMidBin = true) ;
void table(const Hist& h1, const Hist& h2, string fileName,
  bool printOverUnder = false, bool xMidBin = true) ;

// Operator overloading with friends
Hist operator+(double f, const Hist& h1);
Hist operator-(double f, const Hist& h1);
Hist operator*(double f, const Hist& h1);
Hist operator/(double f, const Hist& h1);

//==========================================================================

// HistPlot class.
// Writes a Python program that can generate PDF plots from Hist histograms.

class HistPlot {

public:

  // Constructor requires name of Python program (and adds .py).
  HistPlot(string pythonName, bool useLegacyIn = false)
    : nFrame(), nTable(), useLegacy(useLegacyIn) {
    toPython.open( (pythonName + ".py").c_str() );
    toPython << "from matplotlib import pyplot as plt" << endl
             << "from matplotlib.backends.backend_pdf import PdfPages" << endl;
    nPDF = 0; }

  // Destructor should do final close.
  ~HistPlot() { toPython << "pp.close()" << endl; }

  // New plot frame, with title, x and y labels, x and y sizes..
  void frame( string frameIn, string titleIn = "", string xLabIn = "",
    string yLabIn = "", double xSizeIn = 8., double ySizeIn = 6.) {
    framePrevious = frameName; frameName = frameIn; title = titleIn;
    xLabel = xLabIn; yLabel = yLabIn; xSize = xSizeIn; ySize = ySizeIn;
    histos.clear(); styles.clear(); legends.clear(); files.clear();
    fileStyles.clear(); fileLegends.clear(); filexyerr.clear();}

  // Add a histogram to the current plot, with optional style and legend.
  void add( const Hist& histIn, string styleIn = "h",
    string legendIn = "void") {
    if (histIn.getBinNumber() == 0) {
      cout << " Error: histogram is not booked" << endl;
      return;
    }
    histos.push_back(histIn);
    styles.push_back(styleIn); legends.push_back(legendIn); }

  // Add a file of (x, y) values not from a histogram, e.g. data points.
  void addFile( string fileIn, string styleIn = "o",
    string legendIn = "void", string xyerrIn="") { files.push_back(fileIn);
    fileStyles.push_back(styleIn); fileLegends.push_back(legendIn);
    filexyerr.push_back(xyerrIn);}

  // Plot a frame given the information from the new and add calls.
  void plot( bool logY = false, bool logX = false, bool userBorders = false);
  void plot( double xMinUserIn, double xMaxUserIn,  double yMinUserIn,
     double yMaxUserIn, bool logY = false, bool logX = false) {
     xMinUser = xMinUserIn; xMaxUser = xMaxUserIn; yMinUser = yMinUserIn;
     yMaxUser = yMaxUserIn; plot( logY, logX, true);}

  //  Omnibus single call when only one histogram in the frame.
  void plotFrame( string frameIn, const Hist& histIn, string titleIn = "",
    string xLabIn = "", string yLabIn = "", string styleIn = "h",
    string legendIn = "void",  bool logY = false) {
    frame( frameIn, titleIn, xLabIn, yLabIn);
    add( histIn, styleIn, legendIn); plot( logY); }

private:

  // Initialization code.
  void init( string pythonName);

  // Stored quantities.
  ofstream toPython;
  int      nPDF, nFrame, nTable;
  double   xSize, ySize, xMinUser, xMaxUser, yMinUser, yMaxUser;
  string   frameName, framePrevious, title, xLabel, yLabel, fileName, tmpFig;
  vector<Hist> histos;
  vector<string> styles, legends, files, fileStyles, fileLegends, filexyerr;

  // If true, use old linthreshy matplotlib parameter (removed in 3.5.0)
  bool useLegacy;

};

//==========================================================================

// Methods used for debugging random number sequences.

#ifdef RNGDEBUG
#define flat() flatDebug(__METHOD_NAME__, __FILE__, __LINE__)
#define xexp() xexpDebug(__METHOD_NAME__, __FILE__, __LINE__)
#define gauss() gaussDebug(__METHOD_NAME__, __FILE__, __LINE__)
#define gamma(...) gammaDebug(__METHOD_NAME__, __FILE__, __LINE__, __VA_ARGS__)
#define phaseSpace2(...) phaseSpace2Debug(__METHOD_NAME__, __FILE__, __LINE__,\
    __VA_ARGS__)
#endif

//==========================================================================

// Randomly shuffle a vector, standard Fisher-Yates algorithm.
// This must be defined after possible RNG debugging.

template<typename T> void Rndm::shuffle(vector<T>& vec) {
  for (int i = vec.size() - 1; i > 0; --i)
    swap(vec[i], vec[floor(flat() * (i + 1))]);
}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Basics_H
