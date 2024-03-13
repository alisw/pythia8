// Basics.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Rndm, Vec4,
// RotBstMatrix and Hist classes, and some related global functions.

#include "Pythia8/Basics.h"

// Access time information.
#include <ctime>
#include <limits>

namespace Pythia8 {

//==========================================================================

// Rndm class.
// This class handles random number generation according to the
// Marsaglia-Zaman-Tsang algorithm

//--------------------------------------------------------------------------

// Method to pass in pointer for external random number generation.

bool Rndm::rndmEnginePtr( RndmEnginePtr rndmEngPtrIn) {

  // Save pointer.
  if (rndmEngPtrIn == nullptr) return false;
  rndmEngPtr      = rndmEngPtrIn;
  useExternalRndm = true;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Initialize, normally at construction or in first call.

void Rndm::init(int seedIn) {

  // Pick seed in convenient way. Assure it to be non-negative.
  int seed = seedIn;
  if (seedIn < 0) seed = DEFAULTSEED;
  else if (seedIn == 0) seed = int(time(0));
  if (seed < 0) seed = -seed;

  // Unpack seed.
  int ij = (seed/30082) % 31329;
  int kl = seed % 30082;
  int i  = (ij/177) % 177 + 2;
  int j  = ij % 177 + 2;
  int k  = (kl/169) % 178 + 1;
  int l  =  kl % 169;

  // Initialize random number array.
  for (int ii = 0; ii < 97; ++ii) {
    double s = 0.;
    double t = 0.5;
    for (int jj = 0; jj < 48; ++jj) {
      int m = (( (i*j)%179 )*k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ( (l*m) % 64 >= 32) s += t;
      t *= 0.5;
    }
    stateSave.u[ii] = s;
  }

  // Initialize other variables.
  double twom24 = 1.;
  for (int i24 = 0; i24 < 24; ++i24) twom24 *= 0.5;
  stateSave.c   = 362436. * twom24;
  stateSave.cd  = 7654321. * twom24;
  stateSave.cm  = 16777213. * twom24;
  stateSave.i97 = 96;
  stateSave.j97 = 32;

  // Finished.
  initRndm = true;
  stateSave.seed = seed;
  stateSave.sequence = 0;

}

//--------------------------------------------------------------------------

// Generate random numbers according to exp(-x).
// Must be defined before possible RNG debugging methods.

double Rndm::exp() { return -log(flat()) ;}

//--------------------------------------------------------------------------

// Pick one option among vector of (positive) probabilities.
// Must be defined before possible RNG debugging methods.

int Rndm::pick(const vector<double>& prob) {

  double work = 0.;
  for (int i = 0; i < int(prob.size()); ++i) work += prob[i];
  work *= flat();
  int index = -1;
  do work -= prob[++index];
  while (work > 0. && index < int(prob.size()));
  return index;

}

//--------------------------------------------------------------------------

// Use standard random number generation, rather than debug versions.

#include "Pythia8/RngDebug.h"

//--------------------------------------------------------------------------

// Define debug random number calls.

#ifdef RNGDEBUG
// Flags to control debugging behaviour.
bool        Rndm::debugNow       = true;
bool        Rndm::debugLocation  = true;
bool        Rndm::debugIndex     = false;
int         Rndm::debugPrecision = 4;
int         Rndm::debugCall      = 0;
set<string> Rndm::debugStarts    = {};
set<string> Rndm::debugEnds      = {};
set<string> Rndm::debugContains  = {};
set<string> Rndm::debugMatches   = {};
// Debug random number calls.
double rngDebug(double val, string loc, string call, string file, int line) {
  ++Rndm::debugCall;
  if (!Rndm::debugNow) return val;
  bool print = Rndm::debugStarts.size() + Rndm::debugEnds.size() +
    Rndm::debugContains.size() + Rndm::debugMatches.size() == 0;
  for (auto &exp : Rndm::debugStarts)
    print = print || loc.rfind(exp, 0) == 0;
  for (auto &exp : Rndm::debugEnds)
    print = print || (loc.length() >= exp.length()
      && loc.compare(loc.length() - exp.length(), exp.length(), exp) == 0);
  for (auto &exp : Rndm::debugContains)
    print = print || loc.find(exp) != string::npos;
  for (auto &exp : Rndm::debugMatches) print = print || loc == exp;
  if (!print) return val;
  cout << setw(80) << left << loc + ":" + call
       << setprecision(Rndm::debugPrecision) << scientific
       << setw(Rndm::debugPrecision + 8) << right << val;
  if (Rndm::debugIndex) cout << setw(12) << right << Rndm::debugCall;
  if (Rndm::debugLocation) cout << " " << file << ":" << line;
  cout << "\n";
  return val;
}
double Rndm::flatDebug(string loc, string file, int line) {
  return rngDebug(flat(), loc, "flat", file, line);}
double Rndm::xexpDebug(string loc, string file, int line) {
  return rngDebug(xexp(), loc, "xexp", file, line);}
double Rndm::gaussDebug(string loc, string file, int line) {
  return rngDebug(gauss(), loc, "gauss", file, line);}
pair<double, double> Rndm::gauss2Debug(string loc, string file, int line) {
  pair<double, double> val = gauss2();
  rngDebug(val.first, loc, "gauss2:first", file, line);
  rngDebug(val.second, loc, "gauss2:second", file, line);
  return val;}
double Rndm::gammaDebug(string loc, string file, int line,
  double k0, double r0) {
  return rngDebug(gamma(k0, r0), loc, "gamma", file, line);}
pair<Vec4, Vec4> Rndm::phaseSpace2Debug(string loc, string file, int line,
  double eCM, double m1, double m2) {
  pair<Vec4, Vec4> val = phaseSpace2(eCM, m1, m2);
  for (int i = 0; i < 4; ++i) {
    rngDebug(val.first[i], loc, "phaseSpace2:first:" + to_string(i),
      file, line);
    rngDebug(val.second[i], loc, "phaseSpace2:second:" + to_string(i),
      file, line);
  }
  return val;}
#endif

//--------------------------------------------------------------------------

// Generate next random number uniformly between 0 and 1.

double Rndm::flat() {

  // Use external random number generator if such has been linked.
  if (useExternalRndm) return rndmEngPtr->flat();

  // Ensure that already initialized.
  if (!initRndm) init(DEFAULTSEED);

  // Find next random number and update saved state.
  ++stateSave.sequence;
  double uni;
  do {
    uni = stateSave.u[stateSave.i97] - stateSave.u[stateSave.j97];
    if (uni < 0.) uni += 1.;
    stateSave.u[stateSave.i97] = uni;
    if (--stateSave.i97 < 0) stateSave.i97 = 96;
    if (--stateSave.j97 < 0) stateSave.j97 = 96;
    stateSave.c -= stateSave.cd;
    if (stateSave.c < 0.) stateSave.c += stateSave.cm;
    uni -= stateSave.c;
    if(uni < 0.) uni += 1.;
   } while (uni <= 0. || uni >= 1.);
  return uni;

}

//--------------------------------------------------------------------------

// Generate a random number according to a Gamma-distribution.

double Rndm::gamma(double k0, double r0) {
  int k = int(k0);

  double x = 0.0;
  for ( int i = 0; i < k; ++i ) x += -log(flat());

  double del = k0 - k;
  if ( del == 0.0 ) return x*r0;

  while ( true ) {
    double U = flat();
    double V = flat();
    double W = flat();

    double xi = 0.0;
    if ( U <= M_E / (M_E + del) ) {
      xi = pow(V, 1.0/del);
      if ( W <= std::exp(-xi) ) return r0*(x + xi);
    } else {
      xi = 1.0 - log(V);
      if ( W <= pow(xi, del - 1.0) ) return r0*(x + xi);
    }
  }
  return 0.0;
}

//--------------------------------------------------------------------------

// Generate two random vectors according to the phase space distribution

pair<Vec4, Vec4> Rndm::phaseSpace2(double eCM, double m1, double m2) {

  // Calculate phase space configuration.
  double pAbs = 0.5 * sqrtpos( (eCM - m1 - m2) * (eCM + m1 + m2)
    * (eCM + m1 - m2) * (eCM - m1 + m2) ) / eCM;

  // Isotropic angles give three-momentum.
  double cosTheta = 2. * flat() - 1.;
  double sinTheta = sqrt(1. - cosTheta*cosTheta);
  double phi      = 2. * M_PI * flat();
  double pX       = pAbs * sinTheta * cos(phi);
  double pY       = pAbs * sinTheta * sin(phi);
  double pZ       = pAbs * cosTheta;
  double e1   = sqrt(pAbs * pAbs + m1 * m1);
  double e2   = sqrt(pAbs * pAbs + m2 * m2);

  return { Vec4(pX, pY, pZ, e1), Vec4(-pX, -pY, -pZ, e2) };
}

//--------------------------------------------------------------------------

// Save current state of the random number generator to a binary file.

bool Rndm::dumpState(string fileName) {

  // Open file as output stream.
  const char* fn = fileName.c_str();
  ofstream ofs(fn, ios::binary);

  if (!ofs.good()) {
    cout << " Rndm::dumpState: could not open output file" << endl;
    return false;
  }

  // Write the state of the generator on the file.
  ofs.write((char *) &stateSave.seed,     sizeof(int));
  ofs.write((char *) &stateSave.sequence, sizeof(long));
  ofs.write((char *) &stateSave.i97,      sizeof(int));
  ofs.write((char *) &stateSave.j97,      sizeof(int));
  ofs.write((char *) &stateSave.c,        sizeof(double));
  ofs.write((char *) &stateSave.cd,       sizeof(double));
  ofs.write((char *) &stateSave.cm,       sizeof(double));
  ofs.write((char *) &stateSave.u,        sizeof(double) * 97);

  // Write confirmation on cout.
  cout << " PYTHIA Rndm::dumpState: seed = " << stateSave.seed
       << ", sequence no = " << stateSave.sequence << endl;
  return true;

}

//--------------------------------------------------------------------------

// Read in the state of the random number generator from a binary file.

bool Rndm::readState(string fileName) {

  // Open file as input stream.
  const char* fn = fileName.c_str();
  ifstream ifs(fn, ios::binary);

  if (!ifs.good()) {
    cout << " Rndm::readState: could not open input file" << endl;
    return false;
  }

  // Read the state of the generator from the file.
  ifs.read((char *) &stateSave.seed,     sizeof(int));
  ifs.read((char *) &stateSave.sequence, sizeof(long));
  ifs.read((char *) &stateSave.i97,      sizeof(int));
  ifs.read((char *) &stateSave.j97,      sizeof(int));
  ifs.read((char *) &stateSave.c,        sizeof(double));
  ifs.read((char *) &stateSave.cd,       sizeof(double));
  ifs.read((char *) &stateSave.cm,       sizeof(double));
  ifs.read((char *) &stateSave.u,        sizeof(double) *97);

  // Write confirmation on cout.
  cout << " PYTHIA Rndm::readState: seed " << stateSave.seed
       << ", sequence no = " << stateSave.sequence << endl;
  return true;

}

//==========================================================================

// RndmState class.
// This class describes the configuration of a Rndm object.

//--------------------------------------------------------------------------

// Test whether two random states would generate the same random sequence.
bool RndmState::operator==(const RndmState& other) const {
  if (i97 != other.i97 || j97 != other.j97 || sequence != other.sequence
   || c != other.c || cd != other.cd || cm != other.cm)
    return false;
  for (int i = 0; i < 97; ++i)
    if (u[i] != other.u[i])
      return false;
  return true;
}

//==========================================================================

// Vec4 class.
// This class implements four-vectors, in energy-momentum space.
// (But could also be used to hold space-time four-vectors.)

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Small number to avoid division by zero.
const double Vec4::TINY = 1e-20;

//--------------------------------------------------------------------------

// Rotation (simple).

void Vec4::rot(double thetaIn, double phiIn) {

  double cthe = cos(thetaIn);
  double sthe = sin(thetaIn);
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double tmpx =  cthe * cphi * xx -    sphi * yy + sthe * cphi * zz;
  double tmpy =  cthe * sphi * xx +    cphi * yy + sthe * sphi * zz;
  double tmpz = -sthe *        xx +                cthe *        zz;
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Azimuthal rotation phi around an arbitrary axis (nz, ny, nz).

void Vec4::rotaxis(double phiIn, double nx, double ny, double nz) {

  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx         *= norm;
  ny         *= norm;
  nz         *= norm;
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Azimuthal rotation phi around an arbitrary (3-vector component of) axis.

void Vec4::rotaxis(double phiIn, const Vec4& n) {

  double nx   = n.xx;
  double ny   = n.yy;
  double nz   = n.zz;
  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx         *= norm;
  ny          *=norm;
  nz          *=norm;
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Boost (simple).

void Vec4::bst(double betaX, double betaY, double betaZ) {

  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  if (beta2 >= 1.) return;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx += prod2 * betaX;
  yy += prod2 * betaY;
  zz += prod2 * betaZ;
  tt = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost (simple, given gamma).

void Vec4::bst(double betaX, double betaY, double betaZ, double gamma) {

  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx += prod2 * betaX;
  yy += prod2 * betaY;
  zz += prod2 * betaZ;
  tt = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p.

void Vec4::bst(const Vec4& pIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = pIn.xx / pIn.tt;
  double betaY = pIn.yy / pIn.tt;
  double betaZ = pIn.zz / pIn.tt;
  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  if (beta2 >= 1.) return;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p and double m.

void Vec4::bst(const Vec4& pIn, double mIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = pIn.xx / pIn.tt;
  double betaY = pIn.yy / pIn.tt;
  double betaZ = pIn.zz / pIn.tt;
  double gamma = pIn.tt / mIn;
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p; boost in opposite direction.

void Vec4::bstback(const Vec4& pIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = -pIn.xx / pIn.tt;
  double betaY = -pIn.yy / pIn.tt;
  double betaZ = -pIn.zz / pIn.tt;
  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  if (beta2 >= 1.) return;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p and double m; boost in opposite direction.

void Vec4::bstback(const Vec4& pIn, double mIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = -pIn.xx / pIn.tt;
  double betaY = -pIn.yy / pIn.tt;
  double betaZ = -pIn.zz / pIn.tt;
  double gamma = pIn.tt / mIn;
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Arbitrary combination of rotations and boosts defined by 4 * 4 matrix.

void Vec4::rotbst(const RotBstMatrix& M) {

  double x = xx; double y = yy; double z = zz; double t = tt;
  tt = M.M[0][0] * t + M.M[0][1] * x + M.M[0][2] * y +  M.M[0][3] * z;
  xx = M.M[1][0] * t + M.M[1][1] * x + M.M[1][2] * y +  M.M[1][3] * z;
  yy = M.M[2][0] * t + M.M[2][1] * x + M.M[2][2] * y +  M.M[2][3] * z;
  zz = M.M[3][0] * t + M.M[3][1] * x + M.M[3][2] * y +  M.M[3][3] * z;

}

//--------------------------------------------------------------------------

// Print a four-vector: also operator overloading with friend.

ostream& operator<<(ostream& os, const Vec4& v) {
  os << fixed << setprecision(3) << " " << setw(9) << v.xx << " "
     << setw(9) << v.yy << " " << setw(9) << v.zz << " " << setw(9)
     << v.tt << " (" << setw(9) << v.mCalc() << ")\n";
  return os;
}

//--------------------------------------------------------------------------

// The invariant mass of one or more four-vectors.

double m(const Vec4& v1) {
  double s = m2(v1); return (s >= 0) ? sqrt(s) : -sqrt(-s);}

double m(const Vec4& v1, const Vec4& v2) {
  double s = m2(v1, v2); return (s > 0.) ? sqrt(s) : 0.;}

//--------------------------------------------------------------------------

// The squared invariant mass of one or more four-vectors.

double m2(const Vec4& v1) {
  return pow2(v1.tt) - pow2(v1.xx) - pow2(v1.yy) - pow2(v1.zz);}

double m2(const Vec4& v1, const Vec4& v2) {
  return pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
     - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);}

double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3) {
  return pow2(v1.tt + v2.tt + v3.tt) - pow2(v1.xx + v2.xx + v3.xx)
    - pow2(v1.yy + v2.yy + v3.yy) - pow2(v1.zz + v2.zz + v3.zz);}

double m2(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Vec4& v4) {
  return pow2(v1.tt + v2.tt + v3.tt + v4.tt)
    - pow2(v1.xx + v2.xx + v3.xx + v4.xx) - pow2(v1.yy + v2.yy + v3.yy + v4.yy)
    - pow2(v1.zz + v2.zz + v3.zz + v4.zz);}

//--------------------------------------------------------------------------

// The scalar product of two three-vectors.

double dot3(const Vec4& v1, const Vec4& v2) {
  return v1.xx*v2.xx + v1.yy*v2.yy + v1.zz*v2.zz;
}

//--------------------------------------------------------------------------

// The cross product of two three-vectors.

Vec4 cross3(const Vec4& v1, const Vec4& v2) {
  Vec4 v;
  v.xx = v1.yy * v2.zz - v1.zz * v2.yy;
  v.yy = v1.zz * v2.xx - v1.xx * v2.zz;
  v.zz = v1.xx * v2.yy - v1.yy * v2.xx; return v;
}


//--------------------------------------------------------------------------

// Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c)

Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c) {
  Vec4 v(0.,0.,0.,0.);
  v.tt =   a.xx*b.yy*c.zz + a.yy*b.zz*c.xx + a.zz*b.xx*c.yy
         - a.xx*b.zz*c.yy - a.zz*b.yy*c.xx - a.yy*b.xx*c.zz;
  v.xx = -(- a.tt*b.yy*c.zz - a.yy*b.zz*c.tt - a.zz*b.tt*c.yy
           + a.tt*b.zz*c.yy + a.zz*b.yy*c.tt + a.yy*b.tt*c.zz);
  v.yy = -(- a.xx*b.tt*c.zz - a.tt*b.zz*c.xx - a.zz*b.xx*c.tt
           + a.xx*b.zz*c.tt + a.zz*b.tt*c.xx + a.tt*b.xx*c.zz);
  v.zz = -(- a.xx*b.yy*c.tt - a.yy*b.tt*c.xx - a.tt*b.xx*c.yy
           + a.xx*b.tt*c.yy + a.tt*b.yy*c.xx + a.yy*b.xx*c.tt);
  return v;
}

//--------------------------------------------------------------------------

// Opening angle between two three-vectors.

double theta(const Vec4& v1, const Vec4& v2) {
  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
    / sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
    * (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
  cthe = max(-1., min(1., cthe));
  return acos(cthe);
}

//--------------------------------------------------------------------------

// Cosine of the opening angle between two three-vectors.

double costheta(const Vec4& v1, const Vec4& v2) {
  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
    / sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
    * (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
  cthe = max(-1., min(1., cthe));
  return cthe;
}

//--------------------------------------------------------------------------

// Cosine of the opening angle between two particles.

double costheta(double e1, double e2, double m1, double m2, double s12) {
  return (2.0*e1*e2 - s12)/(2.0*sqrt(e1*e1 - m1*m1)*sqrt(e2*e2 - m2*m2));}

//--------------------------------------------------------------------------

// Azimuthal angle between two three-vectors.

double phi(const Vec4& v1, const Vec4& v2) {
  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( Vec4::TINY,
    (v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
  cphi = max(-1., min(1., cphi));
  return acos(cphi);
}

//--------------------------------------------------------------------------

// Cosine of the azimuthal angle between two three-vectors.

double cosphi(const Vec4& v1, const Vec4& v2) {
  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( Vec4::TINY,
    (v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
  cphi = max(-1., min(1., cphi));
  return cphi;
}

//--------------------------------------------------------------------------

// Azimuthal angle between two three-vectors around a third.

double phi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
  double nx = n.xx; double ny = n.yy; double nz = n.zz;
  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
  nx *= norm; ny *=norm; nz *=norm;
  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
  double cphi = (v1v2 - v1n * v2n) / sqrt( max( Vec4::TINY,
    (v1s - v1n*v1n) * (v2s - v2n*v2n) ));
  cphi = max(-1., min(1., cphi));
  return acos(cphi);
}

//--------------------------------------------------------------------------

// Cosine of the azimuthal angle between two three-vectors around a third.

double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
  double nx = n.xx; double ny = n.yy; double nz = n.zz;
  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
  nx *= norm; ny *=norm; nz *=norm;
  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
  double cphi = (v1v2 - v1n * v2n) / sqrt( max( Vec4::TINY,
    (v1s - v1n*v1n) * (v2s - v2n*v2n) ));
  cphi = max(-1., min(1., cphi));
  return cphi;
}

//--------------------------------------------------------------------------

// Distance in cylindrical (y, phi) coordinates.

double RRapPhi(const Vec4& v1, const Vec4& v2) {
  double dRap = abs(v1.rap() - v2.rap());
  double dPhi = abs(v1.phi() - v2.phi());
  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
  return sqrt(dRap*dRap + dPhi*dPhi);
}

//--------------------------------------------------------------------------

// Distance in cylindrical (eta, phi) coordinates.

double REtaPhi(const Vec4& v1, const Vec4& v2) {
  double dEta = abs(v1.eta() - v2.eta());
  double dPhi = abs(v1.phi() - v2.phi());
  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
  return sqrt(dEta*dEta + dPhi*dPhi);
}

//--------------------------------------------------------------------------

// Shift four-momenta within pair from old to new masses.
// Note that p1Move and p2Move change values during operation.

bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New) {

  // Standard kinematics variables.
  double sH  = (p1Move + p2Move).m2Calc();
  double r1  = p1Move.m2Calc() / sH;
  double r2  = p2Move.m2Calc() / sH;
  double r3  = m1New * m1New / sH;
  double r4  = m2New * m2New / sH;
  double l12 = sqrtpos(pow2(1. - r1 - r2) - 4. * r1 * r2);
  double l34 = sqrtpos(pow2(1. - r3 - r4) - 4. * r3 * r4);

  // Check that shift operation possible.
  if (sH <= pow2(m1New + m2New) || l12 < Vec4::TINY || l34 < Vec4::TINY)
    return false;

  // Calculate needed shift and apply it.
  double c1  = 0.5 * ( (1. - r1 + r2) * l34 / l12 - (1. - r3 + r4) );
  double c2  = 0.5 * ( (1. + r1 - r2) * l34 / l12 - (1. + r3 - r4) );
  Vec4   pSh = c1 * p1Move - c2 * p2Move;
  p1Move    += pSh;
  p2Move    -= pSh;
  return true;
}

//--------------------------------------------------------------------------

// Create two vectors that are perpendicular to both input vectors.

pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2) {

  // One perpendicular vector from three-dimensional cross-product.
  Vec4 nPerp( cross3(v1,v2) );
  double TINY = numeric_limits<double>::epsilon();
  if ( abs(nPerp.pAbs()) < TINY) {
    Vec4 aux;
    if (v1.px() != 0.)      aux.p(v1.yy,v1.px(),v1.pz(),v1.e());
    else if (v1.py() != 0.) aux.p(v1.px(),v1.pz(),v1.py(),v1.e());
    else if (v1.pz() != 0.) aux.p(v1.pz(),v1.py(),v1.px(),v1.e());
    nPerp.p( cross3(v1,aux) );
  }
  nPerp /= abs(nPerp.pAbs());

  // Second perpendicular vector from four-dimensional cross-product.
  Vec4 lPerp( cross4(v1,v2,nPerp) );
  lPerp /= sqrt(abs(lPerp.m2Calc()));
  return make_pair(nPerp,lPerp);
}

//==========================================================================

// RotBstMatrix class.
// This class implements 4 * 4 matrices that encode an arbitrary combination
// of rotations and boosts, that can be applied to Vec4 four-vectors.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Small number to avoid division by zero.
const double RotBstMatrix::TINY = 1e-20;

//--------------------------------------------------------------------------

// Rotate by polar angle theta and azimuthal angle phi.

void RotBstMatrix::rot(double theta, double phi) {

  // Set up rotation matrix.
  double cthe = cos(theta);
  double sthe = sin(theta);
  double cphi = cos(phi);
  double sphi = sin(phi);
  double Mrot[4][4] = {
    {1.,           0.,         0.,          0.},
    {0.,  cthe * cphi,     - sphi, sthe * cphi},
    {0.,  cthe * sphi,       cphi, sthe * sphi},
    {0., -sthe,                0., cthe       } };

  // Rotate current matrix accordingly.
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = Mrot[i][0] * Mtmp[0][j] + Mrot[i][1] * Mtmp[1][j]
            + Mrot[i][2] * Mtmp[2][j] + Mrot[i][3] * Mtmp[3][j];

}

//--------------------------------------------------------------------------

// Rotate so that vector originally along z axis becomes parallel with p.

void RotBstMatrix::rot(const Vec4& p) {

  double theta = p.theta();
  double phi = p.phi();
  rot(0., -phi);
  rot(theta, phi);

}

//--------------------------------------------------------------------------

// Boost with velocity vector (betaX, betaY, betaZ).

void RotBstMatrix::bst(double betaX, double betaY, double betaZ) {

  // Set up boost matrix.
  double gm = 1. / sqrt( max( TINY, 1. - betaX*betaX - betaY*betaY
    - betaZ*betaZ ) );
  double gf = gm*gm / (1. + gm);
  double Mbst[4][4] = {
    { gm,           gm*betaX,           gm*betaY,          gm*betaZ },
    { gm*betaX, 1. + gf*betaX*betaX, gf*betaX*betaY, gf*betaX*betaZ },
    { gm*betaY, gf*betaY*betaX, 1. + gf*betaY*betaY, gf*betaY*betaZ },
    { gm*betaZ, gf*betaZ*betaX, gf*betaZ*betaY, 1. + gf*betaZ*betaZ } };

  // Boost current matrix correspondingly.
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = Mbst[i][0] * Mtmp[0][j] + Mbst[i][1] * Mtmp[1][j]
            + Mbst[i][2] * Mtmp[2][j] + Mbst[i][3] * Mtmp[3][j];

}

//--------------------------------------------------------------------------

// Boost so that vector originally at rest obtains same velocity as p.

void RotBstMatrix::bst(const Vec4& p) {
  double betaX = p.px() / p.e();
  double betaY = p.py() / p.e();
  double betaZ = p.pz() / p.e();
  bst(betaX, betaY, betaZ);
}

//--------------------------------------------------------------------------

// Boost so vector originally with same velocity as p is brought to rest.

void RotBstMatrix::bstback(const Vec4& p) {
  double betaX = -p.px() / p.e();
  double betaY = -p.py() / p.e();
  double betaZ = -p.pz() / p.e();
  bst(betaX, betaY, betaZ);
}

//--------------------------------------------------------------------------

// Boost that transforms p1 to p2, where p1^2 = p2^2 is assumed.

void RotBstMatrix::bst(const Vec4& p1, const Vec4& p2) {
  double eSum = p1.e() + p2.e();
  double betaX = (p2.px() - p1.px()) / eSum;
  double betaY = (p2.py() - p1.py()) / eSum;
  double betaZ = (p2.pz() - p1.pz()) / eSum;
  double fac = 2. / (1. + betaX*betaX + betaY*betaY + betaZ*betaZ);
  betaX *= fac; betaY *= fac; betaZ *= fac;
  bst(betaX, betaY, betaZ);
}

//--------------------------------------------------------------------------

// Boost and rotation that transforms from p1 and p2
// to their rest frame with p1 along +z axis.

void RotBstMatrix::toCMframe(const Vec4& p1, const Vec4& p2) {
  Vec4 pSum = p1 + p2;
  Vec4 dir  = p1;
  dir.bstback(pSum);
  double theta = dir.theta();
  double phi   = dir.phi();
  bstback(pSum);
  rot(0., -phi);
  rot(-theta, phi);
}

//--------------------------------------------------------------------------

// Rotation and boost that transforms from rest frame of p1 and p2
// with p1 by default along +z axis to actual frame of p1 and p2. (Inverse of
// above.) The option flip handles the case when p1 is along the -z axis in
// the rest frame. This is accomplished by performing the rotation for
// p2 and negating the rotational part.

void RotBstMatrix::fromCMframe(const Vec4& p1, const Vec4& p2, bool flip) {
  Vec4 pSum = p1 + p2;
  Vec4 dir  = (flip) ? p2 : p1;
  dir.bstback(pSum);
  double theta = dir.theta();
  double phi   = dir.phi();
  rot(0., -phi);
  rot(theta, phi);
  if (flip)
    for (int i = 1; i <= 3; ++i)
      for (int j = 1; j <= 3; ++j) M[i][j] *= -1;
  bst(pSum);
}

//--------------------------------------------------------------------------

// Boost and rotation that transforms from p1 and p2 to the frame where
// they have  opposite same-magnitude velocities with p1 along +z axis.

void RotBstMatrix::toSameVframe(const Vec4& p1, const Vec4& p2) {

  // Boost and rotate (p1, p2) = (dir, inv) to CM frame along +-z axis.
  Vec4 pSum = p1 + p2;
  Vec4 dir  = p1;
  Vec4 inv  = p2;
  dir.bstback(pSum);
  inv.bstback(pSum);
  double theta = dir.theta();
  double phi   = dir.phi();
  bstback(pSum);
  rot(0., -phi);
  rot(-theta, phi);

  // Final boost to frame with equal velocities oppositely directed.
  double sDir = p1.m2Calc();
  double sInv = p2.m2Calc();
  if ( abs(sDir - sInv) > 1e-6 * (sDir + sInv) ) {
    double beta = (dir.e() * inv.e() - dir.pAbs2() - sqrt(sDir * sInv))
      * (dir.e() + inv.e()) / (dir.pAbs() * (sDir - sInv));
    bst( 0., 0., beta);
  }
}

//--------------------------------------------------------------------------

// Rotation and boost that transforms from equal-velocity frame of p1 and p2
// with p1 along +z axis to actual frame of p1 and p2. (Inverse of above.)

void RotBstMatrix::fromSameVframe(const Vec4& p1, const Vec4& p2) {

  // Boost and rotation to CM frame along +-z axis.
  Vec4 pSum = p1 + p2;
  Vec4 dir  = p1;
  Vec4 inv  = p2;
  dir.bstback(pSum);
  inv.bstback(pSum);
  double theta = dir.theta();
  double phi   = dir.phi();

  // Initial boost from equal-velocity to equal-momentum frame.
  double sDir = p1.m2Calc();
  double sInv = p2.m2Calc();
  if ( abs(sDir - sInv) > 1e-6 * (sDir + sInv) ) {
    double beta = (dir.e() * inv.e() - dir.pAbs2() - sqrt(sDir * sInv))
      * (dir.e() + inv.e()) / (dir.pAbs() * (sDir - sInv));
    bst( 0., 0., -beta);
  }

  // Rotation and boost back to lab frame.
  rot(0., -phi);
  rot(theta, phi);
  bst(pSum);
}

//--------------------------------------------------------------------------

// Combine existing rotation/boost matrix with another one.

void RotBstMatrix::rotbst(const RotBstMatrix& Mrb) {
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = Mrb.M[i][0] * Mtmp[0][j] + Mrb.M[i][1] * Mtmp[1][j]
            + Mrb.M[i][2] * Mtmp[2][j] + Mrb.M[i][3] * Mtmp[3][j];
}

//--------------------------------------------------------------------------

// Invert the rotation and boost.

void RotBstMatrix::invert() {
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = ( (i == 0 && j > 0) || (i > 0 && j == 0) )
      ? - Mtmp[j][i] : Mtmp[j][i];
}

//--------------------------------------------------------------------------

// Reset to diagonal matrix.

void RotBstMatrix::reset() {
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = (i==j) ? 1. : 0.;
}

//--------------------------------------------------------------------------

// Crude estimate deviation from unit matrix.

double RotBstMatrix::deviation() const {
  double devSum = 0.;
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    devSum += (i==j) ? abs(M[i][j] - 1.) : abs(M[i][j]);
  return devSum;
}

//--------------------------------------------------------------------------

// Print a rotation and boost matrix: operator overloading with friend.

ostream& operator<<(ostream& os, const RotBstMatrix& M) {
  os << fixed << setprecision(5) << "    Rotation/boost matrix: \n";
  for (int i = 0; i <4; ++i)
    os << setw(10) << M.M[i][0] << setw(10) << M.M[i][1]
       << setw(10) << M.M[i][2] << setw(10) << M.M[i][3] << "\n";
  return os;
}

//==========================================================================

// Hist class.
// This class handles a single histogram at a time
// (or a vector of histograms).

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of bins in a histogram.
const int    Hist::NBINMAX   = 10000;

// Maximum number of columns that can be printed for a histogram.
const int    Hist::NCOLMAX   = 100;

// Maximum number of lines a histogram can use at output.
const int    Hist::NLINES    = 30;

// Tolerance in deviation of xMin and xMax between two histograms.
const double Hist::TOLERANCE = 0.001;

// Small and large numbers to avoid division by zero and overflow.
const double Hist::TINY      = 1e-20;
const double Hist::LARGE     = 1e20;

// When minbin/maxbin < SMALLFRAC the y scale goes down to zero.
const double Hist::SMALLFRAC = 0.1;

// Constants for printout: fixed steps on y scale; filling characters.
const double DYAC[] = {0.04, 0.05, 0.06, 0.08, 0.10,
  0.12, 0.15, 0.20, 0.25, 0.30};
const char NUMBER[] = {'0', '1', '2', '3', '4', '5',
  '6', '7', '8', '9', 'X' };

//--------------------------------------------------------------------------

// Create a histogram that is the plot of the given function.

Hist Hist::plotFunc(function<double(double)> f, string titleIn,
    int nBinIn, double xMinIn, double xMaxIn, bool logXIn) {
    Hist result(titleIn, nBinIn, xMinIn, xMaxIn, logXIn);
    if (!logXIn) {
      double dx = (xMaxIn - xMinIn) / nBinIn;
      for (double x = xMinIn + 0.5 * dx; x < xMaxIn; x += dx)
        result.fill(x, f(x));
    } else {
      double rx = pow( xMaxIn / xMinIn, 1. / nBinIn);
      for (double x = xMinIn * sqrt(rx); x < xMaxIn; x *= rx)
        result.fill(x, f(x));
    }
    return result;
  }

//--------------------------------------------------------------------------

// Book a histogram.

void Hist::book(string titleIn, int nBinIn, double xMinIn,
  double xMaxIn, bool logXIn, bool doStatsIn) {

  titleSave = titleIn;
  nBin  = nBinIn;
  if (nBinIn < 1) nBin = 1;
  if (nBinIn > NBINMAX) {
    nBin = NBINMAX;
    cout << " Warning: number of bins for histogram " << titleIn
         << " reduced to " << nBin << endl;
  }
  linX    = !logXIn;
  doStats = doStatsIn;
  xMin    = xMinIn;
  xMax    = xMaxIn;
  if (!linX && xMin < TINY) {
    xMin = TINY;
    cout << " Warning: lower x border of histogram " << titleIn
         << " increased to " << xMin << endl;
  }
  if (xMax < xMin + TINY) {
    xMax = 2. * xMin;
    cout << " Warning: upper x border of histogram " << titleIn
         << " increased to " << xMax << endl;
  }
  dx = (linX) ? (xMax - xMin) / nBin : log10(xMax / xMin) / nBin;
  res.resize(nBin);
  res2.resize(nBin);
  null();

}

//--------------------------------------------------------------------------

// Reset bin contents.

void Hist::null() {

  nFill  = 0;
  under  = 0.;
  inside = 0.;
  over   = 0.;
  for (int m = 0; m < nMoments; ++m) sumxNw[m] = 0.;
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix] = 0.;
    res2[ix] = 0.;
  }
}

//--------------------------------------------------------------------------

// Fill bin with weight.

void Hist::fill(double x, double w) {

  if (!isfinite(x) || !isfinite(w)) {nNonFinite += 1; return;}

  ++nFill;
  if (x < xMin) {under += w; return;}
  if (x > xMax) {over  += w; return;}
  int iBin = (linX) ? int( floor( (x - xMin) / dx) )
           : int( floor( log10(x / xMin) / dx) );
  if      (iBin < 0)     under += w;
  else if (iBin >= nBin) over  += w;
  else {
    res[iBin]  += w;
    res2[iBin] += w * w;
    inside     += w;
    sumxNw[0]  += w;
    sumxNw[1]  += x * w;
    // If requested, also compute unbinned higher moments.
    if (doStats) {
      double xm(x);
      for (int m = 2; m < nMoments; ++m) {
        xm *= x;
        sumxNw[m] += w * xm;
      }
    }
  }
}

//--------------------------------------------------------------------------

// Print a histogram: also operator overloading with friend.

ostream& operator<<(ostream& os, const Hist& h) {

  // Do not print empty histograms.
  if (h.nFill <= 0) {
    os << "     Histogram not shown since it is empty \n \n";
    return os;
  }

  // Write time and title.
  time_t t = time(0);
  char date[18];
  strftime(date,18,"%Y-%m-%d %H:%M",localtime(&t));
  os << "\n\n  " << date << "       " << h.titleSave << "\n\n";

  // Group bins, where required, to make printout have fewer columns.
  // Avoid overflow.
  int nGroup = 1 + (h.nBin - 1) / Hist::NCOLMAX;
  int nCol   = 1 + (h.nBin - 1) / nGroup;
  vector<double> resCol(nCol);
  for (int iCol = 0; iCol < nCol; ++iCol) {
    resCol[iCol] = 0.;
    for (int ix = nGroup * iCol; ix < min( h.nBin, nGroup * (iCol + 1)); ++ix)
      resCol[iCol] += h.res[ix];
    resCol[iCol] = max( -Hist::LARGE, min( Hist::LARGE, resCol[iCol] ) );
  }

  // Find minimum and maximum bin content.
  double yMin = resCol[0];
  double yMax = resCol[0];
  for (int iCol = 1; iCol < nCol; ++iCol) {
    if (resCol[iCol] < yMin) yMin = resCol[iCol];
    if (resCol[iCol] > yMax) yMax = resCol[iCol];
  }

  // Determine scale and step size for y axis.
  if (yMax - yMin > Hist::NLINES * DYAC[0] * 1e-9) {
    if (yMin > 0. && yMin < Hist::SMALLFRAC * yMax) yMin = 0.;
    if (yMax < 0. && yMax > Hist::SMALLFRAC * yMin) yMax = 0.;
    int iPowY = int(floor( log10(yMax - yMin) ));
    if (yMax - yMin < Hist::NLINES * DYAC[0] * pow(10.,iPowY))
      iPowY = iPowY - 1;
    if (yMax - yMin > Hist::NLINES * DYAC[9] * pow(10.,iPowY))
      iPowY = iPowY + 1;
    double nLinePow = Hist::NLINES * pow(10.,iPowY);
    double delY = DYAC[0];
    for (int idel = 0; idel < 9; ++idel)
      if (yMax - yMin >= nLinePow * DYAC[idel]) delY = DYAC[idel+1];
    double dy = delY * pow(10.,iPowY);

    // Convert bin contents to integer form; fractional fill in top row.
    vector<int> row(nCol);
    vector<int> frac(nCol);
    for (int iCol = 0; iCol < nCol ; ++iCol) {
      double cta = abs(resCol[iCol]) / dy;
      row[iCol] = int(cta + 0.95);
      if(resCol[iCol] < 0.) row[iCol] = - row[iCol];
      frac[iCol] = int(10. * (cta + 1.05 - floor(cta + 0.95)));
    }
    int rowMin = int(abs(yMin)/dy + 0.95);
    if ( yMin < 0) rowMin = - rowMin;
    int rowMax = int(abs(yMax)/dy + 0.95);
    if ( yMax < 0) rowMax = - rowMax;

    // Print histogram row by row.
    os << fixed << setprecision(2);
    for (int iRow = rowMax; iRow >= rowMin; iRow--) if (iRow != 0) {
      os << "  " << setw(10) << iRow*delY << "*10^"
         << setw(2) << iPowY << "  ";
      for (int iCol = 0; iCol < nCol ; ++iCol) {
        if (iRow == row[iCol])                  os << NUMBER[frac[iCol]];
        else if (iRow * (row[iCol] - iRow) > 0) os << NUMBER[10];
        else                                    os << " ";
      }
      os << "\n";
    }
    os << "\n";

    // Print sign and value of bin contents
    double maxim = log10(max(yMax, -yMin));
    int iPowBin = int(floor(maxim + 0.0001));
    os << "          Contents  ";
    for (int iCol = 0; iCol < nCol ; ++iCol) {
      if (resCol[iCol] < - pow(10., iPowBin - 4)) os << "-";
      else os << " ";
      row[iCol] = int(abs(resCol[iCol]) * pow(10., 3 - iPowBin) + 0.5);
    } os << "\n";
    for (int iRow = 3; iRow >= 0; iRow--) {
      os << "            *10^" << setw(2) << iPowBin + iRow - 3 << "  ";
      int mask = int( pow(10., iRow) + 0.5);
      for (int iCol = 0; iCol < nCol ; ++iCol) {
        os << NUMBER[(row[iCol] / mask) % 10];
      } os << "\n";
    } os << "\n";

    // Print sign and value of lower bin edge.
    maxim = (h.linX) ? log10( max( -h.xMin, h.xMax - h.dx))
          :  log10( max( -log10(h.xMin), log10(h.xMax) - h.dx ) );
    int iPowExp = int(floor(maxim + 0.0001));
    os << "          Low edge  ";
    for (int iCol = 0; iCol < nCol ; ++iCol) {
      double edgeNow = (h.linX) ? h.xMin + iCol * nGroup * h.dx
        : log10(h.xMin) + iCol * nGroup * h.dx;
      os << ( ( edgeNow < - pow(10., iPowExp - 3) ) ? "-" : " " );
      row[iCol] = int( abs(edgeNow) * pow(10., 2 - iPowExp) + 0.5 );
    } os << "\n";
    for (int iRow = 2; iRow >= 0; iRow--) {
      os << "            *10^" << setw(2) << iPowExp + iRow - 2 << "  ";
      int mask = int( pow(10., iRow) + 0.5);
      for (int iCol = 0; iCol < nCol ; ++iCol)
        os << NUMBER[(row[iCol] / mask) % 10];
      os << "\n";
    } os << "\n";

  // Print explanation if histogram cannot be shown.
  } else os << "     Histogram not shown since lowest value" << scientific
       << setprecision(4) << setw(12) << yMin << " and highest value"
       << setw(12) << yMax << " are too close \n \n";

  // Mean determines overall output format for moments.
  double xMean = h.getXMean(false);
  bool doExp   = false;
  int  prec    = 3;
  if ( abs(xMean) >= 10000. || abs(xMean) <= 1. ) doExp = true;
  else if ( abs(xMean) >= 1000. ) prec = 1;
  else if ( abs(xMean) >= 100. ) prec = 2;

  // Fixed or Scientific for ranges.
  bool doExpRange = ( abs(h.xMin) >= 1e5 || abs(h.xMax) >= 1e5 ||
    (abs(h.xMin) != 0.0 && abs(h.xMin) <= 1e-2) ||
    (abs(h.xMin) != 0.0 && abs(h.xMin) <= 1e-2) );
  int rangePrec = 3;
  if (!doExpRange && max(abs(h.xMin),abs(h.xMax)) <= 1e-1) rangePrec = 5;
  else if (!doExpRange && max(abs(h.xMin),abs(h.xMax)) <= 1.) rangePrec = 4;

  // Fixed or Scientific for contents.
  bool contExp = ( abs(h.inside + h.under + h.over) >= 1e5 ||
    (abs(h.inside) <= 1e-2) );
  int contPrec   = 3;
  if (!contExp && abs(h.inside) <= 1e-1) contPrec = 5;
  else if (!contExp && abs(h.inside) <= 1.) contPrec = 4;

  // doStats and nEff determine if uncertainties should be printed.
  double nEff  = h.getNEffective();
  bool   doErr = (h.doStats && nEff > 0.);
  string pad   = "   ";

  // Fixed or Scientific for nEffective.
  bool neffExp = ( abs(nEff) >= 1e5 || (abs(nEff) <= 1e-2) );
  int nEffPrec = 3;
  if (!neffExp && abs(nEff) <= 1e-1) nEffPrec = 5;
  else if (!neffExp && abs(nEff) <= 1.) nEffPrec = 4;

  // ------------------------------------------------------------------------
  // First line.

  // Number of entries.
  os << fixed;
  os << pad << "Entries  =" << setw(10) << h.nFill;

  // xMin.
  if (doExpRange) os << scientific << setprecision(3);
  else os << fixed << setprecision(rangePrec);
  os << pad << "xMin =" << setw(10) << h.xMin;

  // Underflow
  if (contExp) os << scientific << setprecision(3);
  else os << fixed << setprecision(contPrec);
  os << pad << "Underflow =" << setw(10) << h.under;

  // Mean.
  if (doExp) os << scientific << setprecision(3);
  else os << fixed << setprecision(prec);
  os << pad << "Mean   =" << setw(10) << xMean;
  if (doErr) {
    double xMeanErr = h.getXMeanErr(false);
    if (doExp|| xMeanErr > 10 * abs(xMean)) os << setprecision(1);
    os << " +-" << setw(7) << xMeanErr;
  }

  // RMS.
  double xRMS      = h.getXRMS(false);
  if (doExp) os << scientific << setprecision(3);
  else os << fixed << setprecision(prec);
  if (!doErr) os << pad << "RMS  =" << setw(10) << xRMS;
  else {
    os << pad << "RMS  =" << setw(10) << xRMS;
    double xRMSErr = h.getXRMSErr(false);
    if (doExp || xRMSErr > 10 * abs(xRMS)) os << setprecision(1);
    os << " +-" << setw(7) << xRMSErr;
  }

  // End of line.
  os << endl;

  // ------------------------------------------------------------------------
  // Second line.

  // Sum of Weights (inside histogram range).
  if (contExp) os << scientific << setprecision(3);
  else os << fixed << setprecision(contPrec);
  os << pad << "SumW(in) =" << setw(10) << h.inside;


  // xMax.
  if (doExpRange) os << scientific << setprecision(3);
  else os << fixed << setprecision(rangePrec);
  os << pad << "xMax =" << setw(10) << h.xMax;

  // Overflow
  if (contExp) os << scientific << setprecision(3);
  else os << fixed << setprecision(contPrec);
  os << pad << "Overflow  =" << setw(10) << h.over;

  // Median (excluding underflow and overflow bins).
  double xMedian    = h.getXMedian(false);
  if (doExp) os << scientific << setprecision(3);
  else os << fixed << setprecision(prec);
  os << pad << "Median =" << setw(10) << xMedian;
  if (doErr) {
    double xMedianErr = h.getXMedianErr(false);
    if (doExp || xMedianErr > 10 * abs(xMedian)) os << setprecision(1);
    os << " +-" << setw(7) << xMedianErr;
  }

  // nEff: Statistical power = effective number of unweighted entries.
  // If nEff ~ h.inside, use same precision as for SumW.
  string var = (h.doStats ? "nEffective =   " : "nEff =");
  if (nEff <= 0.) {
    os << pad << var << setw(10) << "N/A";
  } else if (h.doStats) {
    os << pad << var;
    if (neffExp) os << scientific << setprecision(3);
    else os << fixed << setprecision(nEffPrec);
    os << setw(10) << nEff;
  } else {
    os << pad << var;
    if (neffExp) os << scientific << setprecision(3);
    else os << fixed << setprecision(nEffPrec);
    os << setw(10) << nEff;
  }

  // Return to standard PYTHIA format.
  os << scientific << setprecision(3) << endl;
  return os;
}

//--------------------------------------------------------------------------

// Print histogram contents as a table (e.g. for Gnuplot).

void Hist::table(ostream& os, bool printOverUnder, bool xMidBin,
  bool printError) const {

  // Print histogram vector bin by bin, with mean x as first column.
  os << scientific << setprecision(4);
  double xBeg = (xMidBin) ? xMin + 0.5 * dx : xMin;
  if (!linX && xMidBin) xBeg = xMin * pow( 10., 0.5 * dx);
  if (printOverUnder) {
    os << setw(12) << (linX ? xBeg - dx : xBeg * pow(10., -dx))
       << setw(12) << under;
    if (printError) os << setw(12) << 0.0 << "\n";
    else os << "\n";
  }
  for (int ix = 0; ix < nBin; ++ix) {
    os << setw(12) << (linX ? xBeg + ix * dx : xBeg * pow(10., ix * dx))
       << setw(12) << res[ix];
      if (printError) os << setw(12) << sqrtpos(res2[ix]) << "\n";
      else os << "\n";
  }
  if (printOverUnder) {
    os << setw(12) << (linX ? xBeg + nBin * dx : xBeg * pow(10., nBin * dx))
       << setw(12) << over;
    if (printError) os << setw(12) << 0.0 << "\n";
    else os << "\n";
  }

}

//--------------------------------------------------------------------------

// Print histogram contents as a table, in Rivet's *.dat style.

void Hist::rivetTable(ostream& os, bool printError) const {

  // Print histogram vector bin by bin, with x range in first two columns
  // and +- error in last two.
  os << scientific << setprecision(4);
  double xBeg = xMin;
  double xEnd = (linX) ? xMin + dx : xMin * pow(10., dx);
  for (int ix = 0; ix < nBin; ++ix) {
    double err = (printError) ? sqrtpos(res2[ix]) : 0.0;
    os << setw(12) << (linX ? xBeg + ix * dx : xBeg * pow(10., ix * dx))
       << setw(12) << (linX ? xEnd + ix * dx : xEnd * pow(10., ix * dx))
       << setw(12) << res[ix] << setw(12) << err << setw(12) << err << "\n";
  }

}

//--------------------------------------------------------------------------

// Print histogram contents as a table, as appropriate for Pyplot.

void Hist::pyplotTable(ostream& os, bool isHist, bool printError) const {

  // Set precision.
  os << scientific << setprecision(4);

  // For plotting as a histogram one needs bin edges as last column.
  double xBeg = (linX) ? xMin + 0.5 * dx : xMin * pow( 10., 0.5 * dx);
  double xNow, xEdge;
  for (int ix = 0; ix < nBin; ++ix) {
    xNow  = (linX) ? xBeg + ix * dx : xBeg * pow(10., ix * dx);
    xEdge = (linX) ? xMin + ix * dx : xMin * pow(10., ix * dx);
    os << setw(12) << xNow << setw(12) << res[ix];
    if (isHist) os << setw(12) << xEdge;
    if (printError) os << setw(12) << sqrtpos(res2[ix]);
    os << "\n";
  }

  // And also an extra no-weights line to give final upper bin edge.
  if (isHist) {
    double xEnd = (linX) ? xMax - 0.5 * dx : xMax * pow( 10., -0.5 * dx);
    os << setw(12) << xEnd << setw(12) << 0. << setw(12) << xMax;
    if (printError) os << setw(12) << 0.;
    os << "\n";
  }

}

//--------------------------------------------------------------------------

// Fill histogram contents from a table, e.g. written by table() above.

void Hist::fillTable(istream& is) {

  // Read in one line at a time.
  string line;
  double xVal, yVal;
  while ( getline(is, line) ) {
    // Read contents of the line and fill in histogram.
    istringstream splitLine(line);
    splitLine >> xVal >> yVal;
    fill( xVal, yVal);
  }
}

//--------------------------------------------------------------------------

// Print a table out of two histograms with same x axis  (e.g. for Gnuplot).

void table(const Hist& h1, const Hist& h2, ostream& os, bool printOverUnder,
  bool xMidBin) {

  // Require histogram x axes to agree.
  int nBin  = h1.nBin;
  double dx = h1.dx;
  if (nBin != h2.nBin || abs(h1.xMin - h2.xMin) > Hist::TOLERANCE * dx
    || abs(h1.xMax - h2.xMax) > Hist::TOLERANCE * dx || h1.linX != h2.linX)
    return;

  // Print histogram vectors bin by bin, with mean x as first column.
  os << scientific << setprecision(4);
  double xBeg = (xMidBin) ? h1.xMin + 0.5 * dx : h1.xMin;
  if (!h1.linX && xMidBin) xBeg = h1.xMin * pow(10., 0.5 * dx);
  if (printOverUnder)
    os << setw(12) << (h1.linX ? xBeg - dx : xBeg * pow(10., -dx))
       << setw(12) << h1.under << setw(12) << h2.under << "\n";
  for (int ix = 0; ix < nBin; ++ix)
    os << setw(12) << (h1.linX ? xBeg + ix * dx : xBeg * pow(10., ix * dx))
       << setw(12) << h1.res[ix] << setw(12) << h2.res[ix] << "\n";
  if (printOverUnder)
    os << setw(12) << (h1.linX ? xBeg + nBin * dx : xBeg * pow(10., nBin * dx))
       << setw(12) << h1.over << setw(12) << h2.over << "\n";

}

void table(const Hist& h1, const Hist& h2, string fileName,
  bool printOverUnder, bool xMidBin) {

  ofstream streamName(fileName.c_str());
  table( h1, h2, streamName, printOverUnder, xMidBin);

}

//--------------------------------------------------------------------------

// Compute mean value in X, by default unbinned.

double Hist::getXMean(bool unbinned) const {
  if (unbinned) return sumxNw[1] / max(sumxNw[0], Hist::TINY);
  else {
    double wtSum(0.);
    double wtxSum(0.);
    for (int ix=0; ix < nBin ; ++ix) {
      double wta  = abs(res[ix]);
      double x = (linX) ? xMin + (ix + 0.5) * dx
        : xMin * pow( 10., (ix + 0.5) * dx);
      wtSum   = wtSum   + wta;
      wtxSum  = wtxSum  + wta * x;
    }
    return wtxSum / max(wtSum, Hist::TINY);
  }
}

//--------------------------------------------------------------------------

// Compute median in X (with linear interpolation inside median bin).
// Note: absolute values of weights are used.

double Hist::getXMedian(bool includeOverUnder) const {
  double wtSumNow = 0.;
  double wtSumTot = 0.;
  for (int ix = 0; ix < nBin ; ++ix) {
    wtSumTot += abs(res[ix]);
  }
  // Include underflow and overflow bins in median definition?
  if (includeOverUnder) {
    wtSumTot += abs(over) + abs(under);
    wtSumNow = abs(under);
    // If excess bins contain more than half, return low or high edge.
    if (abs(under) > 0.5 * wtSumTot) return xMin;
    else if (abs(over) > 0.5 * wtSumTot) return xMax;
  }
  for (int ix = 0; ix < nBin ; ++ix) {
    double wtSumOld = wtSumNow;
    wtSumNow += abs(res[ix]);
    if (wtSumNow > 0.5 * wtSumTot) {
      double frac = (0.5*wtSumTot - wtSumOld)/(wtSumNow - wtSumOld);
      return (linX) ? xMin + (ix + frac) * dx
                    : xMin * pow( 10., (ix + frac) * dx);
    }
  }
  // Should not arrive here. Safe return in any case.
  return 0.;
}

//--------------------------------------------------------------------------

// Compute statistical error on median in X ( + bin granularity estimate).
// Note: absolute values of weights are used.

double Hist::getXMedianErr(bool includeOverUnder) const {

  // If no statistics, cannot compute error.
  if (getNEffective() <= 0.) return 0.;

  double xMedian = getXMedian(includeOverUnder);
  if (xMedian <= xMin || xMedian >= xMax) return 0.;
  double wtSumTot = max( abs(sumxNw[0]), Hist::TINY );

  // Include underflow and overflow bins in median definition?
  if (includeOverUnder) wtSumTot += abs(over) + abs(under);
  // Laplace's formula for variance of median: 1/(4 nEff f(xmedian)^2).
  int iBin = int( (linX) ? (xMedian - xMin)/dx : log10( xMedian / xMin )/dx);
  double fMedian  = (linX) ? abs(res[iBin]) / dx / wtSumTot :
    abs(res[iBin]) / pow( 10. , dx ) / wtSumTot;
  double statFac    = sqrtpos( 1. / max(getNEffective(), Hist::TINY) );
  double xMedianErr = statFac / 2. / max( fMedian, Hist::TINY );

  // Add uncertainty estimate from bin granularity from binned vs unbinned <x>.
  double xBinErr = getXMean(true) - getXMean(false);
  return sqrtpos( xMedianErr * xMedianErr + xBinErr * xBinErr );

}

//--------------------------------------------------------------------------

// Get RMN = n'th root of n'th moment about the mean.
// Lowest ones:
//  n = 2: RMS = sigma = width of (peaked) distribution.
//  n = 3: cube root of third moment about the mean < ( x - <x> )^3 >.
//         Sensitive to one-sided tails, as is characteristic of many
//         particle-physics distributions. Related to skewness.
//  n = 4: fourth root of fourth moment about the mean < ( x - <x> )^4 >.
//         Sensitive to double-sided long tails (eg BW vs Gauss) and
//         4th-power asymptotics for single-sided ones. Related to kurtosis.
//  ...

double Hist::getXRMN(int n, bool unbinned) const {
  // Compute from unbinned stored sums, or from binned histogram distributions.
  if (unbinned && n >= 1 && n <= 6) {
    // 1st central moment is always zero (Mean - Mean).
    if (n == 1) return 0.;
    // Unbinned RMN implemented for N = 2, 3, 4, 5, 6.
    // One-pass algorithm, may suffer from loss of numerical accuracy
    // if there are catastrophic cancellations among terms.
    double sumw = max(sumxNw[0], Hist::TINY);
    double mean = sumxNw[1] / sumw;
    if (n == 2) return sqrtpos( sumxNw[2]/sumw - mean * mean );
    else if (n == 3) {
      return cbrt( (sumxNw[3] - 3. * mean * sumxNw[2]) / sumw
      + 2 * mean * mean * mean );
    } else if (n == 4) {
      double m4d = (sumxNw[4] - 4. * mean * sumxNw[3]
        + 6. * mean * mean * sumxNw[2]) / sumw - 3. * pow4(mean);
      return pow( max(0., m4d), 0.25);
    } else if (n == 5) {
      double m5d = (sumxNw[5] - 5. * mean * sumxNw[4]
        + 10. * mean * mean * sumxNw[3] - 10 * mean * mean * mean * sumxNw[2] )
        / sumw + 4. * pow(mean, 5);
      if (m5d < 0.) return pow( abs(m5d) , 0.2);
      else return pow(m5d, 0.2);
    } else if (n == 6) {
      double m6d = (sumxNw[6] - 6. * mean * sumxNw[5]
        + 15. * mean * mean * sumxNw[4] - 20 * mean * mean * mean * sumxNw[3]
        + 15 * pow4(mean) * sumxNw[2] ) / sumw - 5. * pow5(mean);
      return pow( max(0., m6d), 1./6.);
    }
    // No further unbinned moments defined.
    else return 0.;
  } else {
    // Binned: two-pass algorithm should be numerically more safe, but
    // may be limited by bin granularity. First pass = compute binned mean.
    // Note: absolute values of weights are used.
    double mean = getXMean(false);
    // Second pass. Compute n'th moment about the mean (for any n).
    double mnd  = 0.;
    double sumw = 0.;
    for (int ix=0; ix < nBin ; ++ix) {
      double wta  = abs(res[ix]);
      sumw += wta;
      double x = (linX) ? xMin + (ix + 0.5) * dx
        : xMin * pow( 10., (ix + 0.5) * dx);
      if ( n == 2) mnd += wta * (x - mean) * (x - mean);
      else mnd += wta * pow(x - mean, n);
    }
    mnd  /= max(sumw, Hist::TINY);
    if ( n == 2 ) return sqrtpos( mnd );
    else if ( n == 3 ) return cbrt( mnd );
    else if ( n == 4 ) return sqrt( sqrtpos( mnd) );
    // General moments, for any n. Odd ones can be negative.
    else if ( mnd < 0 && ( n % 2 ) == 1 ) return -pow( abs(mnd), 1./n );
    else return pow( max(mnd, 0.), 1./n );
  }
  // No further moments defined.
  return 0.;
}

//--------------------------------------------------------------------------

// Error on mean value.

double Hist::getXMeanErr(bool unbinned) const {

  // If no statistics, cannot compute error.
  if (getNEffective() <= 0.) return 0.;

  // Get RMS and multiply by statistical power.
  double sigma = getXRMN(2, unbinned);
  double xMeanErr2 = pow2(sigma) / max( getNEffective(), Hist::TINY );

  // If binned, add estimated bin granularity error <x>_binned - <x>_unbinned.
  if (!unbinned) {
    double xBinErr = getXMean(true) - getXMean(false);
    xMeanErr2 += xBinErr * xBinErr;
  }

  return sqrtpos( xMeanErr2 );
}

//--------------------------------------------------------------------------

// Error on n'th root of n'th moment about the mean.
// Computed using sigma2( <x^n> ) = < (x^n - <x>^n)^2 > / nEff
// Calculation always uses binned values, but if flag unbinned == false,
// a granularity uncertainty measure is added.
// Note: absolute values of weights are used.

double Hist::getXRMNErr(int n, bool unbinned) const {

  // Check if error can be computed.
  double nEffective = getNEffective();
  double xRMN      = getXRMN(n, false);
  if (nEffective <= 0. || xRMN == 0.) return 0.;

  // First pass: compute binned mean.
  double mean = getXMean(false);

  // Second pass: compute < (x^n - <x>^n)^2 >
  double msdn = 0.;
  double sumw = 0.;
  for (int ix=0; ix < nBin ; ++ix) {
    double wta = abs(res[ix]);
    double x = (linX) ? xMin + (ix + 0.5) * dx
      : xMin * pow( 10., (ix + 0.5) * dx);
    sumw += wta;
    msdn += wta * pow2( pow(x, n) - pow(mean, n) );
  }
  msdn /= max(sumw, Hist::TINY);
  // Uncertainty on n'th root = 1/n * relative uncertainty.
  double xRMNErr2 = msdn / (n * n) / max(nEffective, Hist::TINY)
    / pow( abs(xRMN), 2*n - 2 );

  // If the moment itself was computed from binned distribution,
  // add granularity uncertainty.
  if (!unbinned) {
    double xBinErr = getXRMN(n, true) - xRMN;
    xRMNErr2 += xBinErr * xBinErr;
  }

  return sqrtpos(xRMNErr2);
}

//--------------------------------------------------------------------------

// Get content of specific bin.
// Special values are bin 0 for underflow and bin nBin+1 for overflow.
// All other bins outside proper histogram range return 0.

double Hist::getBinContent(int iBin) const {

  if (iBin > 0 && iBin <= nBin) return res[iBin - 1];
  else if (iBin == 0)           return under;
  else if (iBin == nBin + 1)    return over;
  else                          return 0.;

}

// Return the lower edge of the bin.

double Hist::getBinEdge(int iBin) const {

  if (iBin > 0 && iBin <= nBin + 1)
    return linX ? xMin + (iBin - 1) * dx : xMin * pow(10., (iBin - 1) * dx);
  else return numeric_limits<double>::quiet_NaN();

}

// Return the width of the bin.

double Hist::getBinWidth(int iBin) const {

  if (iBin > 0 && iBin <= nBin)
    return linX ? dx : xMin * (pow(10., dx) - 1) * pow(10., (iBin - 1) * dx);
  else return numeric_limits<double>::infinity();

}

//--------------------------------------------------------------------------

// Return bin contents and edges.

vector<double> Hist::getBinContents() const {return res;}

vector<double> Hist::getBinEdges() const {

  vector<double> edges(nBin + 1);
  for (int ix = 0; ix <= nBin; ++ix) edges[ix] = getBinEdge(ix + 1);
  return edges;

}

//--------------------------------------------------------------------------

// Check whether another histogram has same size and limits.

bool Hist::sameSize(const Hist& h) const {

  if (nBin == h.nBin && abs(xMin - h.xMin) < TOLERANCE * dx &&
    abs(xMax - h.xMax) < TOLERANCE * dx) return true;
  else return false;

}

//--------------------------------------------------------------------------

// Take an arbitrary function of bin contents.

void Hist::takeFunc(function<double(double)> func) {

  // Result for sumxNw will use the new binned values.
  for (int m = 0; m < nMoments; ++m) sumxNw[m] = 0;
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix] = func(res[ix]);
    double x = (linX) ? xMin + (ix + 0.5) * dx
      : xMin * pow( 10., (ix + 0.5) * dx);
    sumxNw[0] += res[ix];
    sumxNw[1] += x * res[ix];
    for (int m = 2; m < nMoments; ++m) sumxNw[m] += res[ix] * pow(x, m);
  }
  under  = func(under);
  inside = func(inside);
  over   = func(over);

}

//--------------------------------------------------------------------------

// Take 10-logarithm or natural logarithm of contents bin by bin.

void Hist::takeLog(bool tenLog) {

  // Find smallest positive bin content, and put min a bit below.
  double yMin = Hist::LARGE;
  for (int ix = 0; ix < nBin; ++ix)
    if (res[ix] > Hist::TINY && res[ix] < yMin ) yMin = res[ix];
  yMin *= 0.8;

  // Take base 10 or natural logarithm bin by bin, but ensure positivity.
  takeFunc([yMin, tenLog](double x) {
      return tenLog ? log10( max( yMin, x) ) : log( max( yMin, x) );});

}

//--------------------------------------------------------------------------

// Take square root of contents bin by bin; set 0 for negative content.

void Hist::takeSqrt() {takeFunc(sqrtpos);}

//--------------------------------------------------------------------------

// Normalize bin contents to given sum, by default including overflow bins.

void Hist::normalize(double f, bool overflow) {
  *this *= f / (overflow ? (under + inside + over) : inside);
}

//--------------------------------------------------------------------------

// Normalize bin contents to given area, by default including overflow bins.

void Hist::normalizeIntegral(double f, bool overflow) {
  normalizeSpectrum((overflow ? (under + inside + over) : inside) / f);
}

//--------------------------------------------------------------------------

// Scale each bin content by 1 / (wtSum * bin width).

void Hist::normalizeSpectrum(double wtSum) {
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix] /= wtSum * getBinWidth(ix + 1);
    res2[ix] /= pow2( wtSum * getBinWidth(ix + 1) );
  }
  inside /= wtSum;
  over /= wtSum;
  under /= wtSum;
}

//--------------------------------------------------------------------------

// Add histogram to existing one.

Hist& Hist::operator+=(const Hist& h) {
  if (!sameSize(h)) return *this;
  nFill  += h.nFill;
  under  += h.under;
  inside += h.inside;
  over   += h.over;
  doStats = doStats && h.doStats;
  for (int m = 0; m < nMoments; ++m) sumxNw[m] += h.sumxNw[m];
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix]  += h.res[ix];
    res2[ix] += h.res2[ix];
  }
  return *this;
}

//--------------------------------------------------------------------------

// Subtract histogram from existing one.

Hist& Hist::operator-=(const Hist& h) {
  if (!sameSize(h)) return *this;
  nFill  += h.nFill;
  under  -= h.under;
  inside -= h.inside;
  over   -= h.over;
  doStats = doStats && h.doStats;
  for (int m = 0; m < nMoments; ++m) sumxNw[m]  -= h.sumxNw[m];
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix]  -= h.res[ix];
    res2[ix] += h.res2[ix];
  }
  return *this;
}

//--------------------------------------------------------------------------

// Multiply existing histogram by another one.
// (Assumes Gaussian uncertainty propagation.)

Hist& Hist::operator*=(const Hist& h) {
  if (!sameSize(h)) return *this;
  nFill  += h.nFill;
  under  *= h.under;
  inside *= h.inside;
  over   *= h.over;
  doStats = false;
  // Result for sumxNw has to use binned values.
  for (int m = 0; m < nMoments; ++m) sumxNw[m] = 0.;
  for (int ix = 0; ix < nBin; ++ix) {
    res2[ix] = (abs(res[ix]) < Hist::TINY || abs(h.res[ix]) < Hist::TINY) ? 0 :
      pow2(res[ix]*h.res[ix])*(
        res2[ix]/pow2(res[ix]) + h.res2[ix]/pow2(h.res[ix]));
    res[ix] *= h.res[ix];
    double x = (linX) ? xMin + (ix + 0.5) * dx
      : xMin * pow( 10., (ix + 0.5) * dx);
    sumxNw[0] += res[ix];
    sumxNw[1] += x * res[ix];
    for (int m = 2; m < nMoments; ++m) sumxNw[m] += res[ix] * pow(x, m);
  }
  return *this;
}

//--------------------------------------------------------------------------

// Divide existing histogram by another one.
// (Assumes Gaussian uncertainty propagation.)

Hist& Hist::operator/=(const Hist& h) {
  if (!sameSize(h)) return *this;
  nFill  += h.nFill;
  under   = (abs(h.under) < Hist::TINY) ? 0. : under/h.under;
  inside  = (abs(h.inside) < Hist::TINY) ? 0. : inside/h.inside;
  over    = (abs(h.over) < Hist::TINY) ? 0. : over/h.over;
  doStats = false;
  // Result for sumxNw has to use binned values.
  for (int m = 0; m < nMoments; ++m) sumxNw[m] = 0.;
  for (int ix = 0; ix < nBin; ++ix) {
    res2[ix] = (abs(res[ix]) < Hist::TINY || abs(h.res[ix]) < Hist::TINY) ? 0 :
      pow2(res[ix]/h.res[ix])*(
        res2[ix]/pow2(res[ix]) + h.res2[ix]/pow2(h.res[ix]));
    res[ix]  = (abs(h.res[ix]) < Hist::TINY) ? 0. : res[ix]/h.res[ix];
    double x = (linX) ? xMin + (ix + 0.5) * dx
      : xMin * pow( 10., (ix + 0.5) * dx);
    sumxNw[0] += res[ix];
    sumxNw[1] += x * res[ix];
    for (int m = 2; m < nMoments; ++m) sumxNw[m] += res[ix] * pow(x, m);
  }
  return *this;
}

//--------------------------------------------------------------------------

// Add constant offset to histogram.
// Squared weights treated as filling each bin with weight = f,
// For linear x, moments treated as filling with continuous f/dx.
// for logx, moments treated as filling with f at (log) centre of each bin.

Hist& Hist::operator+=(double f) {
  under  += f;
  inside += nBin * f;
  over   += f;
  sumxNw[0] += nBin * f;
  if (linX) {
    // Moments (for linX): analytical integral_xMin^xMax f/dx * x^m dx.
    double xMinM(xMin);
    double xMaxM(xMax);
    for (int m = 1; m < nMoments; ++m) {
      xMinM *= xMin;
      xMaxM *= xMax;
      sumxNw[m] += f * (xMaxM - xMinM) / (m + 1) / dx;
    }
  }
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix]  += f;
    res2[ix] += f * f;
    if (!linX) {
      // Moments (for logX): numerical sum of f * x^m.
      double x = xMin * pow( 10., (ix + 0.5) * dx);
      double xm(1.);
      for (int m = 1; m < nMoments; ++m) {
        xm *= x;
        sumxNw[m] += f * xm;
      }
    }
  }
  return *this;
}

//--------------------------------------------------------------------------

// Subtract constant offset from histogram.
// Squared weights treated as filling each bin with weight = -f.
// Moments treated as filling with continuous -f/dx.

Hist& Hist::operator-=(double f) {
  under  -= f;
  inside -= nBin * f;
  over   -= f;
  sumxNw[0] -= nBin * f;
  if (linX) {
    // Moments (for linX): analytical integral_xMin^xMax f/dx * x^m dx.
    double xMinM(xMin);
    double xMaxM(xMax);
    for (int m = 1; m < nMoments; ++m) {
      xMinM *= xMin;
      xMaxM *= xMax;
      sumxNw[m] -= f * (xMaxM - xMinM) / (m + 1) / dx;
    }
  }
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix]  -= f;
    res2[ix] -= f * f;
    if (!linX) {
      // Moments (for logX): numerical sum of f * x^m.
      double x = xMin * pow( 10., (ix + 0.5) * dx);
      sumxNw[1] -= f * x;
      double xm(x);
      for (int m = 2; m < nMoments; ++m) {
        xm *= x;
        sumxNw[m] -= f * xm;
      }
    }
  }
  return *this;
}

//--------------------------------------------------------------------------

// Multiply histogram by constant.

Hist& Hist::operator*=(double f) {
  under  *= f;
  inside *= f;
  over   *= f;
  for (int m = 0; m < nMoments; ++m) sumxNw[m] *= f;
  for (int ix = 0; ix < nBin; ++ix) {
    res[ix]  *= f;
    res2[ix] *= f * f;
  }
  return *this;
}

//--------------------------------------------------------------------------

// Divide histogram by constant.

Hist& Hist::operator/=(double f) {
  if (abs(f) > Hist::TINY) {
    under  /= f;
    inside /= f;
    over   /= f;
    for (int m = 0; m < nMoments; ++m) sumxNw[m] /= f;
    for (int ix = 0; ix < nBin; ++ix) {
      res[ix] /= f;
      res2[ix] /= f * f;
    }
  // Set empty contents when division by zero.
  } else {
    under  = 0.;
    inside = 0.;
    over   = 0.;
    for (int m = 0; m < nMoments; ++m) sumxNw[m] = 0.;
    for (int ix = 0; ix < nBin; ++ix) {
      res[ix]  = 0.;
      res2[ix] = 0.;
    }
  }
  return *this;
}

Hist Hist::operator+(double f) const {
  Hist h = *this; return h += f;}

Hist Hist::operator+(const Hist& h2) const {
  Hist h = *this; return h += h2;}

Hist Hist::operator-(double f) const {
  Hist h = *this; return h -= f;}

Hist Hist::operator-(const Hist& h2) const {
  Hist h = *this; return h -= h2;}

Hist Hist::operator*(double f) const {
  Hist h = *this; return h *= f;}

Hist Hist::operator*(const Hist& h2) const {
  Hist h = *this; return h *= h2;}

Hist Hist::operator/(double f) const {
  Hist h = *this; return h /= f;}

Hist Hist::operator/(const Hist& h2) const {
  Hist h = *this; return h /= h2;}

//--------------------------------------------------------------------------

// Implementation of operator overloading with friends.

Hist operator+(double f, const Hist& h1) {
  Hist h = h1; return h += f;}

Hist operator-(double f, const Hist& h1) {
  Hist h    = h1;
  h.under   = f - h1.under;
  h.inside  = h1.nBin * f - h1.inside;
  h.over    = f - h1.over;
  h.doStats = h1.doStats;
  for (int m = 0; m < Hist::nMoments; ++m) h.sumxNw[m]  = f - h1.sumxNw[m];
  for (int ix = 0; ix < h1.nBin; ++ix) {
    h.res[ix]  = f - h1.res[ix];
    h.res2[ix] = h1.res2[ix];
  }
  return h;}

Hist operator*(double f, const Hist& h1) {
  Hist h = h1; return h *= f;}

Hist operator/(double f, const Hist& h1) {
  Hist h = h1;
  h.under   = (abs(h1.under)  < Hist::TINY) ? 0. :  f/h1.under;
  h.inside  = (abs(h1.inside) < Hist::TINY) ? 0. :  f/h1.inside;
  h.over    = (abs(h1.over)   < Hist::TINY) ? 0. :  f/h1.over;
  h.doStats = h1.doStats;
  for (int m = 0; m < Hist::nMoments; ++m)
    h.sumxNw[m] = (abs(h1.sumxNw[m]) < Hist::TINY) ? 0. :  f/h1.sumxNw[m];
  for (int ix = 0; ix < h1.nBin; ++ix) {
    h.res[ix]  = (abs(h1.res[ix]) < Hist::TINY) ? 0. : f/h1.res[ix];
    h.res2[ix] = pow2(f) * h1.res2[ix];
  }
  return h;
}

//==========================================================================

// HistPlot class.
// Writes a Python program that can generate PDF plots from Hist histograms.

//--------------------------------------------------------------------------

//  Generate the Python code for plotting a frame.

void HistPlot::plot( bool logY, bool logX, bool userBorders) {

  // Start new file or add to existing one.
  if (frameName != "" && frameName != framePrevious) {
    if (nPDF > 0) toPython << "pp.close()" << endl;
    ++nPDF;
    fileName = frameName;
    toPython << "pp   = PdfPages('" << fileName << ".pdf')" << endl;
    nFrame = 1;
    nTable = 0;
  } else {
    ++nFrame;
  }
  toPython << "tmp" << nFrame << " = plt.figure(" << nFrame << ")" << endl;
  toPython << "tmp" << nFrame << ".set_size_inches(" << fixed
    << setprecision(2) << xSize << "," << ySize << ")" << endl;

  // Loop through the vector of histograms.
  double xMinTot = 1e10, xMaxTot = -1e10, yMinTot = 1e10, yMaxTot = -1e10,
         yAbsMin = 1e10;
  for (int iHist = 0; iHist < int(histos.size()); ++iHist) {

    // Histogram information for plotting, especially x and y borders.
    string legendNow = (legends[iHist] != "void") ? legends[iHist]
      : histos[iHist].getTitle();
    stringstream nBin;
    nBin << histos[iHist].getBinNumber();
    double xMinNow = histos[iHist].getXMin();
    if (iHist == 0 || xMinNow < xMinTot) xMinTot = xMinNow;
    double xMaxNow = histos[iHist].getXMax();
    if (iHist == 0 || xMaxNow > xMaxTot) xMaxTot = xMaxNow;
    double yMinNow = histos[iHist].getYMin();
    if (iHist == 0 || yMinNow < yMinTot) yMinTot = yMinNow;
    double yMaxNow = histos[iHist].getYMax();
    if (iHist == 0 || yMaxNow > yMaxTot) yMaxTot = yMaxNow;
    double yAbsNow = histos[iHist].getYAbsMin();
    if (iHist == 0 || yAbsNow < yAbsMin) yAbsMin = yAbsNow;

    // Split plotting style and potential colour information.
    string styleNow = (styles[iHist] == "") ? "h" : styles[iHist];
    string style1 = styleNow;
    string style2 = "";
    if (styleNow.find(",") != string::npos) {
      int iComma = styleNow.find(",");
      style1 = (iComma > 0) ? styleNow.substr( 0, iComma) : "h";
      if (iComma + 1 < int(styleNow.length()))
        style2 = styleNow.substr( iComma + 1);
    }

    // Write histogram itself to a data file as two columns of (x,y) values.
    stringstream encode;
    encode << fileName << "-" << nTable + iHist << ".dat";
    histos[iHist].pyplotTable( encode.str(), (style1 == "h" || style1 == "e"),
      (style1 == "e"));

    // Write code to plot histogram.
    toPython << "plot = open('" << encode.str() << "')" << endl
             << "plot = [line.split() for line in plot]" << endl
             << "valx = [float(x[0]) for x in plot]" << endl
             << "valy = [float(x[1]) for x in plot]" << endl;
    if (style1 == "h" || style1 == "e")
      toPython  << "vale = [float(x[2]) for x in plot]"
                << endl << "plt.hist( valx, vale, weights = valy,"
                << " histtype='step',";
    else toPython << "plt.plot( valx, valy, '" << style1 << "',";
    if (style2 != "") toPython << " color='" << style2 << "',";
    toPython << " label=r\"" << legendNow << "\")" << endl;
    if (style1 == "e") {
      toPython << "erry = [float(x[3]) for x in plot]" << endl
               << "plt.errorbar( valx, valy, yerr=erry,"
               << " fmt='.', markersize=0, "
               << " color=plt.gca().patches[-1].get_edgecolor())" << endl;
    }
  }

  // Loop through the vector of already existing files, if any.
  if (histos.size() == 0) yAbsMin = 0.;
  for (int iFile = 0; iFile < int(files.size()); ++iFile) {
    string legendNow = (fileLegends[iFile] != "void") ? fileLegends[iFile]
      : files[iFile];

    // Split plotting style and potential colour information.
    string styleNow = (fileStyles[iFile] == "") ? "o" : fileStyles[iFile];
    string style1 = styleNow;
    string style2 = "";
    if (styleNow.find(",") != string::npos) {
      int iComma = styleNow.find(",");
      style1 = (iComma > 0) ? styleNow.substr( 0, iComma) : "o";
      if (iComma + 1 < int(styleNow.length()))
        style2 = styleNow.substr( iComma + 1);
    }

    // Find out whether and what kind of error bars should be plotted.
    int nxErr = 0;
    if (filexyerr[iFile].find("x") != string::npos) nxErr = 1;
    if (filexyerr[iFile].find("X") != string::npos) nxErr = 2;
    int nyErr = 0;
    if (filexyerr[iFile].find("y") != string::npos) nyErr = 1;
    if (filexyerr[iFile].find("Y") != string::npos) nyErr = 2;

    // Write code to plot existing file. Trivial but tedious error bars.
    toPython << "plot = open('" << files[iFile] << "')" << endl
             << "plot = [line.split() for line in plot]" << endl
             << "valx = [float(x[0]) for x in plot]" << endl
             << "valy = [float(x[1]) for x in plot]" << endl;
    if (style1 == "h") toPython
      << "vale = [float(x[2]) for x in plot]" << endl
      << "plt.hist( valx, vale, weights = valy,"  << " histtype='step',";
    else if (nxErr == 0 && nyErr == 0)
      toPython << "plt.plot( valx, valy, '" << style1 << "',";
    else {
      if (nxErr == 1) toPython
        << "errx = [float(x[2]) for x in plot]" << endl;
      else if (nxErr == 2) toPython
        << "errxlow = [float(x[2]) for x in plot]" << endl
        << "errxupp = [float(x[3]) for x in plot]" << endl
        << "errx = [errxlow, errxupp]" << endl;
      if (nyErr == 1) toPython
        << "erry = [float(x[" << nxErr + 2 << "]) for x in plot]" << endl;
      else if (nyErr == 2) toPython
        << "errylow = [float(x[" << nxErr + 2 << "]) for x in plot]" << endl
        << "erryupp = [float(x[" << nxErr + 3 << "]) for x in plot]" << endl
        << "erry = [errylow, erryupp]" << endl;
      toPython << "plt.errorbar( valx, valy,";
      if (nxErr > 0) toPython << " xerr=errx,";
      if (nyErr > 0) toPython << " yerr=erry,";
      toPython << " fmt='" << style1 << "',";
    }
    if (style2 != "") toPython << " color='" << style2 << "',";
    toPython << " label=r\"" << legendNow << "\", zorder=-1)" << endl;
  }

  // Set borders and write axes.
  toPython << "plt.xlim( " << scientific << setprecision(3)
           << (userBorders ? xMinUser : xMinTot) << ", "
           << (userBorders ? xMaxUser : xMaxTot) << ")" << endl;
  if ((histos.size() > 0 && !histos[0].getLinX()) || logX)
    toPython << "plt.xscale('log')" << endl;
  if (logY) {
    if (userBorders) toPython << "plt.ylim( " << scientific << setprecision(3)
      << yMinUser << ", " << yMaxUser << ")" << endl;
    if (useLegacy)
      toPython << "plt.yscale('symlog', linthreshy=";
    else
      toPython << "plt.yscale('symlog', linthresh=";
    toPython << scientific << setprecision(2)
             << (userBorders ? 0.5 * yMinUser : yAbsMin)
             << ")" << endl;
  } else {
    if (userBorders) {
      yMinTot = yMinUser;
      yMaxTot = yMaxUser;
    } else if (yMinTot < -1e-20 || yMinTot > 0.5 * yMaxTot) {
      double yMargin = 0.05 * (yMaxTot - yMinTot);
      yMinTot -= yMargin;
      yMaxTot += yMargin;
    } else {
      yMinTot  = 0.;
      yMaxTot *= 1.05;
    }
    toPython << "plt.ylim( " << scientific << setprecision(3) << yMinTot
             << ", " << yMaxTot << ")" << endl;
    toPython << "plt.ticklabel_format(axis='y', style='sci', "
             << "scilimits=(-2,3))" << endl;
  }

  // Write title and labels, and create plot.
  toPython << "plt.legend(frameon=False,loc='best')" << endl
           << "plt.title(r\"" << title << "\")" << endl
           << "plt.xlabel(r\"" << xLabel << "\")" << endl
           << "plt.ylabel(r\"" << yLabel << "\")" << endl
           << "pp.savefig(tmp" << nFrame << ",bbox_inches='tight')"
           << endl << "plt.clf()" << endl;

  // Update table counter. Done.
  nTable += histos.size();
}

//==========================================================================

} // end namespace Pythia8
