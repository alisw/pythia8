// main158.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.utheim@thep.lu.se>
//          Philip Ilten <philten@cern.ch>

// Keywords: low energy; resonances

// Illustration of how to create new resonance particles that can form
// during rescattering. The resonances are exotic hadrons defined in
// main158.cmnd. The resulting invariant mass spectra are plotted, but the
// exotic hadrons are rare and statistics are low, even with 100k events.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Number of events. In 100k events, we expect a handful of T(2900)
  // tetraquarks, a few chi_1c(3872) tetraquarks, and possibly a few
  // P_c^+ pentaquarks.
  int nEvent = 100000;

  // Initialize Pythia.
  Pythia pythia;
  Event& event = pythia.event;

  // Read file for exotic hadron definitions.
  pythia.readFile("main158.cmnd");

  // Run configuration. Charmonium processes favour production of the
  // relevant exotic hadrons.
  pythia.readString("Charmonium:all = on");
  pythia.readString("Beams:eCM = 13000");

  // Enable rescattering. Retune pT0Ref to get correct charged multiplicity.
  pythia.readString("HadronLevel:rescatter = on");
  pythia.readString("MultipartonInteractions:pT0Ref = 2.345");

  // Initialize.
  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl;
    return -1;
  }

  // Exotic hadron mass spectra.
  Hist mX3872("chi_1c(3872)",  30, 3.82, 3.92);

  Hist mTcs0("T_0cs(2900)0",        30, 2.2, 3.6);
  Hist mTcs1("T_1cs(2900)0",        30, 2.2, 3.6);
  Hist mTcsbar0("T_csbar(2900)0",   30, 2.2, 3.6);
  Hist mTcsbarPP("T_csbar(2900)++", 30, 2.2, 3.6);

  Hist mPc4312("P_c(4312)+", 30, 4.2, 4.6);
  Hist mPc4440("P_c(4440)+", 30, 4.2, 4.6);
  Hist mPc4457("P_c(4457)+", 30, 4.2, 4.6);

  // Generate events.
  int nSuccess = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next())
      continue;
    nSuccess += 1;

    // Loop over event record and find tetraquarks.
    for (Particle& h : event) {
      if      (h.idAbs() == 9044111) mX3872.fill(h.m());
      else if (h.idAbs() == 9043211) mTcs0.fill(h.m());
      else if (h.idAbs() == 9043212) mTcs1.fill(h.m());
      else if (h.idAbs() == 9043213) mTcsbar0.fill(h.m());
      else if (h.idAbs() == 9043214) mTcsbarPP.fill(h.m());
      else if (h.idAbs() == 9422141) mPc4312.fill(h.m());
      else if (h.idAbs() == 9422142) mPc4440.fill(h.m());
      else if (h.idAbs() == 9422143) mPc4457.fill(h.m());
    }

  }

  // Normalize and calculate production cross sections.
  double sigGen = pythia.info.sigmaGen();
  mX3872.normalizeSpectrum(nSuccess); mX3872 *= sigGen;
  mTcs0.normalizeSpectrum(nSuccess); mTcs0 *= sigGen;
  mTcs1.normalizeSpectrum(nSuccess); mTcs1 *= sigGen;
  mTcsbar0.normalizeSpectrum(nSuccess); mTcsbar0 *= sigGen;
  mTcsbarPP.normalizeSpectrum(nSuccess); mTcsbarPP *= sigGen;
  mPc4312.normalizeSpectrum(nSuccess); mPc4312 *= sigGen;
  mPc4440.normalizeSpectrum(nSuccess); mPc4440 *= sigGen;
  mPc4457.normalizeSpectrum(nSuccess); mPc4457 *= sigGen;

  // Plot.
  string plotName = "main158";
  HistPlot plt(plotName);

  // Plot production cross sections for T(2900) tetraquarks.
  plt.frame(plotName, "$T(2900)$ tetraquark production cross sections",
    "$E_{CM}$ [GeV]", "$\\sigma$ [mb]");
  plt.add(Hist::plotFunc([&](double eCM) {
    return pythia.getSigmaPartial(-411, 321, eCM, 9043211);
  }, "", 200, 2.2, 3.6), "-", "$\\bar{D}^- K^+ \\to T_{cs0}(2900)^0$");
  plt.add(Hist::plotFunc([&](double eCM) {
    return pythia.getSigmaPartial(-411, 321, eCM, 9043212);
  }, "", 200, 2.2, 3.6), "-", "$\\bar{D}^- K^+ \\to T_{cs1}(2900)^0$");
  plt.add(Hist::plotFunc([&](double eCM) {
    return pythia.getSigmaPartial(431, -211, eCM, 9043213);
  }, "", 200, 2.2, 3.6), "-", "$D_s^+ \\pi^- \\to T_{c\\bar{s}}^a(2900)^0$");
  plt.add(Hist::plotFunc([&](double eCM) {
    return pythia.getSigmaPartial(431,  211, eCM, 9043214);
  }, "", 200, 2.2, 3.6), "-",
    "$D_s^+ \\pi^+ \\to T_{c\\bar{s}}^a(2900)^{++}$");
  plt.plot();

  // Plot invariant mass spectra.
  plt.frame(plotName, "$T(2900)$ tetraquark invariant mass spectra",
    "$m$ [GeV]", "$d\\sigma/dm$");
  plt.add(mTcs0, "h", "$T_{cs0}(2900)^0$");
  plt.add(mTcs1,  "h", "$T_{cs1}(2900)^0$");
  plt.add(mTcsbar0, "h", "$T_{c\\bar{s}}^a(2900)^0$");
  plt.add(mTcsbarPP, "h", "$T_{c\\bar{s}}^a(2900)^{++}$");
  plt.plot();

  plt.frame(plotName, "$\\chi_{c1}(3872)$ invariant mass spectra",
    "$m$ [GeV]", "$d\\sigma/dm$");
  plt.add(mX3872, "h", "$X(3872)$");
  plt.plot();

  plt.frame(plotName, "$P_c$ invariant mass spectra",
    "$m$ [GeV]", "$d\\sigma/dm$");
  plt.add(mPc4312, "h", "$P_c^+(4312)$");
  plt.add(mPc4440, "h", "$P_c^+(4440)$");
  plt.add(mPc4457, "h", "$P_c^+(4457)$");
  plt.plot();

  // Print statistics and done.
  pythia.stat();
  return 0;
}
