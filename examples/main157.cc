// main157.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim <marius.utheim@thep.lu.se>
//          Philip Ilten <philten@cern.ch>

// Keywords: low energy; resonances

// Illustration of how to create new resonance particles that can form
// during rescattering. In this case, the new resonance is the chi_1c(3872)
// tetraquark, which is produced from an initial D0 Dbar*0.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//--------------------------------------------------------------------------

int main() {

  // Initialize Pythia.
  Pythia pythia;

  // Define the chi_1c(3872)0 tetraquark.
  ParticleData& pd = pythia.particleData;
  pd.addParticle(9044111, "chi_1c(3872)0", "chi_1c(3872)bar0",
    2, 0, 0, 3.8720, 0.0012, 3.82, 3.92, 0., true);

  ParticleDataEntryPtr pde = pythia.particleData.findParticle(9044111);
  // Define decay channels. These can be used for resonance production.
  pde->addChannel(1, 0.3700, 3, 421,  -423); // D0 Dbar*0
  pde->addChannel(1, 0.0430, 3, 223,   443); // omega J/psi
  pde->addChannel(1, 0.0380, 3, 113,   443); // rho0 J/psi
  pde->addChannel(1, 0.0340, 3, 111, 20443); // pi0 chi_1c

  // Resonance production is only possible in hadron-hadron collisions.
  // These decay channels will not be used for resonance production, but
  // are relevant since they can change the particle composition.
  pde->addChannel(1, 0.0080, 3,  22,       443); // gamma J/psi
  pde->addChannel(1, 0.0450, 3,  22,    100443); // gamma psi(2S)
  pde->addChannel(1, 0.4900, 0, 421, -421, 111); // D0 Dbar0 pi0

  // Process specification. Resonance formation will be the dominant process.
  pythia.readString("LowEnergyQCD:all = on");
  pythia.readString("Beams:idA =  421");
  pythia.readString("Beams:idB = -423");
  pythia.readString("Beams:eCM = 3.9");

  // Number of events to generate/print.
  pythia.readString("Main:numberOfEvents  = 10");
  pythia.readString("Next:numberShowEvent = 10");

  // Initialize.
  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl;
    return -1;
  }

  // Generate events. They will be printed to cout.
  int nEvent = pythia.mode("Main:numberOfEvents");
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
    pythia.next();

  // Print statistics.
  pythia.stat();

  // Done.
  return 0;
}
