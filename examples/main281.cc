// main281.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: utility;

// This program illustrates how all settings and especially particle data
// can be written, and read back in after editing, either as plaintext or xml.

// For settings the xml version gives complete information, including min
// and max values, while min/max info is missing in the plaintext version.
// Since permanent editing of settings has to be done in xmldoc/*.xml files,
// the plaintext version may be more convenient for temporary changes.

// The particle data options, both the xml and plaintext ones, can also be
// used for temporary changes, and are then essentially equivalent.
// But here an edited xml file can be directly pasted into relevant parts
// of the share/Pythia8/xmldoc/ParticleData.xml file for a permanent update.

// More options are available; see the "Settings Scheme" and
// "Particle Data Scheme" in the html manual for further details.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Write out or read in. XML or plaintext (= free format).
  bool doSettings  = true;   // Can be combined freely with doParticles.
  bool doParticles = true;   // Can be combined freely with doSettings.
  bool doWrite     = true;   // If false then read instead.
  bool doXML       = true;   // If false then plaintext instead.
  bool checkRead   = true;   // Optional.

  // File names for settings and particle data.
  string settingsFile  = "settings";
  string particlesFile = "particles";

  // Create the Pythia instance, and read current settings and particle data.
  Pythia pythia;

  // Different ways to write settings.
  if (doWrite) {
    if (doSettings && doXML) {
      ofstream settingsStream( (settingsFile + ".xml").c_str() );
      pythia.settings.writeFileXML( settingsStream);
    } else if (doSettings)
      pythia.settings.writeFile( settingsFile + ".dat", true);

    // Different ways to write particle data.
    if (doParticles && doXML)
      pythia.particleData.listXML(particlesFile + ".xml");
    else if (doParticles)
      pythia.particleData.listFF(particlesFile + ".dat");

  // Different ways to read settings.
  } else {
    if (doSettings && doXML)
      pythia.settings.reInit( settingsFile + ".xml");
    else if (doSettings)
      pythia.readFile( settingsFile + ".dat");

    // Different ways to read particle data.
    if (doParticles && doXML)
      pythia.particleData.readXML(particlesFile + ".xml");
    else if (doParticles)
      pythia.particleData.readFF(particlesFile + ".dat");

    // Use normal listing formats to check if the read back in has worked.
    // These lists are formatted for human readability, and cannot be
    // read back in by any methods in Pythia.
    if (checkRead) {
      if (doSettings)  pythia.settings.listAll();
      if (doParticles) pythia.particleData.listAll();
    }
}

  // Done.
  return 0;
}
