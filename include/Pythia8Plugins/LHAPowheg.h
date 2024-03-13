// LHAPowheg.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Philip Ilten, May 2015.

#ifndef Pythia8_LHAPowheg_H
#define Pythia8_LHAPowheg_H

#include "Pythia8/Plugins.h"
#include "Pythia8Plugins/LHAFortran.h"
#include <sys/stat.h>
#include <unistd.h>

namespace Pythia8 {

//==========================================================================

// Give access to the POWHEG commonblocks and subroutines.

extern "C" {

  // The random number common block.
  extern struct {
    int rnd_numseeds, rnd_initialseed, rnd_iwhichseed;
    char rnd_cwhichseed[4];
    int rnd_i1, rnd_i2;
  } pwhg_rnd_;

  // The RANMAR (R48 modification) common block.
  extern struct {
    double u[97];
    double c;
    int i97, j97;
  } r48st1_;

  // Initialize Powheg.
  void pwhginit_();

  // Reset the counters.
  void resetcnt_(const char *string, int length);

  // Generate an event.
  void pwhgevent_();

  // Access Powheg input data.
  double powheginput_(const char *string, int length);

}

//==========================================================================

// A plugin class to generate events with hard processes from
// POWHEGBOX matrix elements. See http://powhegbox.mib.infn.it/ for
// further details on POWHEGBOX.

// WARNING: If one wishes to use LHAPDF with both POWHEGBOX and
// Pythia, only LHAPDF 6 should be used. If not, and differing PDF
// sets are used between POWHEGBOX and Pythia, POWHEGBOX will not
// re-initialize the PDF set and consequently will use the PDF set
// last used by Pythia.

class LHAupPowheg : public LHAupFortran {

public:

  // Constructor.
  LHAupPowheg(Pythia* pythiaPtr, Settings* settingsPtr, Logger* loggerPtr);

  // Call pwhginit and fill the HEPRUP commonblock.
  bool fillHepRup();

  // Call pwhgevent and fill the HEEUP commonblock.
  bool fillHepEup();

private:

  // Read a POWHEG settings string.
  bool readString(string line);

  // Read a POWHEG settings file.
  bool readFile(string name);

  // Write out the input for POWHEG.
  bool init();

  // The external random number generator.
  Rndm* rndmPtr{};

  // Logger.
  Logger* loggerPtr{};

  // The POWHEG run directory and PDF file (if not LHAPDF).
  string dir, pdf;

  // The map of POWHEG settings.
  map<string, string> settings;

  // The current working directory.
  char cwd[FILENAME_MAX];

};

//--------------------------------------------------------------------------

// Constructor.

LHAupPowheg::LHAupPowheg(Pythia *pythiaPtr, Settings *settingsPtr,
  Logger* loggerPtrIn) : loggerPtr(loggerPtrIn), dir("powhegrun"), pdf("") {

  // Load the settings.
  if (settingsPtr != nullptr) {
    dir = settingsPtr->word("POWHEG:dir");
    pdf = settingsPtr->word("POWHEG:pdf");

    // Read the command files and strings.
    for (string cmndFile : settingsPtr->wvec("POWHEG:cmndFiles"))
      readFile(cmndFile);
    for (string cmnd : settingsPtr->wvec("POWHEG:cmnds"))
      readString(cmnd);
  }

  // Set up the random number generator.
  if (pythiaPtr != nullptr) {
    if (settingsPtr->isFlag("POWHEG:pythiaRandom")) rndmPtr = &pythiaPtr->rndm;
  }
  mkdir(dir.c_str(), 0777);
  init();

}

//--------------------------------------------------------------------------

// Call pwhginit and fill the HEPRUP commonblock.

bool LHAupPowheg::fillHepRup() {

  // Set multiple random seeds to none.
  getcwd(cwd, sizeof(cwd));
  chdir(dir.c_str());
  strcpy(pwhg_rnd_.rnd_cwhichseed, "none");

  // Initialize Powheg.
  pwhginit_();

  // Reset all the counters.
  resetcnt_("upper bound failure in inclusive cross section", 46);
  resetcnt_("vetoed calls in inclusive cross section", 39);
  resetcnt_("upper bound failures in generation of radiation", 47);
  resetcnt_("vetoed radiation", 16);
  chdir(cwd);
  return fillHepEup();

}

//--------------------------------------------------------------------------

// Set the random numbers, call pwhgevent, and fill the HEPEUP commonblock.

bool LHAupPowheg::fillHepEup() {

  // Change directory.
  getcwd(cwd, sizeof(cwd));
  chdir(dir.c_str());

  // Reset the random block if requested.
  if (rndmPtr != nullptr) {
    r48st1_.i97 = 97;
    r48st1_.j97 = 33;
    r48st1_.c = rndmPtr->flat();
    for (int i = 0; i < 97; ++i) r48st1_.u[i] = rndmPtr->flat();
  }

  // Generate the event.
  pwhgevent_();
  chdir(cwd);
  return true;

}

//--------------------------------------------------------------------------

// Read a POWHEG settings string. If a setting is repeated a warning
// is printed but the most recent setting is used.

bool LHAupPowheg::readString(string line) {

  // Copy string without initial and trailing blanks.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;
  int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");
  int lastChar  = line.find_last_not_of(" \n\t\v\b\r\f\a");
  line = line.substr(firstChar, lastChar + 1 - firstChar);

  // Find the key.
  firstChar = line.find_first_of("  \t\f\v\n\r");
  string key = toLower( line.substr(0, firstChar), false);

  // Add the setting.
  if (key.size() > 0
    && key.find_first_of("abcdedfghijklmnopqrtsuvwxyz") == 0) {
    map<string, string>::iterator setting = settings.find(key);
    if (setting != settings.end()) {
      if (loggerPtr) loggerPtr->WARNING_MSG(
        "replacing previous POWHEG setting for " + key);
      setting->second = line;
    } else settings[key] = line;
  }
  return true;

}

//--------------------------------------------------------------------------

// Read a POWHEG settings file.

bool LHAupPowheg::readFile(string name) {

  fstream config(name.c_str(), ios::in); string line;
  while (getline(config, line, '\n')) readString(line);
  config.close();
  return true;

}

//--------------------------------------------------------------------------

// Write the input for POWHEG.

bool LHAupPowheg::init() {

  // Copy over the PDF file if needed.
  if (pdf != "") {
    fstream pdfin(pdf.c_str(), ios::in | ios::binary);
    fstream pdfout((dir + "/" + pdf.substr(0, pdf.find_last_of("/"))).c_str(),
      ios::out | ios::binary);
    pdfout << pdfin.rdbuf();
    pdfin.close();
    pdfout.close();
  }

  // Copy the settings to the configuration file.
  fstream config((dir + "/" + "powheg.input").c_str(), ios::out);
  for (map<string, string>::iterator setting = settings.begin();
    setting != settings.end(); ++setting) config << setting->second << "\n";
  config.close();
  return true;

}

//--------------------------------------------------------------------------

// Register POWHEG settings.

void powhegSettings(Settings *settingsPtr) {
  settingsPtr->addWord("POWHEG:dir", "powhegrun");
  settingsPtr->addWord("POWHEG:pdf", "");
  settingsPtr->addWVec("POWHEG:cmnds", {});
  settingsPtr->addWVec("POWHEG:cmndFiles", {});
  settingsPtr->addFlag("POWHEG:pythiaRandom", true);
}

//--------------------------------------------------------------------------

// Declare the plugin.

PYTHIA8_PLUGIN_CLASS(LHAup, LHAupPowheg, true, true, true)
PYTHIA8_PLUGIN_SETTINGS(powhegSettings)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LHAPowheg_H
