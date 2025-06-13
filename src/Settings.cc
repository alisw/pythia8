// Settings.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Settings class.

#include "Pythia8/Settings.h"
#include "Pythia8/Plugins.h"

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {

//==========================================================================

// Settings class.
// This class contains flags, modes, parms and words used in generation.

//--------------------------------------------------------------------------

// Read in database from specific file.

bool Settings::init(string startFile, bool append) {

  // Don't initialize if it has already been done and not in append mode.
  if (isInit && !append) return true;
  int nError = 0;

  // Reset readString history.
  if (!append) {
    readStringHistory.resize(0);
    readStringSubrun.clear();
  }

  // List of files to be checked. Start with input file.
  vector<string> files;
  files.push_back(startFile);

  // If nontrivial startfile path, then use that for other files as well.
  string xmlNow = "";
  if (startFile.rfind("/") != string::npos)
    xmlNow = startFile.substr(0, startFile.rfind("/") + 1);

  // Only update the XML path if reading in the main XML.
  if (!append) xmlPath = xmlNow;

  // Loop over files. Open them for read.
  for (int i = 0; i < int(files.size()); ++i) {
    const char* cstring = files[i].c_str();
    ifstream is(cstring);

    // Check that instream is OK.
    if (!is.good()) {
      loggerPtr->ERROR_MSG("settings file " + files[i] + " not found");
      return false;
    }

    // Read in one line at a time.
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line, to interpret it as tag. Remove "more".
      istringstream getfirst(line);
      string tag;
      getfirst >> tag;
      if (tag.find("more") != string::npos) tag.erase( tag.find("more"), 4);

      // Skip ahead if not interesting.
      if (tag != "<flag" && tag != "<flagfix" && tag != "<mode"
         && tag != "<modeopen" && tag != "<modepick" && tag != "<modefix"
         && tag != "<parm" && tag != "<parmfix" && tag != "<word"
         && tag != "<wordfix" && tag != "<fvec" && tag != "<fvecfix"
         && tag != "<mvec" && tag != "<mvecfix"
         && tag != "<pvec" && tag != "<pvecfix"
         && tag != "<wvec" && tag != "<wvecfix" && tag != "<aidx") continue;

      // Read and append continuation line(s) if line does not contain >.
      while (line.find(">") == string::npos) {
        string addLine;
        getline(is, addLine);
        line += " " + addLine;
      }

      // Remove extra blanks before an = sign.
      while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);

      // Add file also to be read. Avoid duplication.
      if (tag == "<aidx") {
        string name = attributeValue( line, "href");
        if (name == "") {
          cout << " PYTHIA Error: failed to find name attribute in line "
               << line << endl;
          ++nError;
          continue;
        }
        string fullName = xmlNow + name + ".xml";
        bool inBase = false;
        for (int ifile = 0; ifile < int(files.size()); ++ifile)
          if (fullName == files[ifile]) inBase = true;
        if (!inBase) files.push_back(fullName);
        continue;
      }

      // Find name attribute.
      string name = attributeValue( line, "name=");
      if (name == "") {
        cout << " PYTHIA Error: failed to find name attribute in line "
             << line << endl;
        ++nError;
        continue;
      }

      // Check that default value attribute present, and whether max and min.
      if (line.find("default=") == string::npos) {
        cout << " PYTHIA Error: failed to find default value token in line "
             << line << endl;
        ++nError;
        continue;
      }
      bool hasMin = (line.find("min=") != string::npos);
      bool hasMax = (line.find("max=") != string::npos);

      // Check for occurence of a bool and add to flag map.
      if (tag == "<flag" || tag == "<flagfix") {
        bool value = boolAttributeValue( line, "default=");
        addFlag( name, value);

      // Check for occurence of an int and add to mode map.
      } else if (tag == "<mode" || tag == "<modeopen"
        || tag == "<modepick" || tag == "<modefix") {
        int value    = intAttributeValue( line, "default=");
        int minVal   = intAttributeValue( line, "min=");
        int maxVal   = intAttributeValue( line, "max=");
        // Enforce check that only allowed options are accepted.
        bool optOnly = false;
        if (tag == "<modepick" && hasMin && hasMax) optOnly = true;
        if (tag == "<modefix") {
          hasMin  = true;
          hasMax  = true;
          minVal  = value;
          maxVal  = value;
          optOnly = true;
        }
        addMode( name, value, hasMin, hasMax, minVal, maxVal, optOnly);

      // Check for occurence of a double and add to parm map.
      } else if (tag == "<parm" || tag == "<parmfix") {
        double value  = doubleAttributeValue( line, "default=");
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addParm( name, value, hasMin, hasMax, minVal, maxVal);

      // Check for occurence of a string and add to word map.
      } else if (tag == "<word" || tag == "<wordfix") {
        string value = attributeValue( line, "default=");
        addWord( name, value);

      // Check for occurence of a bool vector and add to fvec map.
      } else if (tag == "<fvec" || tag == "<fvecfix") {
        vector<bool> value = boolVectorAttributeValue( line, "default=");
        addFVec( name, value);

      // Check for occurence of an int vector and add to mvec map.
      } else if (tag == "<mvec" || tag == "<mvecfix") {
        vector<int> value = intVectorAttributeValue( line, "default=");
        int minVal = intAttributeValue( line, "min=");
        int maxVal = intAttributeValue( line, "max=");
        addMVec( name, value, hasMin, hasMax, minVal, maxVal);

      // Check for occurence of a double vector and add to pvec map.
      } else if (tag == "<pvec" || tag == "<pvecfix") {
        vector<double> value = doubleVectorAttributeValue( line, "default=");
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addPVec( name, value, hasMin, hasMax, minVal, maxVal);

      // Check for occurence of a string vector and add to wvec map.
      } else if (tag == "<wvec" || tag == "<wvecfix") {
        vector<string> value = stringVectorAttributeValue( line, "default=");
        addWVec( name, value);
      }

    // End of loop over lines in input file and loop over files.
    };
  };

  // Set up default e+e- and pp tunes, if positive.
  // Only initialize the tunes if the main XML init, i.e. append is false.
  if (!append) {
    int eeTune = mode("Tune:ee");
    if (eeTune > 0) initTuneEE(eeTune);
    int ppTune = mode("Tune:pp");
    if (ppTune > 0) initTunePP(ppTune);
    int vinciaTune = mode("Vincia:Tune");
    if (vinciaTune > -1) initTuneVincia(vinciaTune);
  }

  // Done.
  if (nError > 0) return false;
  isInit = true;
  return true;

}


//--------------------------------------------------------------------------

// Read in database from specific stream.

bool Settings::init(istream& is, bool append) {

  // Don't initialize if it has already been done and not in append mode.
  if (isInit && !append) return true;
  int nError = 0;

  // Check that instream is OK.
  if (!is.good()) {
    cout << "\n PYTHIA Error: settings stream not found " << endl;
    return false;
  }

  // Reset readString history.
  if (!append) {
    readStringHistory.resize(0);
    readStringSubrun.clear();
  }

  // Read in one line at a time.
  string line;
  while ( getline(is, line) ) {

    // Get first word of a line, to interpret it as tag. Remove "more".
    istringstream getfirst(line);
    string tag;
    getfirst >> tag;
    if (tag.find("more") != string::npos) tag.erase( tag.find("more"), 4);

    // Skip ahead if not interesting. Only look for new files in startfile.
    if (tag != "<flag" && tag != "<flagfix" && tag != "<mode"
       && tag != "<modeopen" && tag != "<modepick" && tag != "<modefix"
       && tag != "<parm" && tag != "<parmfix" && tag != "<word"
       && tag != "<wordfix" && tag != "<fvec" && tag != "<fvecfix"
       && tag != "<mvec" && tag != "<mvecfix"
       && tag != "<pvec" && tag != "<pvecfix"
       && tag != "<wvec" && tag != "<wvecfix" && tag != "<aidx") continue;

    // Read and append continuation line(s) if line does not contain >.
    while (line.find(">") == string::npos) {
      string addLine;
      getline(is, addLine);
      line += " " + addLine;
    }

    // Remove extra blanks before an = sign.
    while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);

    // Find name attribute.
    string name = attributeValue( line, "name=");
    if (name == "") {
      cout << " PYTHIA Error: failed to find name attribute in line "
           << line << endl;
      ++nError;
      continue;
    }

    // Check that default value attribute present, and whether max and min.
    if (line.find("default=") == string::npos) {
      cout << " PYTHIA Error: failed to find default value token in line "
           << line << endl;
      ++nError;
      continue;
    }
    bool hasMin = (line.find("min=") != string::npos);
    bool hasMax = (line.find("max=") != string::npos);

    // Check for occurence of a bool and add to flag map.
    if (tag == "<flag" || tag == "<flagfix") {
      bool value = boolAttributeValue( line, "default=");
      addFlag( name, value);

    // Check for occurence of an int and add to mode map.
    } else if (tag == "<mode" || tag == "<modeopen"
      || tag == "<modepick" || tag == "<modefix") {
      int value    = intAttributeValue( line, "default=");
      int minVal   = intAttributeValue( line, "min=");
      int maxVal   = intAttributeValue( line, "max=");

      // Enforce check that only allowed options are accepted.
      bool optOnly = false;
      if (tag == "<modepick" && hasMin && hasMax) optOnly = true;
      if (tag == "<modefix") {
        hasMin  = true;
        hasMax  = true;
        minVal  = value;
        maxVal  = value;
        optOnly = true;
      }
      addMode( name, value, hasMin, hasMax, minVal, maxVal, optOnly);

    // Check for occurence of a double and add to parm map.
    } else if (tag == "<parm" || tag == "<parmfix") {
      double value  = doubleAttributeValue( line, "default=");
      double minVal = doubleAttributeValue( line, "min=");
      double maxVal = doubleAttributeValue( line, "max=");
      addParm( name, value, hasMin, hasMax, minVal, maxVal);

    // Check for occurence of a string and add to word map.
    } else if (tag == "<word" || tag == "<wordfix") {
      string value = attributeValue( line, "default=");
      addWord( name, value);

    // Check for occurence of a bool vector and add to fvec map.
    } else if (tag == "<fvec" || tag == "<fvecfix") {
      vector<bool> value = boolVectorAttributeValue( line, "default=");
      addFVec( name, value);

    // Check for occurence of an int vector and add to mvec map.
    } else if (tag == "<mvec" || tag == "<mvecfix") {
      vector<int> value = intVectorAttributeValue( line, "default=");
      int minVal = intAttributeValue( line, "min=");
      int maxVal = intAttributeValue( line, "max=");
      addMVec( name, value, hasMin, hasMax, minVal, maxVal);

    // Check for occurence of a double vector and add to pvec map.
    } else if (tag == "<pvec" || tag == "<pvecfix") {
      vector<double> value = doubleVectorAttributeValue( line, "default=");
      double minVal = doubleAttributeValue( line, "min=");
      double maxVal = doubleAttributeValue( line, "max=");
      addPVec( name, value, hasMin, hasMax, minVal, maxVal);

    // Check for occurence of a string vector and add to word map.
    } else if (tag == "<wvec" || tag == "<wvecfix") {
      vector<string> value = stringVectorAttributeValue( line, "default=");
      addWVec( name, value);
    }

  // End of loop over lines in input file and loop over files.
  };

  // Set up default e+e- and pp tunes, if positive.
  int eeTune = mode("Tune:ee");
  if (eeTune > 0) initTuneEE(eeTune);
  int ppTune = mode("Tune:pp");
  if (ppTune > 0) initTunePP(ppTune);
  int vinciaTune = mode("Vincia:Tune");
  if (vinciaTune > -1) initTuneVincia(vinciaTune);

  // Done.
  if (nError > 0) return false;
  isInit = true;
  return true;

}

//--------------------------------------------------------------------------

// Overwrite existing database by reading from specific file.

bool Settings::reInit(string startFile) {

  // Reset maps to empty.
  flags.clear();
  modes.clear();
  parms.clear();
  words.clear();
  fvecs.clear();
  mvecs.clear();
  pvecs.clear();
  wvecs.clear();

  // Then let normal init do the rest.
  isInit = false;
  return init(startFile, false);

}

//--------------------------------------------------------------------------

// Read in updates from a character string, like a line of a file.
// Is used by readString (and readFile) in Pythia.

bool Settings::readString(string line, bool warn, int subrun) {

  // If empty line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // If unfinished line then add new to existing, else use input line as is.
  string lineNow = (lineSaved) ? savedLine + line : line;
  lineSaved = false;
  int firstChar = lineNow.find_first_not_of(" \n\t\v\b\r\f\a");

  // Send on particle data to the ParticleData database if connected.
  if (isdigit(lineNow[firstChar])) {
    if (pdPtr != nullptr) {
      bool passed = pdPtr->readString(lineNow, warn);
      if (passed && pdbPtr != nullptr) (*pdbPtr) << lineNow << endl;
      return passed;
    } else return false;
  }

  // If first character is not a letter, then taken to be a comment line.
  if (!isalpha(lineNow[firstChar])) return true;

  // Allow the += notation to add to a settings vector.
  bool append(false);

  // Replace an equal sign by a blank to make parsing simpler, except after {.
  size_t iBrace = (lineNow.find_first_of("{") == string::npos) ? lineNow.size()
    : lineNow.find_first_of("{");
  while (lineNow.find("=") != string::npos
    && lineNow.find_first_of("=") < iBrace) {
    int firstEqual = lineNow.find_first_of("=");
    lineNow.replace(firstEqual, 1, " ");
    if ( lineNow[firstEqual-1] == '+' ) {
      lineNow[firstEqual-1] = ' ';
      append = true;
    }
  }

  // Get first word of a line.
  istringstream splitLine(lineNow);
  string name;
  splitLine >> name;
  name = toLower(name);

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");
  }

  // Check the subrun.
  if (name != "main:subrun"
    && subrunNow != subrun && subrunNow != SUBRUNDEFAULT) return true;

  // Check whether this is in the database.
  int inDataBase = 0;
  if      (isFlag(name))      inDataBase = 1;
  else if (isMode(name))      inDataBase = 2;
  else if (isParm(name))      inDataBase = 3;
  else if (isWord(name))      inDataBase = 4;
  else if (isFVec(name))      inDataBase = 5;
  else if (isMVec(name))      inDataBase = 6;
  else if (isPVec(name))      inDataBase = 7;
  else if (isWVec(name))      inDataBase = 8;
  else if (name == "include") inDataBase = 9;

  // Warn and done if not in database.
  if (inDataBase == 0) {
    if (warn) cout << "\n PYTHIA Error: input string not found in settings"
      << " databases::\n   " << line << endl;
    readingFailedSave = true;
    return false;
  }
  if (append && inDataBase < 5  ) {
    if (warn) cout << "\n PYTHIA Error: the += notation is only"
                   <<" valid for vector settings:\n   " << line << endl;
    readingFailedSave = true;
    return false;
  }

  // Find value. Warn if none found, except that a word (string) can be empty.
  string valueString;
  splitLine >> valueString;
  if (!splitLine && inDataBase == 4) valueString = "";
  else if (!splitLine) {
    if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
      << " not meaningful:\n   " << line << endl;
    readingFailedSave = true;
    return false;
  }

  // If value is a ? then echo the current value.
  if (valueString == "?") {
    cout << output(name);
    return true;
  }

  // If value is DEFAULT, then reset to the default value.
  if (valueString == "DEFAULT") {
    if      (inDataBase == 1) resetFlag(name);
    else if (inDataBase == 2) resetMode(name);
    else if (inDataBase == 3) resetParm(name);
    else if (inDataBase == 4) resetWord(name);
    else if (inDataBase == 5) resetFVec(name);
    else if (inDataBase == 6) resetMVec(name);
    else if (inDataBase == 7) resetPVec(name);
    else if (inDataBase == 8) resetWVec(name);
    return true;
  }

  // Check for FORCE= statements (to ignore min/max values)
  bool force = false;
  if (valueString.find("force") != string::npos) {
    force = true;
    // Read value from next word
    splitLine >> valueString;
    if (!splitLine) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
         << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
  }

  // If string begins with { then find matching } and extract contents.
  if (valueString[0] == '{') {
    size_t openBrace  = lineNow.find_first_of("{");
    size_t closeBrace = lineNow.find_first_of("}");
    // If not } on same line then must append next line and try again.
    if (closeBrace == string::npos) {
      lineSaved = true;
      savedLine = lineNow;
      return true;
    }
    valueString = lineNow.substr(openBrace, closeBrace - openBrace + 1);
  }

  // Update flag map; allow many ways to say yes.
  if (inDataBase == 1) {
    bool value = boolString(valueString);
    flag(name, value, force);

  // Update mode map.
  } else if (inDataBase == 2) {
    istringstream modeData(valueString);
    int value;
    modeData >> value;
    if (!modeData) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (name == "main:subrun") {
      subrunNow = value;
      if (subrun != subrunNow) return true;
    }
    if (!mode(name, value, force)) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " is out of range:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }

  // Update parm map.
  } else if (inDataBase == 3) {
    istringstream parmData(valueString);
    double value;
    parmData >> value;
    if (!parmData) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (!parm(name, value, force)) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " is out of range:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }

  // Update word map.
  } else if (inDataBase == 4)  {
    word(name, valueString, force);

  // Update fvec map.
  } else if (inDataBase == 5) {
    istringstream fvecData(valueString);
    vector<bool> value(boolVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!fvecData) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (append) {
      vector<bool> old = fvec(name);
      old.insert(old.end(), value.begin(), value.end());
      value = old;
    }
    fvec(name, value, force);

  // Update mvec map.
  } else if (inDataBase == 6) {
    istringstream mvecData(valueString);
    vector<int> value(intVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!mvecData) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (append) {
      vector<int> old = mvec(name);
      old.insert(old.end(), value.begin(), value.end());
      value = old;
    }
    if (!mvec(name, value, force)) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " is out of range:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }

  // Update pvec map.
  } else if (inDataBase == 7) {
    istringstream pvecData(valueString);
    vector<double> value(doubleVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!pvecData) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (append) {
      vector<double> old = pvec(name);
      old.insert(old.end(), value.begin(), value.end());
      value = old;
    }
    if (!pvec(name, value, force)) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " is out of range:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }

  // Update wvec map.
  } else if (inDataBase == 8) {
    istringstream wvecData(valueString);
    vector<string> value(stringVectorAttributeValue(
      "value=\"" + valueString + "\"", "value="));
    if (!wvecData) {
      if (warn) cout << "\n PYTHIA Error: variable recognized, but its value"
        << " not meaningful:\n   " << line << endl;
      readingFailedSave = true;
      return false;
    }
    if (append) {
      vector<string> old = wvec(name);
      old.insert(old.end(), value.begin(), value.end());
      value = old;
    }
    wvec(name, value, force);

  // Include files directive.
  } else if (inDataBase == 9) {

    // Try normal path first.
    ifstream isUser(valueString.c_str());
    if (!isUser.good()) {

      // Split the paths from PYTHIA8CMND.
      vector<string> paths;
      size_t pos(0);
      const char* envChar = getenv("PYTHIA8CMND");
      string envPath = envChar ? envChar : "";
      while (envPath != "" && pos != string::npos) {
        pos = envPath.find(":");
        paths.push_back(envPath.substr(0, pos));
        envPath = envPath.substr(pos + 1);
      }

      // Add the Pythia share directory.
      paths.push_back(xmlPath + "../");

      // Try the different paths.
      for (string path : paths) {
        ifstream isPath((path + "/" + valueString).c_str());
        if (isPath.good()) return readFile(isPath, warn, subrun);
      }
      loggerPtr->ERROR_MSG("did not find file", valueString);
      loggerPtr->ERROR_MSG("searched along the following paths:");
      for (string path : paths) loggerPtr->ERROR_MSG(path);
      readingFailedSave = true;
      return false;
    } else return readFile(isUser, warn, subrun);
  }

  // Store history of valid readString statements.
  readStringHistory.push_back(lineNow);
  readStringSubrun[max(-1, mode("Main:subrun"))].push_back(lineNow);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Read in updates from a user-defined file.

bool Settings::readFile(string fileName, bool warn, int subrun) {

  // Open file for reading.
  ifstream is(fileName.c_str());
  if (!is.good()) {
    loggerPtr->ERROR_MSG("did not find file", fileName);
    return false;
  }

  // Hand over real work to next method.
  return readFile(is, warn, subrun);

}

//--------------------------------------------------------------------------

// Read in updates from a user-defined stream (or file).

bool Settings::readFile(istream& is, bool warn, int subrun) {

  // Read in one line at a time.
  string line;
  bool isCommented = false;
  bool accepted = true;
  subrunNow = SUBRUNDEFAULT;
  while (getline(is, line)) {

    // Check whether entering, leaving or inside commented-commands section.
    int    pos = line.find_first_not_of(" \n\t\v\b\r\f\a");
    string sub = line.length() - pos > 2 ? line.substr(pos, 2) : "";
    if      (sub == "/*") isCommented = true;
    else if (sub == "*/") isCommented = false;
    else if (!isCommented && !readString(line, warn, subrun)) accepted = false;
  }

  // Reached end of input file.
  return accepted;

}

//--------------------------------------------------------------------------

// Load a plugin library, and register any settings.

bool Settings::registerPluginLibrary(string libName, string startFile) {

  // Check if already loaded.
  if (pluginLibraries.find(libName) != pluginLibraries.end()) return false;
  pluginLibraries.insert(libName);

  // Load the plugin library.
  shared_ptr<void> libPtr = dlopen_plugin(libName, loggerPtr);
  if (libPtr == nullptr) return false;

  // Check if XML index has been specified.
  if (startFile == "") {
    auto xmlIndex = dlsym_plugin<const char*()>(libPtr, "RETURN_XML");
    if (dlerror() == nullptr) startFile = xmlIndex();
  }

  // Find the path to the XML, first PYTHIA8CONTRIB, then Pythia XML path.
  const char* envPath = getenv("PYTHIA8CONTRIB");
  string xmlNow = envPath ? envPath : "";
  if (xmlNow.length() && xmlNow[xmlNow.length() - 1] != '/') xmlNow += "/";
  ifstream xmlFile((xmlNow + startFile).c_str());
  if (!xmlFile.good()) {
    xmlFile.close();
    xmlNow = xmlPath + "../../";
    xmlFile.open((xmlNow + startFile).c_str());
    if (!xmlFile.good()) xmlNow = "";
  }
  xmlFile.close();

  // Load the XML files, if specified.
  if (startFile != "") init(xmlNow + startFile, true);

  // Load the settings registration symbol.
  auto registerSettings =
    dlsym_plugin<void(Settings*)>(libPtr, "REGISTER_SETTINGS");
  if (dlerror() != nullptr) return false;
  registerSettings(this);
  return true;

}

//--------------------------------------------------------------------------

// Write updates or everything to user-defined file.

bool Settings::writeFile(string toFile, bool writeAll) {

  // Open file for writing.
  const char* cstring = toFile.c_str();
  ofstream os(cstring);
  if (!os) {
    loggerPtr->ERROR_MSG("could not open file", toFile);
    return false;
  }

  // Hand over real work to next method.
  return writeFile( os, writeAll);

}

//--------------------------------------------------------------------------

// Write updates or everything to user-defined stream (or file).

bool Settings::writeFile(ostream& os, bool writeAll) {

  // Write simple header as comment.
  if (writeAll) os << "! List of all current PYTHIA ";
  else          os << "! List of all modified PYTHIA ";
  os << fixed << setprecision(3) << parm("Pythia:versionNumber")
     << " settings.\n";

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();
  map<string, FVec>::iterator fvecEntry = fvecs.begin();
  map<string, MVec>::iterator mvecEntry = mvecs.begin();
  map<string, PVec>::iterator pvecEntry = pvecs.begin();
  map<string, WVec>::iterator wvecEntry = wvecs.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end()
      || parmEntry != parms.end() || wordEntry != words.end()
      || fvecEntry != fvecs.end() || mvecEntry != mvecs.end()
      || pvecEntry != pvecs.end() || wvecEntry != wvecs.end() ) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end()
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first )
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || flagEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || flagEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || flagEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || flagEntry->first < wvecEntry->first )
      ) {
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << flagEntry->second.name << " = " + toString(valNow) + "\n";
      ++flagEntry;

    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end()
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || modeEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || modeEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || modeEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || modeEntry->first < wvecEntry->first )
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << modeEntry->second.name << " = " + toString(valNow) + "\n";
      ++modeEntry;

    // Else check if parm is next, and if so print it; fixed or scientific.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || parmEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || parmEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || parmEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || parmEntry->first < wvecEntry->first )
      ) {
      double valNow = parmEntry->second.valNow;
      double valDefault = parmEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os  << parmEntry->second.name << " = " + toString(valNow) + "\n";
      ++parmEntry;

    // Else check if word is next, and if so print it.
    } else  if ( wordEntry != words.end()
      && ( fvecEntry == fvecs.end() || wordEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || wordEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || wordEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || wordEntry->first < wvecEntry->first )
      ) {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << wordEntry->second.name << " = " + valNow + "\n";
      ++wordEntry;

    // Else check if fvec is next, and if so print it.
    } else if ( fvecEntry != fvecs.end()
      && ( mvecEntry == mvecs.end() || fvecEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || fvecEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || fvecEntry->first < wvecEntry->first )
      ) {
      vector<bool> valNow = fvecEntry->second.valNow;
      vector<bool> valDefault = fvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << fvecEntry->second.name << " = {";
        if (valNow.size() > 0) {
          for (vector<bool>::iterator val = valNow.begin();
               val != --valNow.end(); ++val) os << toString(*val) + ",";
          os << toString(*(--valNow.end())) + "}\n";
        } else os << "}\n";
      }
      ++fvecEntry;

    // Else check if mvec is next, and if so print it.
    } else if ( mvecEntry != mvecs.end()
      && ( pvecEntry == pvecs.end() || mvecEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || mvecEntry->first < wvecEntry->first )
      ) {
      vector<int> valNow = mvecEntry->second.valNow;
      vector<int> valDefault = mvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << mvecEntry->second.name << " = {";
        if (valNow.size() > 0) {
          for (vector<int>::iterator val = valNow.begin();
               val != --valNow.end(); ++val) os << toString(*val) + ",";
          os << toString(*(--valNow.end())) + "}\n";
        } else os << "}\n";
      }
      ++mvecEntry;

    // Else check if pvec is next; print fixed or scientific.
    } else if ( pvecEntry != pvecs.end()
      && ( wvecEntry == wvecs.end() || pvecEntry->first < wvecEntry->first )
      ) {
      vector<double> valNow = pvecEntry->second.valNow;
      vector<double> valDefault = pvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        os  << pvecEntry->second.name << " = {";
        if (valNow.size() > 0) {
          for (vector<double>::iterator val = valNow.begin();
               val != --valNow.end(); ++val) os << toString(*val) + ",";
          os << toString(*(--valNow.end())) + "}\n";
        } else os << "}\n";
      }
      ++pvecEntry;

    // Else print wvec.
    } else {
      vector<string> valNow = wvecEntry->second.valNow;
      vector<string> valDefault = wvecEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) {
        if (valNow.size() > 0) {
          os  << wvecEntry->second.name << " = {";
          for (vector<string>::iterator val = valNow.begin();
               val != --valNow.end(); ++val) {
            if ( val != valNow.begin() ) os << "\n      ";
            os << trimString(*val) + ",";
          }
          if ( valNow.size() > 1 )  os << "\n      ";
          os << *(--valNow.end()) + "}\n";
        } else os << "}\n";
      }
      ++wvecEntry;
    }
  } ;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Write updates or everything to user-defined stream (or file).

bool Settings::writeFileXML(ostream& os) {

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();
  map<string, FVec>::iterator fvecEntry = fvecs.begin();
  map<string, MVec>::iterator mvecEntry = mvecs.begin();
  map<string, PVec>::iterator pvecEntry = pvecs.begin();
  map<string, WVec>::iterator wvecEntry = wvecs.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end()
      || parmEntry != parms.end() || wordEntry != words.end()
      || fvecEntry != fvecs.end() || mvecEntry != mvecs.end()
      || pvecEntry != pvecs.end() || wvecEntry != wvecs.end() ) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end()
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first )
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || flagEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || flagEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || flagEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || flagEntry->first < wvecEntry->first )
       ) {
      string state[2] = {"off", "on"};
      bool valDefault = flagEntry->second.valDefault;
      os << "<flag name=\"" << flagEntry->second.name << "\" default=\""
         << state[valDefault] << "\"></flag>" << endl;
      ++flagEntry;

    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end()
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || modeEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || modeEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || modeEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || modeEntry->first < wvecEntry->first )
      ) {
      int valDefault = modeEntry->second.valDefault;
      os << "<mode name=\"" << modeEntry->second.name << "\" default=\""
         << valDefault << "\"";
      if (modeEntry->second.hasMin ) os << " min=\""
        << modeEntry->second.valMin << "\"";
      if (modeEntry->second.hasMax ) os << " max=\""
        << modeEntry->second.valMax << "\"";
      os << "></mode>" << endl;
      ++modeEntry;

    // Else check if parm is next, and if so print it; fixed or scientific.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || parmEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || parmEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || parmEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || parmEntry->first < wvecEntry->first )
      ) {
      double valDefault = parmEntry->second.valDefault;
      os << "<parm name=\"" << parmEntry->second.name << "\" default=\"";
      if ( valDefault == 0. ) os << fixed << setprecision(1);
      else if ( abs(valDefault) < 0.001 ) os << scientific << setprecision(4);
      else if ( abs(valDefault) < 0.1 ) os << fixed << setprecision(7);
      else if ( abs(valDefault) < 1000. ) os << fixed << setprecision(5);
      else if ( abs(valDefault) < 1000000. ) os << fixed << setprecision(3);
      else os << scientific << setprecision(4);
      os << valDefault << "\"";
      if (parmEntry->second.hasMin) {
        os << " min=\"";
        valDefault = parmEntry->second.valMin;
        if ( valDefault == 0. ) os << fixed << setprecision(1);
        else if ( abs(valDefault) < 0.001) os << scientific << setprecision(4);
        else if ( abs(valDefault) < 0.1 ) os << fixed << setprecision(7);
        else if ( abs(valDefault) < 1000. ) os << fixed << setprecision(5);
        else if ( abs(valDefault) < 1000000. ) os << fixed << setprecision(3);
        else os << scientific << setprecision(4);
        os << valDefault << "\"";
      }
      if (parmEntry->second.hasMax) {
        os << " max=\"";
        valDefault = parmEntry->second.valMax;
        if ( valDefault == 0. ) os << fixed << setprecision(1);
        else if ( abs(valDefault) < 0.001) os << scientific << setprecision(4);
        else if ( abs(valDefault) < 0.1 ) os << fixed << setprecision(7);
        else if ( abs(valDefault) < 1000. ) os << fixed << setprecision(5);
        else if ( abs(valDefault) < 1000000. ) os << fixed << setprecision(3);
        else os << scientific << setprecision(4);
        os << valDefault << "\"";
      }
      os << "></parm>" << endl;
      ++parmEntry;

    // Else check if word is next, and if so print it.
    } else  if ( wordEntry != words.end()
      && ( fvecEntry == fvecs.end() || wordEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || wordEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || wordEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || wordEntry->first < wvecEntry->first )
      ) {
      string valDefault = wordEntry->second.valDefault;
      os << "<word name=\"" << wordEntry->second.name << "\" default=\""
         << valDefault << "\"></word>" << endl;
      ++wordEntry;

    // Else check if fvec is next, and if so print it.
    } else if ( fvecEntry != fvecs.end()
      && ( mvecEntry == mvecs.end() || fvecEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || fvecEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || fvecEntry->first < wvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      vector<bool> valDefault = fvecEntry->second.valDefault;
      os << "<fvec name=\"" << fvecEntry->second.name << "\" default=\"{";
      if (valDefault.size() > 0) {
        for (vector<bool>::iterator val = valDefault.begin();
             val != --valDefault.end(); ++val) os << state[*val] << ",";
        os << state[*(--valDefault.end())] << "}\"></fvec>" << endl;
      } else os << "}\"></fvec>" << endl;
      ++fvecEntry;

    // Else check if mvec is next, and if so print it.
    } else if ( mvecEntry != mvecs.end()
      && ( pvecEntry == pvecs.end() || mvecEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || mvecEntry->first < wvecEntry->first )
      ) {
      vector<int> valDefault = mvecEntry->second.valDefault;
      os << "<mvec name=\"" << mvecEntry->second.name << "\" default=\"{";
      if (valDefault.size() > 0) {
        for (vector<int>::iterator val = valDefault.begin();
             val != --valDefault.end(); ++val) os << *val << ",";
        os << *(--valDefault.end()) << "}\"";
      } else os << "}\">";
      if (mvecEntry->second.hasMin ) os << " min=\""
        << mvecEntry->second.valMin << "\"";
      if (mvecEntry->second.hasMax ) os << " max=\""
        << mvecEntry->second.valMax << "\"";
      os << "></mvec>" << endl;
      ++mvecEntry;

    // Else check if pvec is next; print fixed or scientific.
    } else if ( pvecEntry != pvecs.end()
      && ( wvecEntry == wvecs.end() || pvecEntry->first < wvecEntry->first )
      ) {
      vector<double> valDefault = pvecEntry->second.valDefault;
      os << "<pvec name=\"" << pvecEntry->second.name << "\" default=\"{";
      if (valDefault.size() > 0) {
        for (vector<double>::iterator val = valDefault.begin();
             val != --valDefault.end(); ++val) {
          if ( *val == 0. ) os << fixed << setprecision(1);
          else if ( abs(*val) < 0.001 ) os << scientific << setprecision(4);
          else if ( abs(*val) < 0.1 ) os << fixed << setprecision(7);
          else if ( abs(*val) < 1000. ) os << fixed << setprecision(5);
          else if ( abs(*val) < 1000000. ) os << fixed << setprecision(3);
          else os << scientific << setprecision(4);
          os << *val << ",";
        }
        os << *(--valDefault.end()) << "}\"";
      } else os << "}\">";
      if (pvecEntry->second.hasMin ) {
        double valLocal = pvecEntry->second.valMin;
        os << " min=\"";
        if ( valLocal == 0. ) os << fixed << setprecision(1);
        else if ( abs(valLocal) < 0.001 ) os << scientific << setprecision(4);
        else if ( abs(valLocal) < 0.1 ) os << fixed << setprecision(7);
        else if ( abs(valLocal) < 1000. ) os << fixed << setprecision(5);
        else if ( abs(valLocal) < 1000000. ) os << fixed << setprecision(3);
        else os << scientific << setprecision(4);
        os << valLocal << "\"";
      }
      if (pvecEntry->second.hasMax ) {
        double valLocal = pvecEntry->second.valMax;
        os << " max=\"";
        if ( valLocal == 0. ) os << fixed << setprecision(1);
        else if ( abs(valLocal) < 0.001 ) os << scientific << setprecision(4);
        else if ( abs(valLocal) < 0.1 ) os << fixed << setprecision(7);
        else if ( abs(valLocal) < 1000. ) os << fixed << setprecision(5);
        else if ( abs(valLocal) < 1000000. ) os << fixed << setprecision(3);
        else os << scientific << setprecision(4);
        os << valLocal << "\"";
      }
      os << "></pvec>" <<       endl;
      ++pvecEntry;

    // Else print wvec.
    } else {
      vector<string> valDefault = wvecEntry->second.valDefault;
      os << "<wvec name=\"" << wvecEntry->second.name << "\" default=\"{";
      if (valDefault.size() > 0) {
        for (vector<string>::iterator val = valDefault.begin();
             val != --valDefault.end(); ++val) os << *val << ",";
        os << *(--valDefault.end()) << "}\">";
      } else os << "}\">";
      os << "</wvec>" << endl;
      ++wvecEntry;
   }
  } ;

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Print out table of database in lexigraphical order.

void Settings::list(bool doListAll,  bool doListString, string match) {

  // Table header; output for bool as off/on.
  if (doListAll)
    cout << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec "
       << "+ PVec + WVec Settings (all)  ---------------------------* \n";
  else if (!doListString)
    cout << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec "
       << "+ PVec + WVec Settings (changes only)  ------------------* \n" ;
  else
    cout << "\n *-------  PYTHIA Flag + Mode + Parm + Word + FVec + MVec "
       << "+ PVec + WVec Settings (with requested string) ----------* \n" ;
  cout << " |                                                           "
       << "                                                      | \n"
       << " | Name                                          |           "
       << "           Now |      Default         Min         Max | \n"
       << " |                                               |           "
       << "               |                                      | \n";

  // Convert input string to lowercase for match.
  toLowerRep(match);
  if (match == "") match = "             ";

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();
  map<string, FVec>::iterator fvecEntry = fvecs.begin();
  map<string, MVec>::iterator mvecEntry = mvecs.begin();
  map<string, PVec>::iterator pvecEntry = pvecs.begin();
  map<string, WVec>::iterator wvecEntry = wvecs.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end()
      || parmEntry != parms.end() || wordEntry != words.end()
      || fvecEntry != fvecs.end() || mvecEntry != mvecs.end()
      || pvecEntry != pvecs.end() || wvecEntry != wvecs.end() ) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end()
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first )
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || flagEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || flagEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || flagEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || flagEntry->first < wvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( doListAll || (!doListString && valNow != valDefault)
        || (doListString && flagEntry->first.find(match) != string::npos) )
        cout << " | " << setw(45) << left
             << flagEntry->second.name << " | " << setw(24) << right
             << state[valNow] << " | " << setw(12) << state[valDefault]
             << "                         | \n";
      ++flagEntry;

    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end()
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || modeEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || modeEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || modeEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || modeEntry->first < wvecEntry->first )
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( doListAll || (!doListString && valNow != valDefault)
        || (doListString && modeEntry->first.find(match) != string::npos) ) {
        cout << " | " << setw(45) << left
             << modeEntry->second.name << " | " << setw(24) << right
             << valNow << " | " << setw(12) << valDefault;
        if (modeEntry->second.hasMin)
          cout << setw(12) << modeEntry->second.valMin;
        else cout << "            ";
        if (modeEntry->second.hasMax)
          cout << setw(12) << modeEntry->second.valMax;
        else cout << "            ";
        cout << " | \n";
      }
      ++modeEntry;

    // Else check if parm is next, and if so print it; fixed or scientific.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first )
      && ( fvecEntry == fvecs.end() || parmEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || parmEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || parmEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || parmEntry->first < wvecEntry->first )
      ) {
      double valNow = parmEntry->second.valNow;
      double valDefault = parmEntry->second.valDefault;
      if ( doListAll || (!doListString && valNow != valDefault )
        || (doListString && parmEntry->first.find(match) != string::npos) ) {
        cout << " | " << setw(45) << left
             << parmEntry->second.name << right << " |             ";
        for (int i = 0; i < 4; ++i) {
          if (i == 1) valNow = valDefault;
          if (i == 2) valNow = parmEntry->second.valMin;
          if (i == 3) valNow = parmEntry->second.valMax;
          if ( (i == 2 && !parmEntry->second.hasMin)
            || (i == 3 && !parmEntry->second.hasMax) )
            cout << "            ";
          else if ( valNow == 0. )
            cout << fixed << setprecision(1) << setw(12) << valNow;
          else if ( abs(valNow) < 0.001 )
            cout << scientific << setprecision(4) << setw(12) << valNow;
          else if ( abs(valNow) < 0.1 )
            cout << fixed << setprecision(7) << setw(12) << valNow;
          else if ( abs(valNow) < 1000. )
            cout << fixed << setprecision(5) << setw(12) << valNow;
          else if ( abs(valNow) < 1000000. )
            cout << fixed << setprecision(3) << setw(12) << valNow;
          else
            cout << scientific << setprecision(4) << setw(12) << valNow;
          if (i == 0) cout << " | ";
        }
        cout << " | \n";
      }
      ++parmEntry;

    // Else check if word is next, and if so print it.
    } else  if ( wordEntry != words.end()
      && ( fvecEntry == fvecs.end() || wordEntry->first < fvecEntry->first )
      && ( mvecEntry == mvecs.end() || wordEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || wordEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || wordEntry->first < wvecEntry->first )
      ) {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault;
      int blankLeft = max(0, 60 - max(24, int(valNow.length()) )
        - max(12, int(valDefault.length()) ) );
      string blankPad( blankLeft, ' ');
      if ( doListAll || (!doListString && valNow != valDefault)
        || (doListString && wordEntry->first.find(match) != string::npos) )
        cout << " | " << setw(45) << left
             << wordEntry->second.name << " | " << setw(24) << right
             << valNow << " | " << setw(12) << valDefault << blankPad
             << " | \n";
      ++wordEntry;

    // Else check if fvec is next, and if so print it.
    } else if ( fvecEntry != fvecs.end()
      && ( mvecEntry == mvecs.end() || fvecEntry->first < mvecEntry->first )
      && ( pvecEntry == pvecs.end() || fvecEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || fvecEntry->first < wvecEntry->first )
      ) {
      string state[2] = {"off", "on"};
      vector<bool> valsNow = fvecEntry->second.valNow;
      vector<bool> valsDefault = fvecEntry->second.valDefault;
      bool valNow(false), valDefault(false);
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && fvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            cout << " | " << setw(45) << left
                 << fvecEntry->second.name << right << " |             ";
          else
            cout << " | " << setw(45) << " " << right << " |             ";
          for (int j = 0; j < 4; ++j) {
            if (i < valsNow.size()) valNow = valsNow[i];
            if (i < valsDefault.size()) valDefault = valsDefault[i];
            if (j == 1) valNow = valDefault;
            if ( (j == 0 && i >= valsNow.size())
                 || (j == 1 && i >= valsDefault.size()) || (j > 1) )
              cout << "            ";
            else cout << setw(12) << state[valNow];
            if (j == 0) cout << " | ";
          }
          cout << " | \n";
        }
      }
      ++fvecEntry;

    // Else check if mvec is next, and if so print it.
    } else if ( mvecEntry != mvecs.end()
      && ( pvecEntry == pvecs.end() || mvecEntry->first < pvecEntry->first )
      && ( wvecEntry == wvecs.end() || mvecEntry->first < wvecEntry->first )
      ) {
      vector<int> valsNow = mvecEntry->second.valNow;
      vector<int> valsDefault = mvecEntry->second.valDefault;
      int valNow(0), valDefault(0);
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && mvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            cout << " | " << setw(45) << left
                 << mvecEntry->second.name << right << " |             ";
          else
            cout << " | " << setw(45) << " " << right << " |             ";
          for (int j = 0; j < 4; ++j) {
            if (i < valsNow.size()) valNow = valsNow[i];
            if (i < valsDefault.size()) valDefault = valsDefault[i];
            if (j == 1) valNow = valDefault;
            if (j == 2) valNow = mvecEntry->second.valMin;
            if (j == 3) valNow = mvecEntry->second.valMax;
            if ( (j == 0 && i >= valsNow.size())
                 || (j == 1 && i >= valsDefault.size())
                 || (j == 2 && !mvecEntry->second.hasMin)
                 || (j == 3 && !mvecEntry->second.hasMax) )
              cout << "            ";
            else cout << setw(12) << valNow;
            if (j == 0) cout << " | ";
          }
          cout << " | \n";
        }
      }
      ++mvecEntry;

    // Else check if pvec is next; print fixed or scientific.
    } else if ( pvecEntry != pvecs.end()
      && ( wvecEntry == wvecs.end() || pvecEntry->first < wvecEntry->first )
      ) {
      vector<double> valsNow = pvecEntry->second.valNow;
      vector<double> valsDefault = pvecEntry->second.valDefault;
      double valNow(0), valDefault(0);
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && pvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            cout << " | " << setw(45) << left
                 << pvecEntry->second.name << right << " |             ";
          else
            cout << " | " << setw(45) << " " << right << " |             ";
          for (int j = 0; j < 4; ++j) {
            if (i < valsNow.size()) valNow = valsNow[i];
            if (i < valsDefault.size()) valDefault = valsDefault[i];
            if (j == 1) valNow = valDefault;
            if (j == 2) valNow = pvecEntry->second.valMin;
            if (j == 3) valNow = pvecEntry->second.valMax;
            if ( (j == 0 && i >= valsNow.size())
                 || (j == 1 && i >= valsDefault.size())
                 || (j == 2 && !pvecEntry->second.hasMin)
                 || (j == 3 && !pvecEntry->second.hasMax) )
              cout << "            ";
            else if ( valNow == 0. )
              cout << fixed << setprecision(1) << setw(12) << valNow;
            else if ( abs(valNow) < 0.001 )
              cout << scientific << setprecision(4) << setw(12) << valNow;
            else if ( abs(valNow) < 0.1 )
              cout << fixed << setprecision(7) << setw(12) << valNow;
            else if ( abs(valNow) < 1000. )
              cout << fixed << setprecision(5) << setw(12) << valNow;
            else if ( abs(valNow) < 1000000. )
              cout << fixed << setprecision(3) << setw(12) << valNow;
            else
              cout << scientific << setprecision(4) << setw(12) << valNow;
            if (j == 0) cout << " | ";
          }
          cout << " | \n";
        }
      }
      ++pvecEntry;

    // Else print wvec.
    } else {
      vector<string> valsNow = wvecEntry->second.valNow;
      vector<string> valsDefault = wvecEntry->second.valDefault;
      if ( wvecEntry->first.length() > 10 &&
           wvecEntry->first.substr(0, 10) == "init:reuse" )
        valsDefault = valsNow;
      if ( doListAll || (!doListString && valsNow != valsDefault )
        || (doListString && wvecEntry->first.find(match) != string::npos) ) {
        for (unsigned int i = 0; i < valsNow.size() || i < valsDefault.size();
             ++i) {
          if ( i == 0 )
            cout << " | " << setw(45) << left
                 << wvecEntry->second.name << right << " | ";
          else
            cout << " | " << setw(45) << " " << right << " | ";
          string valNow =  (i < valsNow.size()) ? valsNow[i] : " ";
          string valDefault = (i < valsDefault.size()) ? valsDefault[i] : " ";
          int blankLeft = max(0, 60 - max(24, int(valNow.length()) )
            - max(12, int(valDefault.length()) ) );
          string blankPad( blankLeft, ' ');
          cout << setw(24) << right << valNow << " | " << setw(12)
               << valDefault << blankPad << " | \n";
        }
      }
      ++wvecEntry;
    }
  } ;

  // End of loop over database contents.
  cout << " |                                                           "
       << "                                                      | \n"
       << " *-------  End PYTHIA Flag + Mode + Parm + Word + FVec + MVec "
       << "+ PVec + WVec Settings  -----------------------------* " << endl;

}

//--------------------------------------------------------------------------

// Give back current value(s) as a string, whatever the type.

string Settings::output(string keyIn, bool fullLine) {

  // Default string echoes input key =.
  string outVal = (fullLine) ? " " + keyIn + " = " : "";

  // Identify flag, mode, parm or word, and convert to string.
  if (isFlag(keyIn)) {
    outVal += (flag(keyIn)) ? "true" : "false";
  } else if (isMode(keyIn)) {
    ostringstream ostr;
    ostr << mode(keyIn);
    outVal += ostr.str();
  } else if (isParm(keyIn)) {
    ostringstream ostr;
    ostr << scientific << setprecision(5) << parm(keyIn);
    outVal += ostr.str();
  } else if (isWord(keyIn)) {
    outVal += word(keyIn);

  // Identify fvec, mvec, pvec or wvec, and convert to string.
  } else if (isFVec(keyIn)) {
    vector<bool> outVec = fvec(keyIn);
    for (int i = 0; i < int(outVec.size()); ++i) {
      outVal += (outVec[i]) ? "true" : "false";
      if (i != int(outVec.size()) - 1) outVal += "  ";
    }
  } else if (isMVec(keyIn)) {
    vector<int> outVec = mvec(keyIn);
    for (int i = 0; i < int(outVec.size()); ++i) {
      ostringstream ostr;
      ostr << outVec[i];
      outVal +=  ostr.str();
      if (i != int(outVec.size()) - 1) outVal += "  ";
    }
  } else if (isPVec(keyIn)) {
    vector<double> outVec = pvec(keyIn);
    for (int i = 0; i < int(outVec.size()); ++i) {
      ostringstream ostr;
      ostr << scientific << setprecision(5) << outVec[i];
      outVal +=  ostr.str();
      if (i != int(outVec.size()) - 1) outVal += "  ";
    }
  } else if (isWVec(keyIn)) {
    vector<string> outVec = wvec(keyIn);
    for (int i = 0; i < int(outVec.size()); ++i) {
      outVal +=  outVec[i];
      if (i != int(outVec.size()) - 1) outVal += "  ";
    }

  // Default value, possible endline and done.
  } else outVal += "unknown";
  if (fullLine) outVal += "\n";
  return outVal;

}

//--------------------------------------------------------------------------

// Reset all values to their defaults.

void Settings::resetAll() {

  // Loop through the flags table, resetting all entries.
  for (map<string, Flag>::iterator flagEntry = flags.begin();
    flagEntry != flags.end(); ++flagEntry) {
    string name = flagEntry->first;
    resetFlag(name);
  }

  // Loop through the modes table, resetting all entries.
  for (map<string, Mode>::iterator modeEntry = modes.begin();
    modeEntry != modes.end(); ++modeEntry) {
    string name = modeEntry->first;
    resetMode(name);
  }

  // Loop through the parms table, resetting all entries.
  for (map<string, Parm>::iterator parmEntry = parms.begin();
    parmEntry != parms.end(); ++parmEntry) {
    string name = parmEntry->first;
    resetParm(name);
  }

  // Loop through the words table, resetting all entries.
  for (map<string, Word>::iterator wordEntry = words.begin();
    wordEntry != words.end(); ++wordEntry) {
    string name = wordEntry->first;
    resetWord(name);
  }

  // Loop through the fvecs table, resetting all entries.
  for (map<string, FVec>::iterator fvecEntry = fvecs.begin();
    fvecEntry != fvecs.end(); ++fvecEntry) {
    string name = fvecEntry->first;
    resetFVec(name);
  }

  // Loop through the mvecs table, resetting all entries.
  for (map<string, MVec>::iterator mvecEntry = mvecs.begin();
    mvecEntry != mvecs.end(); ++mvecEntry) {
    string name = mvecEntry->first;
    resetMVec(name);
  }

  // Loop through the pvecs table, resetting all entries.
  for (map<string, PVec>::iterator pvecEntry = pvecs.begin();
    pvecEntry != pvecs.end(); ++pvecEntry) {
    string name = pvecEntry->first;
    resetPVec(name);
  }

  // Loop through the wvecs table, resetting all entries.
  for (map<string, WVec>::iterator wvecEntry = wvecs.begin();
    wvecEntry != wvecs.end(); ++wvecEntry) {
    string name = wvecEntry->first;
    resetWVec(name);
  }

}

//--------------------------------------------------------------------------

// Give back current value, with check that key exists.

bool Settings::flag(string keyIn) {
  if (isFlag(keyIn)) return flags[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return false;
}

int Settings::mode(string keyIn) {
  if (isMode(keyIn)) return modes[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return 0;
}

double Settings::parm(string keyIn) {
  if (isParm(keyIn)) return parms[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return 0.;
}

string Settings::word(string keyIn) {
  if (isWord(keyIn)) return words[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return " ";
}

vector<bool> Settings::fvec(string keyIn) {
  if (isFVec(keyIn)) return fvecs[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<bool>(1, false);
}

vector<int> Settings::mvec(string keyIn) {
  if (isMVec(keyIn)) return mvecs[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<int>(1, 0);
}

vector<double> Settings::pvec(string keyIn) {
  if (isPVec(keyIn)) return pvecs[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<double>(1, 0.);
}

vector<string> Settings::wvec(string keyIn) {
  if (isWVec(keyIn)) return wvecs[toLower(keyIn)].valNow;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<string>(1, " ");
}

//--------------------------------------------------------------------------

// Give back default value, with check that key exists.

bool Settings::flagDefault(string keyIn) {
  if (isFlag(keyIn)) return flags[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return false;
}

int Settings::modeDefault(string keyIn) {
  if (isMode(keyIn)) return modes[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return 0;
}

double Settings::parmDefault(string keyIn) {
  if (isParm(keyIn)) return parms[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return 0.;
}

string Settings::wordDefault(string keyIn) {
  if (isWord(keyIn)) return words[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return " ";
}

vector<bool> Settings::fvecDefault(string keyIn) {
  if (isFVec(keyIn)) return fvecs[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<bool>(1, false);
}

vector<int> Settings::mvecDefault(string keyIn) {
  if (isMVec(keyIn)) return mvecs[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<int>(1, 0);
}

vector<double> Settings::pvecDefault(string keyIn) {
  if (isPVec(keyIn)) return pvecs[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<double>(1, 0.);
}

vector<string> Settings::wvecDefault(string keyIn) {
  if (isWVec(keyIn)) return wvecs[toLower(keyIn)].valDefault;
  loggerPtr->ERROR_MSG("unknown key", keyIn);
  return vector<string>(1, " ");
}

//--------------------------------------------------------------------------

// Get a map of entries whose names contain the string "match".

map<string, Flag> Settings::getFlagMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, Flag> flagMap;
  // Loop over the flag map (using iterator).
  for (map<string,Flag>::iterator flagEntry = flags.begin();
       flagEntry != flags.end(); ++flagEntry)
    if (flagEntry->first.find(match) != string::npos)
      flagMap[flagEntry->first] = flagEntry->second;
  return flagMap;
}

map<string, Mode> Settings::getModeMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, Mode> modeMap;
  // Loop over the mode map (using iterator).
  for (map<string,Mode>::iterator modeEntry = modes.begin();
       modeEntry != modes.end(); ++modeEntry)
    if (modeEntry->first.find(match) != string::npos)
      modeMap[modeEntry->first] = modeEntry->second;
  return modeMap;
}

map<string, Parm> Settings::getParmMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, Parm> parmMap;
  // Loop over the parm map (using iterator).
  for (map<string,Parm>::iterator parmEntry = parms.begin();
       parmEntry != parms.end(); ++parmEntry)
    if (parmEntry->first.find(match) != string::npos)
      parmMap[parmEntry->first] = parmEntry->second;
  return parmMap;
}

map<string, Word> Settings::getWordMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, Word> wordMap;
  // Loop over the word map (using iterator).
  for (map<string,Word>::iterator wordEntry = words.begin();
       wordEntry != words.end(); ++wordEntry)
    if (wordEntry->first.find(match) != string::npos)
      wordMap[wordEntry->first] = wordEntry->second;
  return wordMap;
}

map<string, FVec> Settings::getFVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, FVec> fvecMap;
  // Loop over the fvec map (using iterator).
  for (map<string,FVec>::iterator fvecEntry = fvecs.begin();
       fvecEntry != fvecs.end(); ++fvecEntry)
    if (fvecEntry->first.find(match) != string::npos)
      fvecMap[fvecEntry->first] = fvecEntry->second;
  return fvecMap;
}

map<string, MVec> Settings::getMVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, MVec> mvecMap;
  // Loop over the mvec map (using iterator).
  for (map<string,MVec>::iterator mvecEntry = mvecs.begin();
       mvecEntry != mvecs.end(); ++mvecEntry)
    if (mvecEntry->first.find(match) != string::npos)
      mvecMap[mvecEntry->first] = mvecEntry->second;
  return mvecMap;
}

map<string, PVec> Settings::getPVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, PVec> pvecMap;
  // Loop over the pvec map (using iterator).
  for (map<string,PVec>::iterator pvecEntry = pvecs.begin();
       pvecEntry != pvecs.end(); ++pvecEntry)
    if (pvecEntry->first.find(match) != string::npos)
      pvecMap[pvecEntry->first] = pvecEntry->second;
  return pvecMap;
}

map<string, WVec> Settings::getWVecMap(string match) {
  // Make the match string lower case. Start with an empty map.
  toLowerRep(match);
  map<string, WVec> wvecMap;
  // Loop over the wvec map (using iterator).
  for (map<string,WVec>::iterator wvecEntry = wvecs.begin();
       wvecEntry != wvecs.end(); ++wvecEntry)
    if (wvecEntry->first.find(match) != string::npos)
      wvecMap[wvecEntry->first] = wvecEntry->second;
  return wvecMap;
}

//--------------------------------------------------------------------------

// Change current value. Respect limits unless force==true.
// If key not recognised, add new key if force==true, otherwise ignore.

void Settings::flag(string keyIn, bool nowIn, bool force) {
  string keyLower = toLower(keyIn);
  if (isFlag(keyLower)) flags[keyLower].valNow = nowIn;
  else if (force) addFlag( keyIn, nowIn);
  // Print:quiet  triggers a whole set of changes.
  if (keyLower == "print:quiet") printQuiet( nowIn);
}

bool Settings::mode(string keyIn, int nowIn, bool force) {
  if (isMode(keyIn)) {
    string keyLower = toLower(keyIn);
    Mode& modeNow = modes[keyLower];
    // Fail if values are outside range.
    if (!force && ((modeNow.hasMin && nowIn < modeNow.valMin)
                || (modeNow.hasMax && nowIn > modeNow.valMax)) ) {
      loggerPtr->ERROR_MSG("value is out of range", keyIn, true);
      return false;
    }
    else modeNow.valNow = nowIn;
    // Tunes each trigger a whole set of changes.
    if (keyLower == "tune:ee") initTuneEE(modeNow.valNow);
    if (keyLower == "tune:pp") initTunePP(modeNow.valNow);
    if (keyLower == "vincia:tune") initTuneVincia(modeNow.valNow);
  }
  else if (force)
    addMode(keyIn, nowIn, false, false, 0, 0);

  return true;
}

bool Settings::parm(string keyIn, double nowIn, bool force) {
  if (isParm(keyIn)) {
    Parm& parmNow = parms[toLower(keyIn)];
    if (!force && ((parmNow.hasMin && nowIn < parmNow.valMin)
                || (parmNow.hasMax && nowIn > parmNow.valMax)) ) {
      loggerPtr->ERROR_MSG("value is out of range", keyIn, true);
      return false;
    }
    else parmNow.valNow = nowIn;
  }
  else if (force)
    addParm(keyIn, nowIn, false, false, 0., 0.);

  return true;
}

void Settings::word(string keyIn, string nowIn, bool force) {
  if (isWord(keyIn)) words[toLower(keyIn)].valNow = nowIn;
  else if (force) addWord(keyIn, nowIn);
}

void Settings::fvec(string keyIn, vector<bool> nowIn, bool force) {
  if (isFVec(keyIn)) {
    FVec& fvecNow = fvecs[toLower(keyIn)];
    fvecNow.valNow.clear();
    for (vector<bool>::iterator now = nowIn.begin();
        now != nowIn.end(); now++)
      fvecNow.valNow.push_back(*now);
  }
  else if (force)
    addFVec(keyIn, nowIn);
}

bool Settings::mvec(string keyIn, vector<int> nowIn, bool force) {
  if (isMVec(keyIn)) {
    MVec& mvecNow = mvecs[toLower(keyIn)];
    mvecNow.valNow.clear();
    for (vector<int>::iterator now = nowIn.begin();
        now != nowIn.end(); now++) {
      if (!force && ((mvecNow.hasMin && *now < mvecNow.valMin)
                  || (mvecNow.hasMax && *now > mvecNow.valMax)) ) {
        loggerPtr->ERROR_MSG("value is out of range", keyIn, true);
        return false;
      }
      else mvecNow.valNow.push_back(*now);
    }
  }
  else if (force)
    addMVec(keyIn, nowIn, false, false, 0, 0);

  return true;
}

bool Settings::pvec(string keyIn, vector<double> nowIn, bool force) {
  if (isPVec(keyIn)) {
    PVec& pvecNow = pvecs[toLower(keyIn)];
    pvecNow.valNow.clear();
    for (vector<double>::iterator now = nowIn.begin();
        now != nowIn.end(); now++) {
      if (!force && ((pvecNow.hasMin && *now < pvecNow.valMin)
                  || (pvecNow.hasMax && *now > pvecNow.valMax)) ) {
        loggerPtr->ERROR_MSG("value is out of range", keyIn, true);
        return false;
      }
      else pvecNow.valNow.push_back(*now);
    }
  }
  else if (force)
    addPVec(keyIn, nowIn, false, false, 0., 0.);

  return true;
}

void Settings::wvec(string keyIn, vector<string> nowIn, bool force) {
  if (isWVec(keyIn)) {
    WVec& wvecNow = wvecs[toLower(keyIn)];
    wvecNow.valNow.clear();
    for (vector<string>::iterator now = nowIn.begin();
        now != nowIn.end(); now++)
      wvecNow.valNow.push_back(*now);
  }
  else if (force) addWVec(keyIn, nowIn);
  // Register settings from plugin libraries immediately.
  if (toLower(keyIn) == "init:plugins")
    for (string lib : nowIn)
      registerPluginLibrary(lib.substr(0, lib.find("::")));
}

//--------------------------------------------------------------------------

// Restore current value to default.

void Settings::resetFlag(string keyIn) {
  if (isFlag(keyIn)) flags[toLower(keyIn)].valNow
    = flags[toLower(keyIn)].valDefault ;
}

void Settings::resetMode(string keyIn) {
  string keyLower = toLower(keyIn);
  if (isMode(keyIn)) modes[keyLower].valNow
    = modes[toLower(keyIn)].valDefault ;
  // For Tune:ee and Tune:pp must also restore variables involved in tunes.
  if (keyLower == "tune:ee") readString("include = tunes/Reset-ee.cmnd", true);
  if (keyLower == "tune:pp") readString("include = tunes/Reset-pp.cmnd", true);
}

void Settings::resetParm(string keyIn) {
  if (isParm(keyIn)) parms[toLower(keyIn)].valNow
    = parms[toLower(keyIn)].valDefault ;
}

void Settings::resetWord(string keyIn) {
  if (isWord(keyIn)) words[toLower(keyIn)].valNow
    = words[toLower(keyIn)].valDefault ;
}

void Settings::resetFVec(string keyIn) {
  if (isFVec(keyIn)) fvecs[toLower(keyIn)].valNow
    = fvecs[toLower(keyIn)].valDefault ;
}

void Settings::resetMVec(string keyIn) {
  if (isMVec(keyIn)) mvecs[toLower(keyIn)].valNow
    = mvecs[toLower(keyIn)].valDefault ;
}

void Settings::resetPVec(string keyIn) {
  if (isPVec(keyIn)) pvecs[toLower(keyIn)].valNow
    = pvecs[toLower(keyIn)].valDefault ;
}

void Settings::resetWVec(string keyIn) {
  if (isWVec(keyIn)) wvecs[toLower(keyIn)].valNow
    = wvecs[toLower(keyIn)].valDefault ;
}

//--------------------------------------------------------------------------

// Check whether any other processes than SoftQCD and LowEnergyQCD are
// switched on. Note that Les Houches input has to be checked separately.

bool Settings::hasHardProc() {

  // List of (most?) process name groups, in lowercase. Special cases.
  string flagList[26] = { "hardqcd", "promptphoton", "weakbosonexchange",
    "weaksingleboson", "weakdoubleboson", "weakbosonandparton",
    "photoncollision", "photonparton", "onia:all", "charmonium:all",
    "bottomonium:all", "top", "fourthbottom", "fourthtop", "fourthpair",
    "higgssm", "higgsbsm", "susy", "newgaugeboson", "leftrightsymmetry",
    "leptoquark", "excitedfermion", "contactinteractions", "hiddenvalley",
    "extradimensions", "dm:" };
  int sizeList = 26;
  string flagExclude[2] = { "extradimensionsg*:vlvl", "higgssm:nlowidths"};
  int sizeExclude = 2;

  // Loop over the flag map (using iterator), and process names.
  for (map<string,Flag>::iterator flagEntry = flags.begin();
    flagEntry != flags.end(); ++flagEntry) {
    string flagName = flagEntry->first;
    bool doExclude = false;
    for (int i = 0; i < sizeExclude; ++i)
      if (flagName.find( flagExclude[i]) != string::npos) doExclude = true;
    if (doExclude) continue;
    for (int i = 0; i < sizeList; ++i)
      if (flagName.find( flagList[i]) != string::npos
      && flagEntry->second.valNow == true) return true;
  }

  // Done without having found a non-SoftQCD/LowEnergyQCD process on.
  return false;

}

//--------------------------------------------------------------------------

// Regulate level of printout by overall change of settings.

void Settings::printQuiet(bool quiet) {

  // Switch off as much output as possible.
  if (quiet) {
    flag("Init:showProcesses",               false );
    flag("Init:showMultipartonInteractions", false );
    flag("Init:showChangedSettings",         false );
    flag("Init:showAllSettings",             false );
    flag("Init:showChangedParticleData",     false );
    flag("Init:showChangedResonanceData",    false );
    flag("Init:showAllParticleData",         false );
    mode("Init:showOneParticleData",             0 );
    mode("Next:numberCount",                     0 );
    mode("Next:numberShowLHA",                   0 );
    mode("Next:numberShowInfo",                  0 );
    mode("Next:numberShowProcess",               0 );
    mode("Next:numberShowEvent",                 0 );
    flag("Print:errors",                     false );

  // Restore ouput settings to default.
  } else {
    resetFlag("Init:showProcesses");
    resetFlag("Init:showMultipartonInteractions");
    resetFlag("Init:showChangedSettings");
    resetFlag("Init:showAllSettings");
    resetFlag("Init:showChangedParticleData");
    resetFlag("Init:showChangedResonanceData");
    resetFlag("Init:showAllParticleData");
    resetMode("Init:showOneParticleData");
    resetMode("Next:numberCount");
    resetMode("Next:numberShowLHA");
    resetMode("Next:numberShowInfo");
    resetMode("Next:numberShowProcess");
    resetMode("Next:numberShowEvent");
  }

}

//--------------------------------------------------------------------------

// Set the values related to a tune of e+e- data,
// i.e. mainly for final-state radiation and hadronization.

void Settings::initTuneEE(int eeTune) {

  // Map the tune files to integer values.
  vector<string> tunes = {
    "Reset-ee", "", "OldJETSET", "Montull2007", "Hoeth2009", "Skands2013",
    "Fischer2013-1", "Fischer2013-2", "Monash2013-ee"};
  if (eeTune + 1 < (int)tunes.size() && tunes[eeTune + 1] != "")
    readString("include = tunes/" + tunes[eeTune + 1] + ".cmnd", true);

}

//--------------------------------------------------------------------------

// Set the values related to a tune of pp/ppbar data,
// i.e. mainly for initial-state radiation and multiparton interactions.

void Settings::initTunePP(int ppTune) {

  // Map the tune files to integer values.
  vector<string> tunes = {
    "Rest-pp", "", "OldIsrMpi", "Skands2009", "Tune2C", "Tune2M", "Tune4C",
    "Tune4Cx", "ATLAS-MB-A2-CTEQ6L1", "ATLAS-MB-A2-MSTW2008LO",
    "ATLAS-UE-AU2-CTEQ6L1", "ATLAS-UE-AU2-MSTW2008LO", "ATLAS-UE-AU2-CT10",
    "ATLAS-UE-AU2-MRST2007LOx", "ATLAS-UE-AU2-MRST2007LOxx", "Monash2013",
    "CMS-CUETP8S1-CTEQ6L1", "CMS-CUETP8S1-HERAPDF1", "ATLAS-AZ",
    "CMS-CUETP8M1-NNPDF23LO", "ATLAS-A14-CTEQL1", "ATLAS-A14-MSTW2008LO",
    "ATLAS-A14-NNPDF23LO", "ATLAS-A14-HERAPDF15LO", "ATLAS-A14-v+1",
    "ATLAS-A14-v-1", "ATLAS-A14-v+2", "ATLAS-A14-v-2", "ATLAS-A14-v+3a",
    "ATLAS-A14-v-3a", "ATLAS-A14-v+3b", "ATLAS-A14-v-3b", "ATLAS-A14-v+3c",
    "ATLAS-A14-v-3c"};
  if (ppTune + 1 < (int)tunes.size() && tunes[ppTune + 1] != "")
    readString("include = tunes/" + tunes[ppTune + 1] + ".cmnd", true);

}

//--------------------------------------------------------------------------

// Set the values related to a tune of Vincia.

void Settings::initTuneVincia(int vinciaTune) {

  // Currently only a single tune.
  if (vinciaTune == 0)
    readString("include = tunes/VinciaDefault.cmnd", true);

}

//--------------------------------------------------------------------------

// Allow several alternative inputs for true/false.

bool Settings::boolString(string tag) {

  string tagLow = toLower(tag);
  return ( tagLow == "true" || tagLow == "1" || tagLow == "on"
  || tagLow == "yes" || tagLow == "ok" );

}

//--------------------------------------------------------------------------

// Extract XML value string following XML attribute.

string Settings::attributeValue(string line, string attribute) {

  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

//--------------------------------------------------------------------------

// Extract XML bool value following XML attribute.

bool Settings::boolAttributeValue(string line, string attribute) {

  string valString = attributeValue(line, attribute);
  if (valString == "") return false;
  return boolString(valString);

}

//--------------------------------------------------------------------------

// Extract XML int value following XML attribute.

int Settings::intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;

}

//--------------------------------------------------------------------------

// Extract XML double value following XML attribute.

double Settings::doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;

}

//--------------------------------------------------------------------------

// Extract XML bool vector value following XML attribute.

vector<bool> Settings::boolVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<bool>();
  vector<bool> vectorVal;
  size_t       stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    vectorVal.push_back(boolString(valStream.str()));
  }
  return vectorVal;

}

//--------------------------------------------------------------------------

// Extract XML int vector value following XML attribute.

vector<int> Settings::intVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<int>();
  int         intVal;
  vector<int> vectorVal;
  size_t      stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> intVal;
    vectorVal.push_back(intVal);
  }
  return vectorVal;

}

//--------------------------------------------------------------------------

// Extract XML double vector value following XML attribute.

vector<double> Settings::doubleVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<double>();
  double         doubleVal;
  vector<double> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    istringstream  valStream(valString.substr(0, stringPos));
    valString = valString.substr(stringPos + 1);
    valStream >> doubleVal;
    vectorVal.push_back(doubleVal);
  }
  return vectorVal;

}

//--------------------------------------------------------------------------

// Extract XML string vector value following XML attribute.

vector<string> Settings::stringVectorAttributeValue(string line,
  string attribute) {
  string valString = attributeValue(line, attribute);
  size_t openBrace  = valString.find_first_of("{");
  size_t closeBrace = valString.find_last_of("}");
  if (openBrace != string::npos)
    valString = valString.substr(openBrace + 1, closeBrace - openBrace - 1);
  if (valString == "") return vector<string>();
  string         stringVal;
  vector<string> vectorVal;
  size_t         stringPos(0);
  while (stringPos != string::npos) {
    stringPos = valString.find(",");
    if (stringPos != string::npos) {
      vectorVal.push_back(valString.substr(0, stringPos));
      valString = valString.substr(stringPos + 1);
    } else vectorVal.push_back(valString);
  }
  return vectorVal;

}

//==========================================================================

} // end namespace Pythia8
