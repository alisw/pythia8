// Basics.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Logger class.

#include "Pythia8/Logger.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

void Logger::init(Settings& settings) {
  isQuietSave    = settings.flag("Print:quiet");
  printNextSave  = settings.flag("Print:next");
  printInitSave  = settings.flag("Print:init");
  printErrors    = settings.flag("Print:errors");
  verbosity      = settings.mode("Print:verbosity");
  useErrorStream = settings.flag("Print:useErrorStream");
}

//--------------------------------------------------------------------------

void Logger::reportMsg(string loc, string m, string extra, bool showAlways) {
  msg(REPORT, "Report from " + loc + ": " + m, extra, showAlways); }

void Logger::infoMsg(string loc, string m, string extra, bool showAlways) {
  msg(NORMAL, "Info from " + loc + ": " + m, extra, showAlways); }

void Logger::warningMsg(string loc, string m, string extra, bool showAlways) {
  msg(NORMAL, "Warning in " + loc + ": " + m, extra, showAlways); }

void Logger::errorMsg(string loc, string m, string extra, bool showAlways) {
  msg(NORMAL, "Error in " + loc + ": " + m, extra, showAlways); }

void Logger::abortMsg(string loc, string m, string extra, bool showAlways) {
  msg(ABORT, "Abort from " + loc + ": " + m, extra, showAlways); }

//--------------------------------------------------------------------------

void Logger::msg(int verbosityLevel, string message, string extraInfo,
  bool showAlways) {

  // Ignore messages that don't satisfy the verbosity level.
  if (verbosity < verbosityLevel)
    return;

  // Only one thread can write messages at the same time
  lock_guard<mutex> lock(writeMutex);

  // Update message counter.
  int times = messages[message]++;

  // Print message.
  if ( (times == 0 || showAlways) && mayPrintErrors() ) {
    string messageNow = " PYTHIA " + message;
    if (extraInfo != "") messageNow += " " + extraInfo;
    errorStream() << messageNow + "\n";
  }
}

//--------------------------------------------------------------------------

// Add all errors from the other Logger object to the counts of this object.

void Logger::errorCombine(const Logger& other, string prefix) {
  for (pair<string, int> messageEntry : other.messages) {
    string name = messageEntry.first;
    if ( !prefix.empty() ) {
      if ( name.substr(0, 5) == "Info " )
        name = "Info " + prefix + " " + name.substr(5);
      else if ( name.substr(0, 8) == "Warning " )
        name = "Warning " + prefix + " " + name.substr(8);
      else if ( name.substr(0, 6) == "Error " )
        name = "Error " + prefix + " " + name.substr(6);
      else if ( name.substr(0, 6) == "Abort " )
        name = "Abort " + prefix + " " + name.substr(6);
      else
        name = prefix + " " + name;
    }
    messages[name] += messageEntry.second;
  }
}

//--------------------------------------------------------------------------

void Logger::errorStatistics(ostream& stream) const {
  // Header.
  stream << "\n *-------  PYTHIA Error and Warning Messages Statistics  "
            "----------------------------------------------------------* \n"
            " |                                                       "
            "                                                          | \n"
            " |  times   message                                      "
            "                                                          | \n"
            " |                                                       "
            "                                                          | \n";

  // Loop over all messages
  map<string, int>::const_iterator messageEntry = messages.begin();
  if (messageEntry == messages.end())
    stream << " |      0   no errors or warnings to report              "
           << "                                                          | \n";
  while (messageEntry != messages.end()) {
    // Message printout.
    string temp = messageEntry->first;
    int len = temp.length();
    temp.insert( len, max(0, 102 - len), ' ');
    stream << " | " << setw(6) << messageEntry->second << "   "
           << temp << " | \n";
    ++messageEntry;
  }

  // Done.
  stream << " |                                                       "
            "                                                          | \n"
            " *-------  End PYTHIA Error and Warning Messages Statistics"
            "  ------------------------------------------------------* "
         << endl;
}

//==========================================================================

} // end namespace Pythia8
