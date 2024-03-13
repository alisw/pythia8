// Logger.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the Logger class for diagnostic messages.

#ifndef Pythia8_Logger_H
#define Pythia8_Logger_H

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

class Settings;

//==========================================================================

// The Logger class prints diagnostic messages and counts error messages.
// Objects of the Logger class can be used as ostream object,
// supporting the << operator.

class Logger : private streambuf, public ostream {

  // This struct defines the ordering of log messages.
  // If the two messages have different severity (abort, error, warning,
  // info, debug), order by most severe first. Otherwise, order alphabetically.
  struct LogComparer {
    bool operator() (const string& a, const string& b) const {
      char a1 = a[0], b1 = b[0];
      int severityA, severityB;
      if      (a1 == 'A') severityA = 0;
      else if (a1 == 'E') severityA = 1;
      else if (a1 == 'W') severityA = 2;
      else if (a1 == 'I') severityA = 3;
      else if (a1 == 'R') severityA = 4;
      else severityA = 5;
      if      (b1 == 'A') severityB = 0;
      else if (b1 == 'E') severityB = 1;
      else if (b1 == 'W') severityB = 2;
      else if (b1 == 'I') severityB = 3;
      else if (b1 == 'R') severityB = 4;
      else severityB = 5;
      if (severityA != severityB) return severityA < severityB;
      return a < b;
    }
  };

public:

  // Construct with default values.
  Logger() : ostream(this), infoStreamSave(cout), errStreamSave(cerr),
    verbosity(2), printInitSave(true), printNextSave(true), printErrors(true),
    isQuietSave(false), useErrorStream(false) { }
  virtual ~Logger() {}

  void init(Settings& settings);

  // Methods for providing different levels of error messaging.
  // These methods are thread safe.

  // Report messages contain information not relevant in normal runs.
  void reportMsg(string loc, string message, string extraInfo = "",
    bool showAlways = false);

  // Info messages are diagnostic messages that don't indicate an issue.
  void infoMsg(string loc, string message, string extraInfo = "",
    bool showAlways = false);

  // Warnings indicate that there might be an issue, but the run will continue.
  void warningMsg(string loc, string message, string extraInfo = "",
    bool showAlways = false);

  // Errors indicate an issue that might cause the current event to fail.
  void errorMsg(string loc, string message, string extraInfo = "",
    bool showAlways = false);

  // Aborts indicate critical issues that prevent further event generation.
  void abortMsg(string loc, string message, string extraInfo = "",
    bool showAlways = false);

  // If quiet, the logger will not print any messages.
  bool isQuiet()      const { return isQuietSave; }

  // Indicates whether information will be printed during initialization.
  bool mayPrintInit() const { return printInitSave && !isQuietSave; }

  // Indicates whether information will be printed during event generation.
  bool mayPrintNext() const { return printNextSave && !isQuietSave; }

  // Indicates whether error messages will be printed.
  bool mayPrintErrors() const { return printErrors && !isQuietSave; }

  // 0: no messages | 1: aborts only | 2: default | 3: debug
  void setVerbosity(int verbosityIn) { verbosity = verbosityIn; }
  int getVerbosity() const { return verbosity; }

  // Add all errors from the other Logger object to the counts of this object.
  void errorCombine(const Logger& other, string prefix = "");

  // Reset to empty map of error messages.
  void errorReset() {messages.clear();}

  // Total number of errors/aborts/warnings logged.
  int errorTotalNumber() const {
    int nTot = 0;
    for (pair<string, int> messageEntry : messages)
      nTot += messageEntry.second;
    return nTot;
  }

  // Print error statistics.
  void errorStatistics() const { errorStatistics(infoStreamSave); }
  void errorStatistics(ostream& stream) const;

  ostream& infoStream() { return infoStreamSave; }
  ostream& errorStream() {
    return useErrorStream ? errStreamSave : infoStreamSave; }

  // Iterators over error messages.
  map<string, int, LogComparer>::iterator begin() {
    return messages.begin(); }
  map<string, int, LogComparer>::iterator end() {
    return messages.end(); }
  map<string, int, LogComparer>::const_iterator begin() const {
    return messages.begin(); }
  map<string, int, LogComparer>::const_iterator end() const {
    return messages.end(); }

  // Abort verbosity level: print only abort messages.
  static constexpr int ABORT  = 1;
  // Normal verbosity level: print errors, warnings and info messages.
  static constexpr int NORMAL = 2;
  // Report verbosity level: print everything, including report messages.
  static constexpr int REPORT = 3;

private:

  // Map for all error messages.
  map<string, int, LogComparer> messages;

  // Override method from std::streambuf.
  // This makes the operator<< write to infoStreamSave.
  int overflow(int c) override {
    infoStreamSave.put(c);
    return 0;
  }

  // Streams where messages are written.
  ostream& infoStreamSave;
  // Optional separate error stream.
  ostream& errStreamSave;

  // Configuration.
  int verbosity;
  bool printInitSave, printNextSave, printErrors, isQuietSave, useErrorStream;

  // Mutual exclusive access to writing messages.
  mutex writeMutex;

  // Method to print the diagnostic message. This method is thread safe.
  void msg(int verbosityLevel, string message, string extraInfo = "",
    bool showAlways = false);

};

//--------------------------------------------------------------------------

// Shorthand macros for logger messages.

#define INFO_MSG(...) infoMsg(__METHOD_NAME__, __VA_ARGS__)
#define WARNING_MSG(...) warningMsg(__METHOD_NAME__, __VA_ARGS__)
#define ERROR_MSG(...) errorMsg(__METHOD_NAME__, __VA_ARGS__)
#define ABORT_MSG(...) abortMsg(__METHOD_NAME__, __VA_ARGS__)
#define REPORT_MSG(...) reportMsg(__METHOD_NAME__, __VA_ARGS__)

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Logger_H
