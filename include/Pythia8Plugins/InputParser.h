// InputParser.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Christian Bierlich, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef Pythia8_InputParser_H
#define Pythia8_InputParser_H

// Includes.
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Convenience class for parsing C++ command-line arguments on the form
// ./mainXX -a A -b B etc.
//
// Usage: (1) Create an InputParser object in your main program. (2)
// Specify the command-line options that should be available, including
// default values, using add. (3) Call the parse(int& argc, char** argv)
// method, with the two aruments forwarded from argc and argv from your
// main. (4) Extract given command-line values of type T, using get<T>.

class InputParser {
public:

  // Constructor.
  // usage: Text to print when describing usage with help.
  // examples: Vector of examples to print with help.
  // extra: Extra text to print with help.
  // stream: Optional pointer to stream to print messages to.
  // optName: Option name for the help flag.
  // aliases: Aliases for the help flag.
  InputParser(string usage = "", vector<string> examples = {},
    string extra = "", ostream* stream = &cout,
    string optName = "h", set<string> aliases = {"-help","H"}) :
    usageText(usage), examplesText(examples), extraText(extra),
      streamPtr(stream), helpFlag(optName) {
    add(optName, "false", "Show this help message and exit.", aliases);
  }

  // Add an optional command line option to the parser. Do this before parsing.
  // Parameters:
  // optName: The name of the option.
  // defString: The default value for the option.
  // helpText: The help text describing the option (optional).
  // aliases: A set of aliases for the option (optional).
  bool add(const string& optName, const string& defString,
    const string& helpText = "", set<string> aliases = {}) {
    // Check for name conflicts with existing options or aliases.
    if (options.find(optName) != options.end() ||
      aliasMap.find(optName) != aliasMap.end()) {
      print("Name conflict for '" + optName + "'.\n");
      return false;
    }
    // Create an OptionInfo object and add it to the options map.
    OptionInfo optInfo(optName, defString, helpText, aliases);
    options[optName] = optInfo;
    // Add aliases to the alias map.
    for (const string& alias : aliases) {
      if (options.find(alias) != options.end() ||
        aliasMap.find(alias) != aliasMap.end()) {
        print("Name conflict for alias '" + alias + "'.\n");
        return false;
      }
      aliasMap[alias] = optName;
    }
    return true;
  }

  // Add required command line option to the parser. Do this before parsing.
  bool require(const string& optName, const string& helpText = "",
    set<string> aliases = {}) {
    if (!add(optName, "", helpText, aliases)) return false;
    options[optName].required = true;
    return true;
  }

  // Initialize the parser. This includes parsing, printing help if
  // requested, and checking the required options.
  enum Status {Valid = -1, Invalid = EXIT_FAILURE, Help = EXIT_SUCCESS};
  Status init(int& argc, char** argv) {

    // Parse the arguments.
    bool valid = parse(argc, argv);
    if (!valid) return Status::Invalid;

    // Print the help if requested.
    if (options.find(helpFlag) != options.end() && get<bool>(helpFlag)) {
      string out;
      if (usageText != "") out = "Usage: " + usageText;
      if (examplesText.size() > 0) {
        out += "\nExample" + string(examplesText.size() > 1 ? "s:" : ":");
        for (const string& text: examplesText) out += "\n\t" + text;
      }
      if (options.size() > 0) {
        out += "\nOption" + string(options.size() > 1 ? "s:\n" : ":\n")
          + help() + "\n";
      }
      print(out + extraText);
      return Status::Help;
    }

    // Check the required options are set.
    for (map<string, OptionInfo>::iterator opt = options.begin();
         opt != options.end(); opt++) {
      if (opt->second.required && !opt->second.provided) {
        print("Option '-" + opt->first + "' is required but was not set.\n"
          "\t -" + opt->first + " " + opt->second.helpText + "\n");
        valid = false;
      }
    }
    return valid ? Status::Valid : Status::Invalid;

  }

  // Method to parse command line arguments.
  // Returns true if parsing was successful, false otherwise.
  // Print error messages to stream.
  // The hflag option specifies the help option.
  bool parse(int& argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
      string arg = argv[i];
      // Check if the argument is an option (starts with '-').
      if (arg[0] == '-') {
        string optName = arg.substr(1);
        // Check for aliases and get the actual option name.
        if (aliasMap.find(optName) != aliasMap.end())
          optName = aliasMap[optName];
         // Check if the option is defined.
        if (options.find(optName) != options.end()) {
          // If the next argument is not an option, set the value for
          // this option.
          if (i + 1 < argc && argv[i + 1][0] != '-') {
            options[optName].stringValues.push_back(argv[++i]);
            options[optName].provided = true;
          } else {
            // Treat as boolean flag and flip its value.
            const string& def = options[optName].defaultValue;
            if (def == "false" || def == "False" || def == "0"
              || def == "off") {
              options[optName].stringValues.push_back("true");
              options[optName].provided = true;
            } else if (def == "true" || def == "True" || def == "1"
              || def == "on") {
              options[optName].stringValues.push_back("false");
              options[optName].provided = true;
            } else {
              print("Failed to parse command line arguments.\n"
                "No value passed for option '" + string(arg) + "'.\n");
              return false;
            }
          }
        } else {
          print("Failed to parse command line arguments.\n"
            "Unknown option '" + string(arg) + "'.\n");
          return false;
        }
      }
    }
    return true;
  }

  // Check if an option is defined.
  // Returns true if the option is defined, false otherwise.
  bool has(const string& optName) const {
    return options.find(optName) != options.end();
  }

  // Templated method to get the last value of an option.
  // Returns the last value of the option converted to the specified type.
  template<typename T>
  T get(const string& optName) {
    if (!has(optName)) {
      print("Failed to find option '" + optName + "'.\n");
      return T();
    }
    // Retrieve the OptionInfo object for the given option name.
    const OptionInfo& optInfo = options.at(optName);
    // Return default-constructed T if string is empty.
    if (optInfo.stringValues.empty()) return T();
    if (optInfo.stringValues.back().empty()) return T();
    // Convert the string value to the specified type.
    stringstream conv(optInfo.stringValues.back());
    T value;
    conv >> std::boolalpha >> value;
    // Error message and default constructed T() is conversion failed.
    if (conv.fail()) {
      print("Failed to convert '" + optInfo.optName + "'.\n");
      return T();
    }
    return value;
  }

  // Templated method to get all the values of an option.
  template<typename T>
  vector<T> getVector(const string& optName) {
    vector<T> values;
    if (!has(optName)) {
      print("Failed to find option '" + optName + "'.\n");
      return values;
    }
    // Retrieve the OptionInfo object for the given option name.
    const OptionInfo& optInfo = options.at(optName);
    // Return empty vector if no values.
    if (optInfo.stringValues.empty()) return values;
    // Convert the string value to the specified type.
    for (const string& stringValue: optInfo.stringValues) {
      if (stringValue.empty()) {
        values.push_back(T());
        continue;
      }
      stringstream conv(stringValue);
      T value;
      conv >> std::boolalpha >> value;
      // Error message and default constructed T() is conversion failed.
      if (conv.fail()) {
        print("Failed to convert '" + optInfo.optName + "'.\n");
        return values;
      }
      values.push_back(value);
    }
    return values;
  }

  // Method to generate the help text for all options.
  // Returns a formatted string containing the help text.
  const string help() const {
    stringstream out;
    auto oItr = options.cbegin();
    while (oItr != options.cend()) {
      out << "\t-" << oItr->second.optName;
      if (!oItr->second.aliases.empty()) {
        out << " (";
        auto aItr = oItr->second.aliases.cbegin();
        while(aItr != oItr->second.aliases.cend()) {
          out << "-" << *aItr;
          if (++aItr != oItr->second.aliases.cend()) out << ", ";
        }
        out << ")";
      }
      else out << "\t";
      out << "\t" << oItr->second.helpText;
      if (oItr->second.required)
        out << " (required)";
      else if (oItr->second.defaultValue != "")
        out << " (default: " << oItr->second.defaultValue << ")";
      if (++oItr != options.cend()) out << "\n";
    }
    return out.str();
  }

private:

  // Struct to hold information about each option.
  struct OptionInfo {
    // Default constructor must exist.
    OptionInfo() {}
    // Constructor to initialize all fields.
    OptionInfo(const string& n, const string& v, const string& h,
      set<string> a) : optName(n), defaultValue(v),
      stringValues(1, v), helpText(h), aliases(a) {}
    string optName;              // Name of the option.
    string defaultValue;         // Default value of the option.
    vector<string> stringValues; // Values of the option.
    string helpText;             // Help text for the option.
    set<string> aliases;         // Aliases for the option.
    bool required{false};        // Flag if this option must be provided.
    bool provided{false};        // Flag that the user provided this option.
  };

  // Members.
  string usageText{};            // Text for usage.
  vector<string> examplesText{}; // Text for examples.
  string extraText{};            // Extra text to print.
  ostream* streamPtr{};          // Optional stream to print messages.
  string helpFlag{};             // The flag used to access help.

  // Method to print to stream.
  void print(string out) {if (streamPtr != nullptr) *streamPtr << out;}

  // Map to store options with their names as keys.
  map<string, OptionInfo> options;
  // Map to store aliases with aliases as keys and option names as values.
  map<string, string> aliasMap;

};

} // end namespace Pythia8

#endif // end Pythia8_InputParser_H
