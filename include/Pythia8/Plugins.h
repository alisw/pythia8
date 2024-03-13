// Plugins.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Philip Ilten, Manuel Szewc, and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the runtime loading of plugins.

#ifndef Pythia8_Plugins_H
#define Pythia8_Plugins_H

#include "Pythia8/Pythia.h"

// Allow the use of dlopen without warnings for some GCC versions.
#if defined (__GNUC__) && ((__GNUC__ + 0) < 5)
#pragma GCC system_header
#endif

namespace Pythia8 {

//==========================================================================

// Demangle a symbol name, if the necessary demangling libraries are present.

string demangle(string name);

//==========================================================================

// Determine a plugin type.

string type_plugin(string libName, string className,
  Logger* loggerPtr = nullptr);

//==========================================================================

// Load a plugin library with dlopen.

shared_ptr<void> dlopen_plugin(string libName, Logger* loggerPtr);

//==========================================================================

// Load a symbol from a plugin library.

template <typename T> function<T> dlsym_plugin(void* libPtr, string symbol) {
  return (T*)dlsym(libPtr, symbol.c_str());}

template <typename T> function<T> dlsym_plugin(shared_ptr<void> libPtr,
  string symbol) {
  return (T*)dlsym( static_cast<void*>(libPtr.get()), symbol.c_str());}

//==========================================================================

// Load a plugin, given a full set of arguments.

template <typename T> shared_ptr<T> make_plugin(
  string libName, string className, Pythia* pythiaPtr,
  Settings* settingsPtr, Logger* loggerPtr) {

  // Set up the available pointers.
  if (loggerPtr == nullptr && pythiaPtr != nullptr)
    loggerPtr = &pythiaPtr->logger;
  if (settingsPtr == nullptr && pythiaPtr != nullptr)
    settingsPtr = &pythiaPtr->settings;

  // Load the library.
  shared_ptr<void> libPtr = dlopen_plugin(libName, loggerPtr);
  if (libPtr == nullptr) return shared_ptr<T>(nullptr);

  // Check the plugin object type.
  string objType = type_plugin(libName, className, loggerPtr);
  if (objType != typeid(T).name()) {
    string msg = "class " + className + " from library " + libName +
      " must be loaded as type " + demangle(objType);
    if (loggerPtr != nullptr) loggerPtr->errorMsg("make_plugin", msg);
    else cout << msg << "\n";
    return shared_ptr<T>(nullptr);
  }

  // Check the pointers required by the plugin object.
  for (string ptr : {"PYTHIA", "SETTINGS", "LOGGER"}) {
    auto ptrsObject =
      dlsym_plugin<bool()>(libPtr, "REQUIRE_" + ptr + "_" + className);
    if (dlerror() != nullptr) continue;
    if (!ptrsObject()) continue;
    if (ptr == "PYTHIA" && pythiaPtr != nullptr) continue;
    if (ptr == "SETTINGS" && settingsPtr != nullptr) continue;
    if (ptr == "LOGGER" && loggerPtr != nullptr) continue;
    string msg = "class " + className + " requires a " + ptr + " pointer";
    if (loggerPtr != nullptr) loggerPtr->errorMsg("make_plugin", msg);
    else cout << msg << "\n";
    return shared_ptr<T>(nullptr);
  }

  // Load the symbol to construct the plugin object.
  auto newObject =
    dlsym_plugin<T*(Pythia*, Settings*, Logger*)>(libPtr, "NEW_" + className);
  const char* error = dlerror();
  if (error != nullptr) {
    string msg = "class " + className + " not available from library " +
      libName;
    if (loggerPtr != nullptr) loggerPtr->errorMsg("make_plugin", msg);
    else cout << msg << "\n";
    return shared_ptr<T>(nullptr);
  }

  // Construct the plugin object shared pointer.
  return shared_ptr<T>(
    newObject(pythiaPtr, settingsPtr, loggerPtr),
    // Use a custom destructor.
    [libPtr, className](T* objectPtr) {

      // Destroy the plugin object.
      auto deleteObject =
        dlsym_plugin<void(T*)>(libPtr, "DELETE_" + className);
      if (dlerror() == nullptr && deleteObject != nullptr)
        deleteObject(objectPtr);});

}

//==========================================================================

// Load a plugin, given no pointers.

template <typename T> shared_ptr<T> make_plugin(
  string libName, string className) {return make_plugin<T>(
    libName, className, nullptr, nullptr, nullptr);}

//==========================================================================

// Load a plugin, given a Pythia pointer.
template <typename T> shared_ptr<T> make_plugin(
  string libName, string className, Pythia* pythiaPtr) {
  return make_plugin<T>(
    libName, className, pythiaPtr, nullptr, nullptr);}

//==========================================================================

// Load a plugin, given a Pythia pointer and command vector pointer.

template <typename T> shared_ptr<T> make_plugin(
  string libName, string className, Pythia* pythiaPtr,
  const vector<string>& cmnds) {
  pythiaPtr->settings.registerPluginLibrary(libName);
  for (string cmnd : cmnds) pythiaPtr->readString(cmnd);
  return make_plugin<T>(libName, className, pythiaPtr);
}

//==========================================================================

// Load a plugin, given a Pythia pointer, command file, and subrun.

template <typename T> shared_ptr<T> make_plugin(
  string libName, string className, Pythia* pythiaPtr,
  string fileName, int subrun = SUBRUNDEFAULT) {
  pythiaPtr->settings.registerPluginLibrary(libName);
  if (fileName != "") pythiaPtr->readFile(fileName, subrun);
  return make_plugin<T>(libName, className, pythiaPtr);
}

//==========================================================================

// Macro to declare a plugin class.

#define PYTHIA8_PLUGIN_CLASS(BASE, CLASS, PYTHIA, SETTINGS, LOGGER) \
  extern "C" {                                                      \
    bool REQUIRE_PYTHIA_##CLASS() {return PYTHIA;}                  \
    bool REQUIRE_SETTINGS_##CLASS() {return SETTINGS;}              \
    bool REQUIRE_LOGGER_##CLASS() {return LOGGER;}                  \
    const char* TYPE_##CLASS() {return typeid(BASE).name();}        \
    CLASS* NEW_##CLASS(Pythia* pythiaPtr, Settings* settingsPtr,    \
      Logger* loggerPtr) {                                          \
      return new CLASS(pythiaPtr, settingsPtr, loggerPtr);}         \
    void DELETE_##CLASS(CLASS* ptr) {delete ptr;}}

//==========================================================================

// Macro to register settings for a plugin library.

#define PYTHIA8_PLUGIN_SETTINGS(METHOD) \
  extern "C" {                          \
    void REGISTER_SETTINGS(Settings* settingsPtr) {METHOD(settingsPtr);}}

//==========================================================================

// Macro to register an XML settings index file.

#define PYTHIA8_PLUGIN_XML(INDEX) \
  extern "C" {const char* RETURN_XML() {return INDEX;}}

//==========================================================================

// Macro to return compatible Pythia versions and the compiled version.

#define PYTHIA8_PLUGIN_VERSIONS(...)                                       \
  extern "C" {                                                             \
    bool CHECK_COMPATIBLE_VERSION(int ver) {set<int> vers = {__VA_ARGS__}; \
      return vers.find(ver) != vers.end();}                                \
    bool CHECK_COMPILED_VERSION(int ver) {                                 \
      return ver == PYTHIA_VERSION_INTEGER;}}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_Plugins_H
