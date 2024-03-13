// Plugins.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Philip Ilten, Manuel Szewc, and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for plugin manipulation.

#include "Pythia8/Plugins.h"

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {

//==========================================================================

// Demangle a symbol name, if the necessary demangling libraries are present.

// The following method was partially taken from
// boost/core/demangle.hpp
//
// Copyright 2014 Peter Dimov
// Copyright 2014 Andrey Semashev
//
// Distributed under the Boost Software License, Version 1.0.
// See include/Pythia8/PythiaStdLib.h or http://www.boost.org/LICENSE_1_0.txt

#if defined( __has_include ) && ((__GNUC__ + 0) >= 5)
# if __has_include(<cxxabi.h>)
#  define HAS_CXXABI_H
# endif
#elif defined( __GLIBCXX__ ) || defined( __GLIBCPP__ )
# define HAS_CXXABI_H
#endif
#if defined( HAS_CXXABI_H )
# include <cxxabi.h>
# if defined( __GABIXX_CXXABI_H__ )
#  undef HAS_CXXABI_H
# endif
#endif
#if defined( HAS_CXXABI_H )
string demangle(string name) {
  unique_ptr<char, void(*)(void*)> demangle(
    abi::__cxa_demangle(name.c_str(), nullptr, nullptr, nullptr), free);
  return demangle.get();}
#else
string demangle(string name) {return name;}
#endif

//==========================================================================

// Load a plugin library with dlopen. The shared pointer destructor
// calls dlcose automatically.

shared_ptr<void> dlopen_plugin(string libName, Logger* loggerPtr) {

  // Load the plugin library.
  void* libPtr = dlopen(libName.c_str(), RTLD_LAZY);
  const char* error = dlerror();
  if (error != nullptr) {
    if (loggerPtr != nullptr) loggerPtr->ERROR_MSG(string(error));
    else cout << string(error) << "\n";
    return shared_ptr<void>(nullptr);
  }

  // Check the Pythia version is compatible.
  auto versionObject =
    dlsym_plugin<bool(int)>(libPtr, "CHECK_COMPATIBLE_VERSION");
  error = dlerror();
  if (error != nullptr) {
    string msg("could not determine compatible Pythia versions for "
      + libName);
    if (loggerPtr != nullptr) loggerPtr->ERROR_MSG(msg);
    else cout << msg << "\n";
    return shared_ptr<void>(nullptr);
  } else if (!versionObject(PYTHIA_VERSION_INTEGER)) {
    stringstream ver;
    ver << fixed << setprecision(3) << PYTHIA_VERSION;
    string msg(libName + " is not compatible with Pythia version "
      + ver.str());
    if (loggerPtr != nullptr) loggerPtr->ERROR_MSG(msg);
    else cout << msg << "\n";
    return shared_ptr<void>(nullptr);
  }

  // Check the compiled Pythia version is the same (warning, not error).
  auto compiledVersionObject =
    dlsym_plugin<bool(int)>(libPtr, "CHECK_COMPILED_VERSION");
  error = dlerror();
  if (error != nullptr) {
    string msg("could not determine the version of Pythia used when compiling "
      + libName);
    if (loggerPtr != nullptr) loggerPtr->ERROR_MSG(msg);
    else cout << msg << "\n";
    return shared_ptr<void>(nullptr);
  } else if (!compiledVersionObject(PYTHIA_VERSION_INTEGER)) {
    stringstream ver;
    ver << fixed << setprecision(3) << PYTHIA_VERSION;
    string msg(libName + " was not compiled with Pythia version " + ver.str());
    if (loggerPtr != nullptr) loggerPtr->WARNING_MSG(msg);
    else cout << msg << "\n";
  }

  // Return the shared pointer.
  return shared_ptr<void>(
    libPtr,
    // Use a custom destructor.
    [](void* ptr) {
      dlclose(ptr);
      dlerror();});

}

//==========================================================================

// Determine a plugin type.

string type_plugin(string libName, string className, Logger* loggerPtr) {

  // Load the plugin library.
  shared_ptr<void> libPtr = dlopen_plugin(libName, loggerPtr);
  if (libPtr == nullptr) return "";

  // Grab the plugin object type.
  auto typeObject = dlsym_plugin<const char*()>(libPtr, "TYPE_" + className);
  if (dlerror() != nullptr) {
    string msg = "class " + className + " not available from library " +
      libName;
    if (loggerPtr) loggerPtr->ERROR_MSG(msg);
    else cout << msg << "\n";
    return "";
  } else return typeObject();

}

//==========================================================================

} // end namespace Pythia8
