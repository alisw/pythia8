// FragmentationModel.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the classes involved in the modelling of
// fragmentation. A FragmentationModel has one main methd, fragment,
// which takes the color subsystem, color configuration, and
// event. Optional flags can be passed if the event is diffractive and
// if system recoil should be used in some cases. The class hold
// pointers to flavor, pT, and z generators, although these do not
// need to be used by the model. The LundFragmentation model is
// defined in StringFragmentation, and provides the standard Pythia
// fragmentation.

#ifndef Pythia8_FragmentationModel_H
#define Pythia8_FragmentationModel_H

#include "Pythia8/PhysicsBase.h"
#include "Pythia8/FragmentationSystems.h"

namespace Pythia8 {

//==========================================================================

// FragmentationModel is the base class for handling fragmentation algorithms.

class FragmentationModel : public PhysicsBase {

public:

  // Empty constructor.
  FragmentationModel() = default;

  // Empty virtual destructor.
  virtual ~FragmentationModel() {}

  // Initialize and save pointers.
  virtual bool init(StringFlav* flavSelPtrIn = nullptr,
    StringPT* pTSelPtrIn = nullptr, StringZ* zSelPtrIn = nullptr,
    FragModPtr fragModPtrIn = nullptr) = 0;

  // Do the fragmentation: driver routine.
  virtual bool fragment(int iSub, ColConfig& colConfig, Event& event,
    bool isDiff = false, bool systemRecoil = true) = 0;

protected:

  // Pointers to classes for flavour, pT and z generation.
  StringFlav*   flavSelPtr{};
  StringPT*     pTSelPtr{};
  StringZ*      zSelPtr{};

};

//==========================================================================

// Forward reference to StringFragmentation and
// MiniStringFragmentation classes; needed in LundFragmentation class.
class StringFragmentation;
class MiniStringFragmentation;

//--------------------------------------------------------------------------

// The LundFragmentation class handles the default Pythia fragmentation,
// using both the StringFragmentation and MiniStringFragmentation classes.

class LundFragmentation : public FragmentationModel {

public:

  // Constructor (creates string and mini-string pointers, needed for
  // forward declaration and factorization).
  LundFragmentation();

  // Destructor (deletes string and mini-string pointers).
  ~LundFragmentation() override;

  // Initialize and save pointers.
  bool init(StringFlav* flavSelPtrIn = nullptr,
    StringPT* pTSelPtrIn = nullptr, StringZ* zSelPtrIn = nullptr,
    FragModPtr fragModPtrIn = nullptr) override;

  // Do the fragmentation: driver routine.
  bool fragment(int iSub, ColConfig& colConfig, Event& event,
    bool isDiff = false, bool systemRecoil = true) override;

  // Internal StringFragmentation and MiniStringFragmentation objects.
  StringFragmentation* stringFragPtr{};
  MiniStringFragmentation* ministringFragPtr{};

private:

  // Parameters controlling the fragmentation.
  double mStringMin{};
  bool tryMiniAfterFailedFrag{};

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_FragmentationModel_H
