// FragmentationModel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// FragmentationModel (currently none) and LundFragmentation classes.

#include "Pythia8/StringFragmentation.h"
#include "Pythia8/MiniStringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The LundFragmentation class.

//--------------------------------------------------------------------------

// Constructor.

LundFragmentation::LundFragmentation() {

  stringFragPtr     = new StringFragmentation();
  ministringFragPtr = new MiniStringFragmentation();

}

//--------------------------------------------------------------------------

// Destructor.

LundFragmentation::~LundFragmentation() {

  delete stringFragPtr;
  delete ministringFragPtr;

}

//--------------------------------------------------------------------------

// Initialize and save pointers.

bool LundFragmentation::init(StringFlav* flavSelPtrIn, StringPT* pTSelPtrIn,
  StringZ* zSelPtrIn, FragModPtr fragModPtrIn) {

  // Register the string and ministring fragmentation models.
  registerSubObject(*stringFragPtr);
  registerSubObject(*ministringFragPtr);

  // Set the pointers.
  stringFragPtr->init(flavSelPtrIn, pTSelPtrIn, zSelPtrIn, fragModPtrIn);
  ministringFragPtr->init(flavSelPtrIn, pTSelPtrIn, zSelPtrIn, fragModPtrIn);

  // Boundary mass between string and ministring handling.
  mStringMin = parm("HadronLevel:mStringMin");

  // Try ministring fragmentation also if normal fails.
  tryMiniAfterFailedFrag = flag("MiniStringFragmentation:tryAfterFailedFrag");

  // Return successful.
  return true;

}

//--------------------------------------------------------------------------

// Do the fragmentation: driver routine.

bool LundFragmentation::fragment(int iSub, ColConfig& colConfig, Event& event,
  bool isDiff, bool) {

  // String fragmentation of each colour singlet (sub)system.
  // If fails optionally try ministring fragmentation.
  if (iSub == -1) return true;
  if (colConfig[iSub].massExcess > mStringMin) {
    if (!stringFragPtr->fragment( iSub, colConfig, event)) {
      if (!tryMiniAfterFailedFrag) return false;
      loggerPtr->ERROR_MSG("string fragmentation failed, "
        "trying ministring fragmetation instead");
      if (!ministringFragPtr->fragment(iSub, colConfig, event, isDiff)) {
        loggerPtr->ERROR_MSG("also ministring fragmentation failed "
          "after failed normal fragmentation");
        return false;
      }
    }

  // Low-mass string treated separately.
  } else {
    if (!ministringFragPtr->fragment( iSub, colConfig, event, isDiff)) {
      loggerPtr->ERROR_MSG("ministring fragmentation failed");
      return false;
    }
  }

  // Return successful.
  return true;

}

//==========================================================================

} // end namespace Pythia8
