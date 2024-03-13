// HIBasics.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the EventInfo and
// HIUserHooks classes, as well as the HIUnits namespace.
//
// EventInfo: stores full nucleon-nucleon events with corresponding Info.
// HIUserHooks: User hooks for HeavyIons models.

#ifndef Pythia8_HIBasics_H
#define Pythia8_HIBasics_H

#include "Pythia8/Pythia.h"

namespace Pythia8 {

// Forward declarations.
class Nucleon;
class SubCollision;

//==========================================================================

// The HIUnits namespace defines the unitsystem used by the heavy ion
// machinery in Pythia8. In particular all lengths are in femtometer
// and cross sections are therefore in squared femtometer.
namespace HIUnits {

// Lengths
const double femtometer = 1.0;
const double millimeter = 1.0e12;

// Cross sections
const double femtometer2 = 1.0;
const double millibarn = 0.1;
const double nanobarn = 1.0e-7;

}

using namespace HIUnits;

//==========================================================================

// Class for storing Events and Info objects.

class EventInfo {

public:

  // Empty constructor.
  EventInfo(): code(0), ordering(-1.0), coll(0), ok(false) {}

  // The Event object.
  Event event;

  // The corresponding Info object.
  Info info;

  // The code for the subprocess.
  int code;

  // The ordering variable of this event.
  double ordering;
  bool operator<(const EventInfo & ei) const {
    return ordering < ei.ordering;
  }

  // The associated SubCollision object.
  const SubCollision* coll;

  // Is the event properly generated?
  bool ok;

  // Which projectile and target nucleons are included and where are
  // they placed?
  map<Nucleon*, pair<int,int> > projs, targs;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HIBasics_H
