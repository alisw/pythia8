// PartonVertex.h is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for vertex information during the parton-level evolution.

#ifndef Pythia8_PartonVertex_H
#define Pythia8_PartonVertex_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/PhysicsBase.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The PartonVertex class sets parton-level vertex information.

class PartonVertex : public PhysicsBase {

public:

  // Constructor.
  PartonVertex() : doVertex(), modeVertex(), rProton(), pTmin(),
    widthEmission(), bScale(), bNow(), xMax(), yMax(), mux() {}

  // Destructor.
  virtual ~PartonVertex() {}

  // Initialize a few parameters from Settings.
  virtual void init();

  // Select vertex for a Beam Remnant particle.
  virtual void vertexBeam( int iNow, int iBeam, Event& event);

  // Select vertex for an MPI.
  virtual void vertexMPI( int iBeg, int nAdd, double bNowIn, Event& event);

  // Select vertex for an FSR branching.
  virtual void vertexFSR( int iNow, Event& event);

  // Select vertex for an ISR branching.
  virtual void vertexISR( int iNow, Event& event);

private:

  // Data related to currently implemented models.
  bool   doVertex;
  int    modeVertex;
  double rProton, pTmin, widthEmission, bScale;

  // Current values.
  double bNow, xMax, yMax, mux;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PartonVertex_H
