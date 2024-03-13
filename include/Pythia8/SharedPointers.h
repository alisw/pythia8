// SharedPointers.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header files for forward declaring of some common Pythia classes
// and typedefs for shared pointers to these.

#ifndef Pythia8_SharedPointers_H
#define Pythia8_SharedPointers_H

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

class BeamShape;
typedef shared_ptr<BeamShape> BeamShapePtr;

class Brancher;
typedef shared_ptr<Brancher> BrancherPtr;

class ColourDipole;
typedef shared_ptr<ColourDipole> ColourDipolePtr;

class ColourReconnectionBase;
typedef shared_ptr<ColourReconnectionBase> ColRecPtr;

class DecayHandler;
typedef shared_ptr<DecayHandler> DecayHandlerPtr;

class DipoleSwingBase;
typedef shared_ptr<DipoleSwingBase> DipSwingPtr;

class ExternalMEs;
typedef shared_ptr<ExternalMEs> ExternalMEsPtr;

class FragmentationModifierBase;
typedef shared_ptr<FragmentationModifierBase> FragModPtr;

class HeavyIons;
typedef shared_ptr<HeavyIons> HeavyIonsPtr;

class HIUserHooks;
typedef shared_ptr<HIUserHooks> HIUserHooksPtr;

class LHAup;
typedef shared_ptr<LHAup> LHAupPtr;

class LHEF3FromPythia8;
typedef shared_ptr<LHEF3FromPythia8> LHEF3FromPythia8Ptr;

class Merging;
typedef shared_ptr<Merging> MergingPtr;

class MergingHooks;
typedef shared_ptr<MergingHooks> MergingHooksPtr;

class PartonVertex;
typedef shared_ptr<PartonVertex> PartonVertexPtr;

class ParticleDataEntry;
typedef shared_ptr<ParticleDataEntry> ParticleDataEntryPtr;

class PDF;
typedef shared_ptr<PDF> PDFPtr;

class PhaseSpace;
typedef shared_ptr<PhaseSpace> PhaseSpacePtr;

class ResonanceWidths;
typedef shared_ptr<ResonanceWidths> ResonanceWidthsPtr;

class RndmEngine;
typedef shared_ptr<RndmEngine> RndmEnginePtr;

class ShowerModel;
typedef shared_ptr<ShowerModel> ShowerModelPtr;

class SigmaProcess;
typedef shared_ptr<SigmaProcess> SigmaProcessPtr;

class SpaceShower;
typedef shared_ptr<SpaceShower> SpaceShowerPtr;

class SplitOnia;
typedef shared_ptr<SplitOnia> SplitOniaPtr;

class StringInteractions;
typedef shared_ptr<StringInteractions> StringIntPtr;

class StringRepulsionBase;
typedef shared_ptr<StringRepulsionBase> StringRepPtr;

class TimeShower;
typedef shared_ptr<TimeShower> TimeShowerPtr;

class UserHooks;
typedef shared_ptr<UserHooks> UserHooksPtr;

class VinciaModule;
typedef shared_ptr<VinciaModule> VinciaModulePtr;

class ZetaGenerator;
typedef shared_ptr<ZetaGenerator> ZetaGeneratorPtr;

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SharedPointers_H
