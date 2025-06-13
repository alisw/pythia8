#ifndef Pythia8_main144Dct_H
#define Pythia8_main144Dct_H

#include <vector>
#include "TObject.h"
#include "TTree.h"

// Forward declare Pythia objects (here, no implementation is required).
namespace Pythia8 {
  class Particle;
  class Info;
}

//==========================================================================

// RootParticle class.

class RootParticle : public TObject {

public:

  // Data members. Any additional members should be added here.
  double phi{0}, eta{0}, y{0}, pT{0};
  int pid{0};

  // Default constructor is necessary for ROOT.
  RootParticle() = default;

  // Constructor from a Pythia8 Particle. This constructor is not
  // available in the ROOT interpreter, or when linking against this
  // library. The definition for this constructor is in main144.cc and
  // can be modified there
  RootParticle(Pythia8::Particle &prt);

  // Macro to declare this class to ROOT.
  // https://root.cern/manual/io_custom_classes/#the-classdef-macro
  ClassDef(RootParticle, 1)

};

//==========================================================================

// RootEvent class, which contains a RootParticle vector and event Info.

class RootEvent : public TObject {

public:

  // Data members. Any additional members should be added here.
  double weight{1};
  std::vector<RootParticle> particles;

  // Default constructor is necessary for ROOT.
  RootEvent() = default;

  // Fill the event. This method is not available in the ROOT
  // interpreter, or when linking against this library. The definition
  // for this method is in main144.cc and can be modified there.
  void fill(const Pythia8::Info &infoIn, std::vector<RootParticle> &prtsIn,
    TTree *treeIn);

  // Macro to declare this class to ROOT.
  // https://root.cern/manual/io_custom_classes/#the-classdef-macro
  ClassDef(RootEvent, 1)

};

//==========================================================================

// The following pre-processor statements allow ROOT to generate a
// streamer library so the classes above can be stored in ROOT files.
// The following __CLING__ statemement should be changed to __CINT__
// if using __CINT__ to generate the library with ROOT 5.

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class RootParticle+;
#pragma link C++ class RootEvent+;
#endif

#endif // Pythia8_main144Dct_H
