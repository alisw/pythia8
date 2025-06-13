// PythiaBatch.h is a part of the PYTHIA event generator.
// Copyright (C) 2025 Aartem Havryliuk and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file to provide batch events in Python via the Awkward package.

#ifndef awkward_PythiaEvent_H
#define awkward_PythiaEvent_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <limits> 
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "awkward/LayoutBuilder.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"

// Define alias declarations for common Awkward templates used
// (typedefs are not used because they do not allow for templates).
using UserDefinedMap = std::map<std::size_t, std::string>;
template<class... BUILDERS> using RecordBuilder =
  awkward::LayoutBuilder::Record<UserDefinedMap, BUILDERS...>;
template<std::size_t field_name, class BUILDER> using RecordField =
  awkward::LayoutBuilder::Field<field_name, BUILDER>;
template<class PRIMITIVE, class BUILDER> using ListOffsetBuilder =
  awkward::LayoutBuilder::ListOffset<PRIMITIVE, BUILDER>;
template<class PRIMITIVE> using NumpyBuilder =
  awkward::LayoutBuilder::Numpy<PRIMITIVE>;
template<class INDEX, class BUILDER> using IndexedOptionBuilder =
  awkward::LayoutBuilder::IndexedOption<INDEX, BUILDER>;

// Define the ordering of the relevant class fields.
enum VecField : size_t {px, py, pz, e};
enum PrtField : size_t {id, status, mother1, mother2, daughter1, daughter2,
    col, acol,p, m, scale, pol, vProd, tau};
enum InfoField : size_t {id1, id2, x1, x2, pdf1, pdf2, alphaS, alphaEM, Q2Fac,
    Q2Ren, mHat, sHat, tHat, uHat, pT2Hat, weights};
enum EventField : size_t {prt, info};

// Define alias declarations for the corresponding builders.
using VecBuilder = RecordBuilder<
  RecordField<VecField::px, NumpyBuilder<double> >,
  RecordField<VecField::py, NumpyBuilder<double> >,
  RecordField<VecField::pz, NumpyBuilder<double> >,
  RecordField<VecField::e, NumpyBuilder<double> > >;
using PrtBuilder = RecordBuilder<
  RecordField<PrtField::id, NumpyBuilder<int32_t> >,
  RecordField<PrtField::status, NumpyBuilder<int32_t> >,
  RecordField<PrtField::mother1, NumpyBuilder<int32_t> >,
  RecordField<PrtField::mother2, NumpyBuilder<int32_t> >,
  RecordField<PrtField::daughter1, NumpyBuilder<int32_t> >,
  RecordField<PrtField::daughter2, NumpyBuilder<int32_t> >,
  RecordField<PrtField::col, NumpyBuilder<int32_t> >,
  RecordField<PrtField::acol, NumpyBuilder<int32_t> >,
  RecordField<PrtField::p, VecBuilder>,
  RecordField<PrtField::m, NumpyBuilder<double> >,
  RecordField<PrtField::scale, NumpyBuilder<double> >,
  RecordField<PrtField::pol, NumpyBuilder<double> >,
  RecordField<PrtField::vProd, IndexedOptionBuilder<int32_t, VecBuilder> >,
  RecordField<PrtField::tau, NumpyBuilder<double> > >;
using InfoBuilder = RecordBuilder<
  RecordField<InfoField::id1, NumpyBuilder<int32_t> >,
  RecordField<InfoField::id2, NumpyBuilder<int32_t> >,
  RecordField<InfoField::x1, NumpyBuilder<double> >,
  RecordField<InfoField::x2, NumpyBuilder<double> >,
  RecordField<InfoField::pdf1, NumpyBuilder<double> >,
  RecordField<InfoField::pdf2, NumpyBuilder<double> >,
  RecordField<InfoField::alphaS, NumpyBuilder<double> >,
  RecordField<InfoField::alphaEM, NumpyBuilder<double> >,
  RecordField<InfoField::Q2Fac, NumpyBuilder<double> >,
  RecordField<InfoField::Q2Ren, NumpyBuilder<double> >,
  RecordField<InfoField::mHat, NumpyBuilder<double> >,
  RecordField<InfoField::sHat, NumpyBuilder<double> >,
  RecordField<InfoField::tHat, NumpyBuilder<double> >,
  RecordField<InfoField::uHat, NumpyBuilder<double> >,
  RecordField<InfoField::pT2Hat, NumpyBuilder<double> >,
  RecordField<InfoField::weights,
              ListOffsetBuilder<int64_t, NumpyBuilder<double> > > >;
using EventBuilder = RecordBuilder<
  RecordField<EventField::prt, ListOffsetBuilder<int64_t, PrtBuilder> >,
  RecordField<EventField::info, InfoBuilder> >;

// Define the keys and maps for the builders.
const std::string vecKey("\"__record__\": \"Momentum4D\"");
const UserDefinedMap vecMap({
    {VecField::px, "px"},
    {VecField::py, "py"},
    {VecField::pz, "pz"},
    {VecField::e, "e"}
  });
const std::string prtKey("\"__record__\": \"PythiaParticle\"");
const UserDefinedMap prtMap({
    {PrtField::id, "id"},
    {PrtField::status, "status"},
    {PrtField::mother1, "mother1"},
    {PrtField::mother2, "mother2"},
    {PrtField::daughter1, "daughter1"},
    {PrtField::daughter2, "daughter2"},
    {PrtField::col, "col"},
    {PrtField::acol, "acol"},
    {PrtField::m, "m"},
    {PrtField::scale, "scale"},
    {PrtField::pol, "pol"},
    {PrtField::tau, "tau"},
    {PrtField::vProd, "vProd"},
    {PrtField::p, "p"},
  });
const std::string infoKey("\"__record__\": \"PythiaInfo\"");
const UserDefinedMap infoMap({
    {InfoField::id1, "id1"},
    {InfoField::id2, "id2"},
    {InfoField::x1, "x1"},
    {InfoField::x2, "x2"},
    {InfoField::pdf1, "pdf1"},
    {InfoField::pdf2, "pdf2"},
    {InfoField::alphaS, "alphaS"},
    {InfoField::alphaEM, "alphaEM"},
    {InfoField::Q2Fac, "Q2Fac"},
    {InfoField::Q2Ren, "Q2Ren"},
    {InfoField::mHat, "mHat"},
    {InfoField::sHat, "sHat"},
    {InfoField::tHat, "tHat"},
    {InfoField::uHat, "uHat"},
    {InfoField::pT2Hat, "pT2Hat"},
    {InfoField::weights, "weights"}
  });
const std::string eventKey("\"__record__\": \"PythiaEvent\"");
const UserDefinedMap eventMap({
    {EventField::prt, "prt"},
    {EventField::info, "info"}
  });

//==========================================================================

// Class which constructs an Awkward array from a run.

class Run {

 public:

  // Constructor.
  Run(int nEventsIn, pybind11::object errorMode);

  // Destructor.
  ~Run();
  
  // Generate the next run of events for a single Pythia instance.
  pybind11::object next(Pythia8::Pythia* pythiaPtr);
  
  // Fill an event, info, or particles.
  void fillEvent(Pythia8::Pythia* pythiaPtr);
  void fillInfo(InfoBuilder& builder, const Pythia8::Info& info);
  void fillPrt(PrtBuilder &builder, const Pythia8::Event& event);
  
  // Fill an array.
  template<typename T> pybind11::object fillArray(const T &builder) {
    
    // Determine necessary memory.
    std::map<std::string, std::size_t> bytes = {};
    builder.buffer_nbytes(bytes);
    
    // Allocate the memory.
    std::map<std::string, void*> buffers = {};
    for (auto it: bytes) buffers[it.first] = (void*)(new uint8_t[it.second]);
    
    // Write non-contiguous contents to memory.
    builder.to_buffers(buffers);
    
    // Build Python dictionary containing arrays. The dtypes are not
    // important here as long as they match the underlying buffer as
    // Awkward Array calls `frombuffer` to convert to the correct type
    pybind11::dict container;
    for (auto it: buffers)
      container[pybind11::str(it.first)] = pybind11::array_t<uint8_t>(
        {bytes[it.first]}, {sizeof(uint8_t)},
        reinterpret_cast<uint8_t*>(it.second),
        pybind11::capsule(it.second, [](void *data) {
            uint8_t *dataPtr = reinterpret_cast<uint8_t *>(data);
            delete[] dataPtr;})
        );
    return pybind11::module::import("awkward").attr("from_buffers")
      (builder.form(), builder.length(), container);
    
  }
  
  // Builders.
  IndexedOptionBuilder<int64_t, EventBuilder>* noneBuilder{};
  EventBuilder* eventBuilder{};
  PrtBuilder* prtBuilder{};
  InfoBuilder* infoBuilder{};
  
  // Error mode flags.
  bool valid{false},
    isNone{false}, isSkip{false}, isFail{false}, isFloat{false};

  // Event counters.
  int nEvents{0}, nTries{0}, nAcc{0}, nTry{0};
  
};

//==========================================================================

// Method to generate next batch of events from a single Pythia
// instance. errorMode defines how failed events are handled and can
// take on the following values.
//     "none": failed events are included in the record as None.
//     "skip": failed events are not included in the record.
//     "fail": this method fails if any single event fails.
//      <fac>: if events fail, additional events are generated until nEvents
//             is reached. A maximum of <fac>*nEvents events are tried,
//             where <fac> is >= 1. If -1, then infinite events are tried.

pybind11::object nextBatch(Pythia8::Pythia* pythiaPtr, int nEvents,
  pybind11::object errorMode = pybind11::str("skip"));

//==========================================================================

// Method to generate next batch of events from a PythiaParallel
// instance. This instance only runs in the "skip" error mode.

pybind11::object nextBatchParallel(
  Pythia8::PythiaParallel* pythiaPtr, int nEvents);

//==========================================================================

#endif // awkward_PythiaEvent_H
