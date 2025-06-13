// PythiaBatch.cpp is a part of the PYTHIA event generator.
// Copyright (C) 2025 Aartem Havryliuk and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Definitions for batch events in Python via the Awkward package.
// This file is NOT automatically generated.

#include "awkward/PythiaBatch.h"

//==========================================================================

// Class which constructs an Awkward array from a run.

//--------------------------------------------------------------------------

// Constructor.

Run::Run(int nEventsIn, pybind11::object errorMode) {
  
  // Check and define the error mode.
  if (!pybind11::isinstance<pybind11::str>(errorMode) &&
    !pybind11::isinstance<pybind11::float_>(errorMode)) {
    std::cerr << "errorMode must be str or float.\n";
    valid = false;
    return;
  }
  isNone = pybind11::isinstance<pybind11::str>(errorMode) &&
    errorMode.cast<std::string>() == "none";
  isSkip = pybind11::isinstance<pybind11::str>(errorMode) &&
    errorMode.cast<std::string>() == "skip";
  isFail = pybind11::isinstance<pybind11::str>(errorMode) &&
    errorMode.cast<std::string>() == "fail";
  isFloat = pybind11::isinstance<pybind11::float_>(errorMode);
  
  // Determine the maximum number of trials.
  nEvents = nEventsIn;
  nTries  = nEventsIn;
  if (isFloat) {
    float fac = errorMode.cast<float>();
    if (fac == -1) nTries = std::numeric_limits<int>::max();
    else if (fac < 1) throw std::invalid_argument(
      "errorMode as float must be >= 1, or exactly -1 for infinite trials"); 
    else nTries = static_cast<int>(fac * nEvents);
  } else if (!isSkip && !isFail && !isNone) throw std::invalid_argument(
    "errorMode must be 'skip', 'fail', 'none', or a float.");
  
  // Set the event counters.
  nAcc = 0;
  nTry = 0;
  
  // Set the event builders.
  if (isNone) {
    noneBuilder = new IndexedOptionBuilder<int64_t, EventBuilder>();
  } else {
    eventBuilder = new EventBuilder();
    eventBuilder->set_fields(eventMap);
    eventBuilder->set_parameters(eventKey);
    prtBuilder = &eventBuilder->content<EventField::prt>().begin_list();
    infoBuilder = &eventBuilder->content<EventField::info>();
    infoBuilder->set_fields(infoMap);
    infoBuilder->set_parameters(infoKey);
  }
  valid = true;
  
}

//--------------------------------------------------------------------------

// Destructor.

Run::~Run() {if (isNone) delete noneBuilder; else delete eventBuilder;}

//--------------------------------------------------------------------------

// Generate the next run of events for a single Pythia instance.

pybind11::object Run::next(Pythia8::Pythia* pythiaPtr) {
  
  // The event loop.
  if (!valid) return pybind11::none();
  while (nAcc < nEvents && nTry < nTries) {
    ++nTry;
    if (!pythiaPtr->next()) {
      if (isNone) noneBuilder->append_invalid();
      else if (isFail) throw std::runtime_error(
        "failed to generate event " + std::to_string(nTry));
    } else fillEvent(pythiaPtr);
  }
  
  // Exit if maximum tries reached.
  if (isFloat && nAcc < nEvents) throw std::runtime_error(
    "exceeded allowed trials " + std::to_string(nTries));
  
  // Fill the Awkward array.
  if (isNone) return fillArray(*noneBuilder);
  else return fillArray(*eventBuilder);
  
}

//--------------------------------------------------------------------------

// Fill an event.

void Run::fillEvent(Pythia8::Pythia* pythiaPtr) {

  // Get the builders if "none" error mode.
  if (isNone) {
    eventBuilder = &noneBuilder->append_valid();
    eventBuilder->set_fields(eventMap);
    eventBuilder->set_parameters(eventKey);
    prtBuilder  = &eventBuilder->content<EventField::prt>().begin_list();
    infoBuilder = &eventBuilder->content<EventField::info>();
    infoBuilder->set_fields(infoMap);
    infoBuilder->set_parameters(infoKey);
  }
  
  // Set the particle builder.
  prtBuilder->set_fields(prtMap);
  prtBuilder->set_parameters(prtKey);
  
  // Fill the info.
  fillInfo(*infoBuilder, pythiaPtr->info);
  fillPrt(*prtBuilder, pythiaPtr->event);
  eventBuilder->content<EventField::prt>().end_list();
  ++nAcc;
  
}

//--------------------------------------------------------------------------

// Fill info.

void Run::fillInfo(InfoBuilder& builder, const Pythia8::Info& info) {

  // Get the info attribute builders.
  auto &id1     = builder.content<InfoField::id1>();
  auto &id2     = builder.content<InfoField::id2>();
  auto &x1      = builder.content<InfoField::x1>();
  auto &x2      = builder.content<InfoField::x2>();
  auto &pdf1    = builder.content<InfoField::pdf1>();
  auto &pdf2    = builder.content<InfoField::pdf2>();
  auto &alphaS  = builder.content<InfoField::alphaS>();
  auto &alphaEM = builder.content<InfoField::alphaEM>();
  auto &Q2Fac   = builder.content<InfoField::Q2Fac>();
  auto &Q2Ren   = builder.content<InfoField::Q2Ren>();
  auto &mHat    = builder.content<InfoField::mHat>();
  auto &sHat    = builder.content<InfoField::sHat>();
  auto &tHat    = builder.content<InfoField::tHat>();
  auto &uHat    = builder.content<InfoField::uHat>();
  auto &pT2Hat  = builder.content<InfoField::pT2Hat>();
  auto &weights = builder.content<InfoField::weights>().begin_list();
  
  // Set the info.
  id1.append(info.id1());
  id2.append(info.id2());
  x1.append(info.x1());
  x2.append(info.x2());
  pdf1.append(info.pdf1());
  pdf2.append(info.pdf2());
  alphaS.append(info.alphaS());
  alphaEM.append(info.alphaEM());
  Q2Fac.append(info.Q2Fac());
  Q2Ren.append(info.Q2Ren());
  mHat.append(info.mHat());
  sHat.append(info.sHat());
  tHat.append(info.tHat());
  uHat.append(info.uHat());
  pT2Hat.append(info.pT2Hat());
  for (int iWgt = 0; iWgt < info.nWeightGroups(); ++iWgt)
    weights.append(info.getGroupWeight(iWgt));
  builder.content<InfoField::weights>().end_list();
  
}

//--------------------------------------------------------------------------

// Fill particles.

void Run::fillPrt(PrtBuilder &builder, const Pythia8::Event& event) {
  
  // Get the attribute builders.
  auto &id        = builder.content<PrtField::id>();
  auto &status    = builder.content<PrtField::status>();
  auto &mother1   = builder.content<PrtField::mother1>();
  auto &mother2   = builder.content<PrtField::mother2>();
  auto &daughter1 = builder.content<PrtField::daughter1>();
  auto &daughter2 = builder.content<PrtField::daughter2>();
  auto &col       = builder.content<PrtField::col>();
  auto &acol      = builder.content<PrtField::acol>();
  auto &m         = builder.content<PrtField::m>();
  auto &scale     = builder.content<PrtField::scale>();
  auto &pol       = builder.content<PrtField::pol>();
  auto &tau       = builder.content<PrtField::tau>();
  auto &p         = builder.content<PrtField::p>();
  
  // Get the momentum attribute builders.
  p.set_fields(vecMap);
  p.set_parameters(vecKey);
  auto &px = p.content<VecField::px>();
  auto &py = p.content<VecField::py>();
  auto &pz = p.content<VecField::pz>();
  auto &e  = p.content<VecField::e>();
  
  // Get the optional production vertex builders.
  auto &vProd = builder.content<PrtField::vProd>();
  auto &vProdSub = vProd.append_valid();
  vProdSub.set_fields(vecMap);
  vProdSub.set_parameters(vecKey);
  auto &xProd = vProdSub.content<VecField::px>();
  auto &yProd = vProdSub.content<VecField::py>();
  auto &zProd = vProdSub.content<VecField::pz>();
  auto &tProd = vProdSub.content<VecField::e>();
  
  // Loop over the particles.
  for (int iPrt = 0; iPrt < event.size(); ++iPrt) {
    const Pythia8::Particle& prt = event[iPrt];
    status.append(prt.status());
    id.append(prt.id());
    mother1.append(prt.mother1());
    mother2.append(prt.mother2());
    daughter1.append(prt.daughter1());
    daughter2.append(prt.daughter2());
    col.append(prt.col());
    acol.append(prt.acol());
    m.append(prt.m());
    scale.append(prt.scale());
    pol.append(prt.pol());
    tau.append(prt.tau());

    // Set the momentum.
    px.append(prt.px());
    py.append(prt.py());
    pz.append(prt.pz());
    e.append(prt.e());

    // Set the vertex.
    if (prt.hasVertex()) {
      vProd.append_valid();
      xProd.append(prt.xProd());
      yProd.append(prt.yProd());
      zProd.append(prt.zProd());
      tProd.append(prt.tProd());
    } else vProd.append_invalid();
  }
  
}

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
  pybind11::object errorMode) {
  Run run(nEvents, errorMode);
  return run.next(pythiaPtr);
}

//==========================================================================

// Method to generate next batch of events from a PythiaParallel
// instance. This instance only runs in the "skip" error mode.

pybind11::object nextBatchParallel(
  Pythia8::PythiaParallel* pythiaPtr, int nEvents) {
  Run run(nEvents, pybind11::str("skip"));
  pythiaPtr->run(nEvents, [&run](Pythia8::Pythia* pythiaPtr) {
      run.fillEvent(pythiaPtr);});
  return run.fillArray(*run.eventBuilder);
}

//==========================================================================
