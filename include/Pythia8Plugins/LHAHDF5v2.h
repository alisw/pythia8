// LHAHDF5v2.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Authors: Christian Preuss, Stefan Hoeche December 2023.

#ifndef Pythia8_LHAHDF5v2_H
#define Pythia8_LHAHDF5v2_H

// Interface includes.
#include "Pythia8Plugins/LHEH5v2.h"

// Generator includes.
#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// HDF5 version 2 file reader.
// Converts to Pythia-internal events by acting as replacement
// Les Houches Event reader.

class LHAupH5v2 : public Pythia8::LHAup {

 public:

  LHAupH5v2(HighFive::File* h5fileIn, size_t firstEventIn, size_t readSizeIn,
    bool normalize) :
    lhefPtr(new LHEH5::LHEFile()), readSizeSav(readSizeIn),
      nReadSav(0), nTrialsSav(0) {

    // Read event-file header and events.
    lhefPtr->ReadHeader(*h5fileIn);
    lhefPtr->ReadEvents(*h5fileIn, firstEventIn, readSizeIn);

    // This reads the init information.
    _weightnames = lhefPtr->WeightNames();
    if (normalize) lhefPtr->Scale(readSizeIn/lhefPtr->SumTrials());
    std::vector<double> info;
    h5fileIn->getDataSet("init").read(info);
    setBeamA(info[0],info[2],info[4],info[6]);
    setBeamB(info[1],info[3],info[5],info[7]);
    setStrategy(-4);
    int numProcesses = info[9];
    vector<int> procId(numProcesses);
    vector<double> xSection(numProcesses);
    vector<double> error(numProcesses);
    vector<double> unitWeight(numProcesses);
    for (int i = 0; i<numProcesses; ++i) {
      LHEH5::ProcInfo pi(lhefPtr->GetProcInfo(i));
      procId[i]     = pi.pid;
      xSection[i]   = pi.xsec;
      error[i]      = pi.error;
      unitWeight[i] = pi.unitwgt;
    }
    for (int np = 0; np<numProcesses; ++np) {
      addProcess(procId[np], xSection[np], error[np], unitWeight[np]);
      xSecSumSave += xSection[np];
      xErrSumSave += pow2(error[np]);
    }
  }

  ~LHAupH5v2() {delete lhefPtr;}

  bool setInit() override {return true;}
  bool setEvent(int idProc=0) override;
  void forceStrategy(int strategyIn) {setStrategy(strategyIn);}
  size_t nTrials() {return nTrialsSav;}
  size_t nRead()   {return nReadSav;}

private:

  // HDF5 event file.
  LHEH5::LHEFile *lhefPtr;

  // Info for reader.
  size_t readSizeSav, nReadSav, nTrialsSav;

  // Multiweight vector.
  vector<double> _eventweightvalues;
  vector<string> _weightnames;

  // Particle production scales.
  LHAscales scalesNow;

};

//--------------------------------------------------------------------------

// Read an event.

bool LHAupH5v2::setEvent(int) {

  // Equivalent of end of file.
  if (nReadSav >= readSizeSav) return false;

  // Read event.
  LHEH5::Event evt(lhefPtr->GetEvent(nReadSav));
  if (evt[0].pz<0 && evt[1].pz>0) swap<LHEH5::Particle>(evt[0], evt[1]);

  setProcess(evt.pinfo.pid, evt.wgts[0], evt.mur, evt.aqed, evt.aqcd);
  nupSave    = evt.size();
  idprupSave = evt.pinfo.pid;
  xwgtupSave = evt.wgts[0];
  scalupSave = evt.mur;
  aqedupSave = evt.aqed;
  aqcdupSave = evt.aqcd;
  double scalein = -1.;

  // Communicate event weight to Info.
  _eventweightvalues=evt.wgts;
  // infoPtr->weights_compressed_names = &_weightnames;
  infoPtr->weights_compressed = &_eventweightvalues;

  // Set particles.
  for (unsigned int ip=0; ip<evt.size(); ++ip) {
    const LHEH5::Particle& p = evt[ip];
    if (ip < 2) addParticle(p.id, p.st, 0, 0,
      p.cl1, p.cl2, p.px, p.py, p.pz, p.e, p.m,
      p.lt, p.sp, scalein);
    else addParticle(p.id, p.st, p.mo1, p.mo2,
      p.cl1, p.cl2, p.px, p.py, p.pz, p.e, p.m,
      p.lt, p.sp, scalein);
  }

  // Scale setting
  scalesNow.clear();
  scalesNow.muf   = evt.muf;
  scalesNow.mur   = evt.mur;
  scalesNow.mups  = evt.muq;
  infoPtr->scales = &scalesNow;
  infoPtr->setEventAttribute("npLO",std::to_string(evt.pinfo.nplo));
  infoPtr->setEventAttribute("npNLO",std::to_string(evt.pinfo.npnlo));

  // Update counters.
  ++nReadSav;
  nTrialsSav += evt.trials;
  return true;
}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LHAHDF5v2_H
