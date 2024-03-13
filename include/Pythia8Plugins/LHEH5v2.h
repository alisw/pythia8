// The LHEH5 code below has been adapted from
// https://gitlab.com/hpcgen/tools authored by Holger Schulz
// and Stefan Hoeche, and developed under the Scientific Discovery
// through Advanced Computing (SciDAC) program funded by
// the U.S. Department of Energy, Office of Science, Advanced Scientific
// Computing Research.
//
// Note, this header can be used in conjuction with LHAHF5v2.h.
//
// Fermilab Software Legal Information (BSD License)
// Copyright (c) 2009, FERMI NATIONAL ACCELERATOR LABORATORY
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the
// distribution.
//
// Neither the name of the FERMI NATIONAL ACCELERATOR LABORATORY, nor
// the names of its contributors may be used to endorse or promote
// products derived from this software without specific prior written
// permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef Pythia8_LHEH5v2_H
#define Pythia8_LHEH5v2_H

// Standard includes.
#include <iostream>
#include <string>
#include <vector>

// HighFive includes.
#include "H5File.hpp"
#include "H5FileDriver.hpp"
#include "H5DataSet.hpp"

using namespace HighFive;

namespace LHEH5 {

//==========================================================================

// Process info struct.

//--------------------------------------------------------------------------

struct ProcInfo {
  int pid, nplo, npnlo;
  double unitwgt, xsec, error;

  inline ProcInfo(int _pid,int _nplo,int _npnlo,
    double _unitwgt,double _xsec,double _error):
    pid(_pid), nplo(_nplo), npnlo(_npnlo),
      unitwgt(_unitwgt), xsec(_xsec), error(_error) {}
};

//--------------------------------------------------------------------------

// Print process info.

std::ostream& operator<<(std::ostream& s,const ProcInfo& p) {
  return s<<"[pid="<<p.pid<<",nplo="<<p.nplo<<",npnlo="<<p.npnlo<<",xsec="
          <<p.xsec<<",err="<<p.error<<",unitwgt="<<p.unitwgt<<"]";
}

//==========================================================================

// Particle struct.

//--------------------------------------------------------------------------

struct Particle {
  int id, st, mo1, mo2, cl1, cl2;
  double px, py, pz, e, m, lt, sp;

  inline Particle(int _id,int _st,int _mo1,int _mo2,int _cl1,int _cl2,
    double _px,double _py,double _pz,double _e,double _m,
    double _lt,double _sp):
      id(_id), st(_st), mo1(_mo1), mo2(_mo2), cl1(_cl1), cl2(_cl2),
      px(_px), py(_py), pz(_pz), e(_e), m(_m), lt(_lt), sp(_sp) {}
  };// end of struct Particle

//--------------------------------------------------------------------------

// Print a particle

std::ostream& operator<<(std::ostream& s,const Particle& p) {
  return s << "{id=" << p.id << ",st=" << p.st
           << ",mo=[" << p.mo1 << "," << p.mo2 << "]"
           << ",cl=[" << p.cl1 << "," << p.cl2 << "]"
           << ",p=(" << p.e << "," << p.px << ","
           << p.py << "," << p.pz << ")}";
}

//==========================================================================

// Event struct.

//--------------------------------------------------------------------------

struct Event : public std::vector<Particle> {
  ProcInfo pinfo;
  size_t trials;
  std::vector<double> wgts;
  double mur, muf, muq, aqed, aqcd;
  std::vector<Particle> ctparts;
  int ijt, kt, i, j, k;
  double z1, z2, bbpsw, psw;

  inline Event(const ProcInfo &_pinfo,
    size_t _trials, std::vector<double> _wgts,
    double _mur, double _muf, double _muq,
    double _aqed, double _aqcd):
    pinfo(_pinfo), trials(_trials), wgts(_wgts),
      mur(_mur), muf(_muf), muq(_muq), aqed(_aqed), aqcd(_aqcd),
      ijt(-1), kt(-1), i(-1), j(-1), k(-1),
      z1(0), z2(0), bbpsw(0), psw(0) {}
  inline void AddCTInfo(int _ijt, int _kt, int _i, int _j, int _k,
    double _z1, double _z2, double _bbw, double _w) {
    ijt=_ijt; kt=_kt, i=_i; j=_j; k=_k;
    z1=_z1; z2=_z2; bbpsw=_bbw; psw=_w;
  }

};

//--------------------------------------------------------------------------

// Print an event.

std::ostream& operator<<(std::ostream& s, const Event& e) {
  s << "Event " << e.pinfo << " {\n"
    << "  trials=" << e.trials << ",weights=(" << e.wgts[0];
  for (size_t i(1); i<e.wgts.size(); ++i) s << "," << e.wgts[i];
  s << ")\n  mur=" << e.mur << ", muf=" << e.muf << ", muq=" << e.muq
    << ",aqed=" << e.aqed << ",aqcd=" << e.aqcd << "\n";
  for (size_t i(0); i<e.size(); ++i) s << "  " << e[i] << "\n";
  if (!e.ctparts.empty() || e.psw) {
    s << "  (" << e.ijt << "," << e.kt << ")->(" << e.i << "," << e.j
      << "," << e.k << "), z1=" << e.z1 << ", z2=" << e.z2
      << ", bbpsw=" << e.bbpsw << ", psw=" << e.psw << "\n";
    for (size_t i(0); i<e.ctparts.size(); ++i) s<<"  "<<e.ctparts[i]<<"\n";
  }
  return s<<"}";
}

//==========================================================================

// LHEFile class.

//--------------------------------------------------------------------------

class LHEFile {

 public:

  inline const std::vector<int>& Version() const { return version; }

  inline double TotalXS() const {
    double xs(0.);
    for (size_t i(0); i<pinfo.size(); ++i) xs += pinfo[i][3];
    return xs;
  }

  inline double SumTrials() const {
    double trials(0.);
    for (size_t i(0); i<evts.size(); ++i) trials += evts[i][3];
    return trials;
  }

  inline double UnitWeight() const {
    double wgt(0.);
    for (size_t i(0); i<pinfo.size(); ++i) wgt += pinfo[i][5];
    // For version 2.0.0 multiply by pb/GeV^2.
    if (version[0]==2 && version[1]==0 && version[2]==0) wgt*=3.89379656e8;
    return wgt;
  }

  inline const std::vector<std::string>& WeightNames() const {
    return wgtnames;
  }

  inline size_t NProcesses() const { return pinfo.size(); }
  inline ProcInfo GetProcInfo(const size_t pid) const {
    return ProcInfo(pid, pinfo[pid][1], pinfo[pid][2],
      pinfo[pid][5], pinfo[pid][3], pinfo[pid][4]);
  }

  inline size_t NEvents() const { return evts.size(); }
  inline Event GetEvent(size_t i) const {
    std::vector<double> wgts(evts[i].begin()+9, evts[i].end());
    Event e(GetProcInfo(evts[i][0] ? evts[i][0]-1 : 0), evts[i][3], wgts,
      evts[i][6], evts[i][5], evts[i][4], evts[i][7], evts[i][8]);
    double wgt(0.);
    for (std::vector<double>::const_iterator it(wgts.begin());
         it!=wgts.end(); ++it) wgt += std::abs(*it);
    if (!wgt) return e;
    for (int n(0); n<evts[i][1]; ++n)
      e.push_back(GetParticle(evts[i][2]-evts[0][2]+n));
    if (!ctevts.empty()) {
      e.AddCTInfo(ctevts[i][0], ctevts[i][1], ctevts[i][2],
        ctevts[i][3], ctevts[i][4], ctevts[i][5],
        ctevts[i][6], ctevts[i][7], ctevts[i][8]);
      if (ctevts[i][0]>=0 && ctevts[i][1]>=0)
        for (int n(0); n<evts[i][1]+(ctevts[i][0]>=0 ? 1 : 0); ++n)
          e.ctparts.push_back(GetCTParticle(evts[i][2]-evts[0][2]+n));
    }
    return e;
  }

  void Scale(const double s) {
    for (size_t i(0); i<evts.size(); ++i)
      for (size_t j(9); j<evts[i].size(); ++j) evts[i][j] *= s;
  }

  void ReadHeader(File &file) {
    auto xfer_props = DataTransferProps{};
#ifdef USING__MPI
    xfer_props.add(UseCollectiveIO{});
#endif
    file.getDataSet("version").read(version,xfer_props);
    file.getDataSet("procInfo").read(pinfo,xfer_props);
    DataSet events(file.getDataSet("events"));
    auto attr_keys(events.listAttributeNames());
    Attribute a(events.getAttribute(attr_keys[0]));
    a.read(wgtnames);
    for (int i(0); i<9; ++i) wgtnames.erase(wgtnames.begin());
  }

  void ReadEvents(File &file, size_t first_event, size_t n_events) {
    auto xfer_props = DataTransferProps{};
#ifdef USING__MPI
    xfer_props.add(UseCollectiveIO{});
#endif
    DataSet events(file.getDataSet("events"));
    std::vector<size_t> eoffsets{first_event, 0};
    std::vector<size_t> ecounts{n_events, 9+wgtnames.size()};
    evts.resize(n_events,std::vector<double>(9+wgtnames.size()));
    events.select(eoffsets,ecounts).read(evts, xfer_props);
    DataSet particles(file.getDataSet("particles"));
    std::vector<size_t> poffsets{(size_t)evts.front()[2], 0};
    size_t nmax(0);
    for (size_t i(0); i<pinfo.size(); ++i)
      nmax = std::max((size_t)std::max(pinfo[i][1],pinfo[i][2]+1), nmax);
    std::vector<size_t> pcounts{n_events*nmax, 13};
    parts.resize(n_events*nmax, std::vector<double>(13));
    particles.select(poffsets, pcounts).read(parts, xfer_props);
    if (file.exist("ctevents")) {
      DataSet ctevents(file.getDataSet("ctevents"));
      std::vector<size_t> cteoffsets{first_event, 0};
      std::vector<size_t> ctecounts{n_events, 9};
      ctevts.resize(n_events, std::vector<double>(9));
      ctevents.select(cteoffsets, ctecounts).read(ctevts, xfer_props);
      DataSet ctparticles(file.getDataSet("ctparticles"));
      std::vector<size_t> ctpoffsets{(size_t)evts.front()[2],0};
      std::vector<size_t> ctpcounts{n_events*nmax, 4};
      ctparts.resize(n_events*nmax, std::vector<double>(4));
      ctparticles.select(ctpoffsets, ctpcounts).read(ctparts, xfer_props);
    }
  }

 private:

  std::vector<int> version;
  std::vector<std::vector<double> > evts, parts, pinfo;
  std::vector<std::vector<double> > ctevts, ctparts;
  std::vector<std::string> wgtnames;

  inline Particle GetParticle(size_t i) const {
    return Particle(parts[i][0], parts[i][1], parts[i][2], parts[i][3],
      parts[i][4], parts[i][5], parts[i][6], parts[i][7],
      parts[i][8], parts[i][9], parts[i][10],
      parts[i][11], parts[i][12]);
  }

  inline Particle GetCTParticle(size_t i) const {
    return Particle(-1, -1, -1, -1, -1, -1, ctparts[i][0], ctparts[i][1],
      ctparts[i][2], ctparts[i][3], -1, -1, -1);
  }

};

//==========================================================================

}// end of namespace LHEH5

#endif // Pythia8_LHEH5v2_H
