// LHAPDF6.h is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the LHAPDF6 PDF plugin class.

#ifndef Pythia8_LHAPDF6_H
#define Pythia8_LHAPDF6_H

#include "Pythia8/PartonDistributions.h"
#include "Pythia8/Plugins.h"
#include "LHAPDF/LHAPDF.h"
#include <mutex>

namespace Pythia8 {

//==========================================================================

// Containers for PDF sets.

//--------------------------------------------------------------------------

// Class to hold a PDF set, its information, and its uncertainty sets.

class PdfSets {

public:

  // Constructors.
  PdfSets() {;}
  PdfSets(string setName) : info(::LHAPDF::PDFSet(setName)),
    pdfs(vector< ::LHAPDF::PDF* >(info.size(), 0)) {;}

  // Access a PDF set.
  ::LHAPDF::PDF *operator[](unsigned int member) {
    if (!pdfs[member]) {
      lock_guard<mutex> lck (mtx);
      pdfs[member] = info.mkPDF(member);
    }
    return pdfs[member];
  }

  // Get number of PDF sets.
  int size() {return pdfs.size();}

  // PDF sets and info.
  ::LHAPDF::PDFSet info;
  vector< ::LHAPDF::PDF* > pdfs;
  static mutex mtx;

};

mutex PdfSets::mtx;

//==========================================================================

// Provide interface to the LHAPDF6 library of parton densities.

class LHAPDF6 : public PDF {

public:

  // Constructor.
  LHAPDF6(Pythia*, Settings* settingsPtr, Logger*) :
    PDF(), pdf(nullptr), extrapol(false) {
    if (settingsPtr == nullptr) return;
    sSymmetric(settingsPtr->flag("LHAPDF:sSymmetric"));
    cSymmetric(settingsPtr->flag("LHAPDF:cSymmetric"));
    bSymmetric(settingsPtr->flag("LHAPDF:bSymmetric"));
  }

  // Initialization of PDF set.
  bool init(int idBeamIn, string setName, int member, Logger* loggerPtr)
    override;

  // Allow extrapolation beyond boundaries (not implemented).
  void setExtrapolate(bool extrapolIn) {extrapol = extrapolIn;}

private:

  // The LHAPDF objects.
  PdfSets pdfs;
  ::LHAPDF::PDF *pdf;
  ::LHAPDF::Extrapolator *ext;
  bool extrapol;

  // Update parton densities.
  void xfUpdate(int id, double x, double Q2);

  // Check whether x and Q2 values fall inside the fit bounds.
  bool insideBounds(double x, double Q2) {
    return (x > xMin && x < xMax && Q2 > q2Min && Q2 < q2Max);}

  // Return the running alpha_s shipped with the LHAPDF set.
  double alphaS(double Q2) { return pdf->alphasQ2(Q2); }

  // Return quark masses used in the PDF fit.
  double muPDFSave, mdPDFSave, mcPDFSave, msPDFSave, mbPDFSave,
         xMin, xMax, q2Min, q2Max;
  double mQuarkPDF(int id) {
    switch(abs(id)){
      case 1: return mdPDFSave;
      case 2: return muPDFSave;
      case 3: return msPDFSave;
      case 4: return mcPDFSave;
      case 5: return mbPDFSave;
    }
    return -1.;
 }

  // Calculate uncertainties using the LHAPDF prescription.
  void calcPDFEnvelope(int, double, double, int);
  void calcPDFEnvelope(pair<int,int>, pair<double,double>, double, int);
  PDFEnvelope pdfEnvelope;
  PDFEnvelope getPDFEnvelope() {return pdfEnvelope;}
  static const double PDFMINVALUE;

  int nMembersSave;
  int nMembers() { return nMembersSave; }

};

//--------------------------------------------------------------------------

// Constants.

const double LHAPDF6::PDFMINVALUE = 1e-10;

//--------------------------------------------------------------------------

// Initialize a parton density function from LHAPDF6.

bool LHAPDF6::init(int idBeamIn, string setName, int member,
  Logger* loggerPtr) {
  idBeam = idBeamIn;
  idBeamAbs = abs(idBeamIn);
  isSet = false;

  // Find the PDF set. Note, LHAPDF aborts if the PDF does not exist,
  // which we avoid with this try/catch statement. Ideally, we would
  // check with :LHAPDF::lookupLHAPDFID, but this is not thread safe.
  try {
    pdfs = PdfSets(setName);
  } catch (const std::exception &e) {
    loggerPtr->ERROR_MSG("unknown PDF " + setName);
    return false;
  }

  // Find the PDF member.
  if (pdfs.size() == 0) {
    loggerPtr->ERROR_MSG("could not initialize PDF " + setName);
    return false;
  } else if (member >= pdfs.size()) {
    loggerPtr->ERROR_MSG(setName + " does not contain requested member");
    return false;
  }
  pdf = pdfs[member];
  isSet = true;

  // Save x and Q2 limits.
  xMax  = pdf->xMax();
  xMin  = pdf->xMin();
  q2Max = pdf->q2Max();
  q2Min = pdf->q2Min();

  // Store quark masses used in PDF fit.
  muPDFSave = pdf->info().get_entry_as<double>("MUp");
  mdPDFSave = pdf->info().get_entry_as<double>("MDown");
  mcPDFSave = pdf->info().get_entry_as<double>("MCharm");
  msPDFSave = pdf->info().get_entry_as<double>("MStrange");
  mbPDFSave = pdf->info().get_entry_as<double>("MBottom");
  nMembersSave  = pdf->info().get_entry_as<int>("NumMembers");
  return true;

}

//--------------------------------------------------------------------------

// Give the parton distribution function set from LHAPDF6.

void LHAPDF6::xfUpdate(int, double x, double Q2) {
  if (!isSet) return;

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  if (x < xMin && !extrapol) x = xMin;
  if (x > xMax)    x = xMax;
  if (Q2 < q2Min) Q2 = q2Min;
  if (Q2 > q2Max) Q2 = q2Max;

  // Update values.
  xg     = pdf->xfxQ2(21, x, Q2);
  xd     = pdf->xfxQ2(1,  x, Q2);
  xu     = pdf->xfxQ2(2,  x, Q2);
  xdbar  = pdf->xfxQ2(-1, x, Q2);
  xubar  = pdf->xfxQ2(-2, x, Q2);
  xs     = pdf->xfxQ2(3,  x, Q2);
  xc     = pdf->xfxQ2(4,  x, Q2);
  xb     = pdf->xfxQ2(5,  x, Q2);
  xsbar  = sSymmetricSave ? xs : pdf->xfxQ2(-3, x, Q2);
  xcbar  = cSymmetricSave ? xc : pdf->xfxQ2(-4, x, Q2);
  xbbar  = bSymmetricSave ? xb : pdf->xfxQ2(-5, x, Q2);
  xgamma = pdf->xfxQ2(22, x, Q2);

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// Calculate uncertainties using the LHAPDF prescription.

void LHAPDF6::calcPDFEnvelope(int idNow, double xNow, double Q2NowIn,
  int valSea) {
  if (!isSet) return;

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  double x1 = (xNow < xMin && !extrapol) ? xMin : xNow;
  if (x1 > xMax) x1 = xMax;
  double Q2Now = (Q2NowIn < q2Min) ? q2Min : Q2NowIn;
  if (Q2Now > q2Max) Q2Now = q2Max;

  // Loop over the members.
  vector<double> xfCalc(pdfs.size());
  for(int iMem = 0; iMem < pdfs.size(); ++iMem) {
    if (valSea==0 || (idNow != 1 && idNow != 2)) {
      xfCalc[iMem] = pdfs[iMem]->xfxQ2(idNow, x1, Q2Now);
    } else if (valSea==1 && (idNow == 1 || idNow == 2 )) {
      xfCalc[iMem] = pdfs[iMem]->xfxQ2(idNow, x1, Q2Now) -
        pdfs[iMem]->xfxQ2(-idNow, x1, Q2Now);
    } else if (valSea==2 && (idNow == 1 || idNow == 2 )) {
      xfCalc[iMem] = pdfs[iMem]->xfxQ2(-idNow, x1, Q2Now);
    }
  }

  // Calculate the uncertainty.
  ::LHAPDF::PDFUncertainty xfErr = pdfs.info.uncertainty(xfCalc);
  pdfEnvelope.centralPDF = xfErr.central;
  pdfEnvelope.errplusPDF = xfErr.errplus;
  pdfEnvelope.errminusPDF = xfErr.errminus;
  pdfEnvelope.errsymmPDF = xfErr.errsymm;
  pdfEnvelope.scalePDF = xfErr.scale;
}

//--------------------------------------------------------------------------

// Calculate uncertainties using the LHAPDF prescription.

void LHAPDF6::calcPDFEnvelope(pair<int,int> idNows, pair<double,double> xNows,
  double Q2NowIn, int valSea) {
  if (!isSet) return;

  // Freeze at boundary value if PDF is evaluated outside the fit region.
  double x1 = (xNows.first < xMin && !extrapol) ? xMin : xNows.first;
  if (x1 > xMax) x1 = xMax;
  double x2 = (xNows.second < xMin && !extrapol) ? xMin : xNows.second;
  if (x2 > xMax) x2 = xMax;
  double Q2Now = (Q2NowIn < q2Min) ? q2Min : Q2NowIn;
  if (Q2Now > q2Max) Q2Now = q2Max;

  // Loop over the members.
  vector<double> xfCalc(pdfs.size());
  pdfEnvelope.pdfMemberVars.resize(pdfs.size());
  for(int iMem = 0; iMem < pdfs.size(); ++iMem) {
    if        (valSea == 0 || (idNows.first != 1 && idNows.first != 2 ) ) {
      xfCalc[iMem] = pdfs[iMem]->xfxQ2(idNows.first, x1, Q2Now);
    } else if (valSea == 1 && (idNows.first == 1 || idNows.first == 2)) {
      xfCalc[iMem] = pdfs[iMem]->xfxQ2(idNows.first, x1, Q2Now)
        - pdfs[iMem]->xfxQ2(-idNows.first, x1, Q2Now);
    } else if (valSea == 2 && (idNows.first == 1 || idNows.first == 2 )) {
      xfCalc[iMem] = pdfs[iMem]->xfxQ2(-idNows.first, x1, Q2Now);
    }
    xfCalc[iMem] = max(0.0, xfCalc[iMem]);
    if        (valSea == 0 || (idNows.second != 1 && idNows.second != 2)) {
      xfCalc[iMem] /= max
        (PDFMINVALUE, pdfs[iMem]->xfxQ2(idNows.second, x2, Q2Now));
    } else if (valSea == 1 && (idNows.second == 1 || idNows.second == 2 )) {
      xfCalc[iMem] /= max
        (pdfs[iMem]->xfxQ2(idNows.second, x2, Q2Now) - pdfs[iMem]->xfxQ2
         (-idNows.second, x2, Q2Now), PDFMINVALUE);
    } else if (valSea == 2 && (idNows.second == 1 || idNows.second == 2 )) {
      xfCalc[iMem] /= max
        (pdfs[iMem]->xfxQ2(-idNows.second, x2, Q2Now), PDFMINVALUE);
    }
    pdfEnvelope.pdfMemberVars[iMem] = xfCalc[iMem];
  }

  // Calculate the uncertainty.
  ::LHAPDF::PDFUncertainty xfErr = pdfs.info.uncertainty(xfCalc);
  pdfEnvelope.centralPDF = xfErr.central;
  pdfEnvelope.errplusPDF = xfErr.errplus;
  pdfEnvelope.errminusPDF = xfErr.errminus;
  pdfEnvelope.errsymmPDF = xfErr.errsymm;
  pdfEnvelope.scalePDF = xfErr.scale;

}

//--------------------------------------------------------------------------

// Declare the plugin.

PYTHIA8_PLUGIN_CLASS(PDF, LHAPDF6, false, false, false)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_LHAPDF6_H
