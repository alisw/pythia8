// main202.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: parton distribution; LHAPDF; minimum bias; tuning

// Studies of hadron-level and parton-level minimum-bias quantities,
// comparing the internal default PDF with an external one from LHAPDF.
// Major differences indicate the need for major retuning, e.g. pT0Ref.

#include "Pythia8/Pythia.h"
#include "Pythia8/Plugins.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Machine: 1 = Tevatron, 2 = LHC. Statistics.
  int machine = 1;
  int nEvent  = 10000;

  // Results as inline histograms or pyplot ones.
  bool usepyplot = true;

  // Select new PDF set; obsolete LHAPDF5 file name conventions.
  //string pdfSet = "LHAPDF5:cteq5l.LHgrid";
  //string pdfSet = "LHAPDF5:cteq61.LHpdf";
  //string pdfSet = "LHAPDF5:cteq61.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2004nlo.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2001lo.LHgrid";

  // Select new PDF set; LHAPDF6 file name conventions.
  // (Bad/unoptimized choice, to illustrate that the PDF matters.)
  string pdfSet = "LHAPDF6:PDF4LHC15_nlo_asvar";

  // Histograms for hadron-level quantities.
  double nMax = (machine == 1) ? 199.5 : 599.5;
  Hist nChargedOld("n_charged old PDF", 100, -0.5, nMax);
  Hist nChargedNew("n_charged new PDF", 100, -0.5, nMax);
  Hist nChargedRat("n_charged new/old PDF", 100, -0.5, nMax);
  Hist ySpecOld("y charged distribution old PDF", 100, -10., 10.);
  Hist ySpecNew("y charged distribution new PDF", 100, -10., 10.);
  Hist ySpecRat("y charged distribution new/old PDF", 100, -10., 10.);
  Hist pTSpecOld("pT charged distribution old PDF", 100, 0., 20.);
  Hist pTSpecNew("pT charged distribution new PDF", 100, 0., 20.);
  Hist pTSpecRat("pT charged distribution new/old PDF", 100, 0., 20.);
  Hist avgPTnChOld("<pT>(n_charged) old PDF", 100, -0.5, nMax);
  Hist avgPTnChNew("<pT>(n_charged) new PDF", 100, -0.5, nMax);
  Hist avgPTnChRat("<pT>(n_charged) new/old PDF", 100, -0.5, nMax);

  // Histograms for parton-level quantities.
  Hist xDistOld("MPI log(x) distribution old PDF", 80, -8., 0.);
  Hist xDistNew("MPI log(x) distribution new PDF", 80, -8., 0.);
  Hist xDistRat("MPI log(x) distribution new/old PDF", 80, -8., 0.);
  Hist pTDistOld("MPI pT (=Q) distribution old PDF", 100, 0., 20.);
  Hist pTDistNew("MPI pT (=Q) distribution new PDF", 100, 0., 20.);
  Hist pTDistRat("MPI pT (=Q) distribution new/old PDF", 100, 0., 20.);

  // PDF path.
  string pdfPath;

  // Loop over one default run and one with new PDF.
  for (int iRun = 0; iRun < 2; ++iRun) {

    // Get starting time in seconds.
    clock_t tBegin = clock();

    // Generator.
    Pythia pythia;
    Event& event = pythia.event;
    pdfPath = pythia.settings.word("xmlPath") + "../pdfdata";

    // Generate minimum-bias events, with or without double diffraction.
    pythia.readString("SoftQCD:nonDiffractive = on");
    //pythia.readString("SoftQCD:doubleDiffractive = on");

    // Alternatively generate QCD jet events, above some threshold.
    //pythia.readString("HardQCD:all = on");
    //pythia.readString("PhaseSpace:pTHatMin = 50.");

    // Reduce output.
    pythia.readString("Next:numberShowEvent = 0");

    // In second run pick new PDF set.
    if (iRun == 1) {
      pythia.readString("PDF:pSet = " + pdfSet);

      // Need to change at least pT0Ref depending on choice of PDF.
      // One possibility: retune to same <n_charged>.
      //pythia.readString("MultipartonInteractions:pT0Ref = 2.17");
    }

    // Allow extrapolation of PDF's beyond x and Q2 boundaries, at own risk.
    // Default behaviour is to freeze PDF's at boundaries.
    pythia.readString("PDF:extrapolate = on");

    // Tevatron/LHC initialization.
    double eCM =  (machine == 1) ? 1960. : 13000.;
    pythia.settings.parm("Beams:eCM", eCM);
    if (machine == 1) pythia.readString("Beams:idB = -2212");
    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;


    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events.  Skip if error.
      if (!pythia.next()) continue;

      // Statistics on multiplicity and pT.
      int    nCh   = 0;
      double pTsum = 0.;
      for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) {
        ++nCh;
        pTsum += event[i].pT();

        // Fill histograms for charged y and pT spectra.
        if (iRun == 0) {
          ySpecOld.fill( event[i].y() );
          pTSpecOld.fill( event[i].pT() );
        } else {
          ySpecNew.fill( event[i].y() );
          pTSpecNew.fill( event[i].pT()  );
        }
      }

      // Fill histograms for summed quantities.
      if (iRun == 0) {
        nChargedOld.fill( nCh );
        avgPTnChOld.fill( nCh, pTsum / max(1, nCh) );
      } else {
        nChargedNew.fill( nCh );
        avgPTnChNew.fill( nCh, pTsum / max(1, nCh) );
      }

      // Loop through event record and fill x of all incoming partons.
      for (int i = 1; i < event.size(); ++i)
      if (event[i].status() == -21 || event[i].status() == -31) {
        double x = 2. * event[i].e() / eCM;
        if (iRun == 0) xDistOld.fill( log10(x) );
        else           xDistNew.fill( log10(x) );
      }

      // Loop through multiparton interactions list and fill pT of all MPI's.
      for (int i = 0; i < pythia.info.nMPI(); ++i) {
        double pT = pythia.info.pTMPI(i);
        if (iRun == 0) pTDistOld.fill( pT );
        else           pTDistNew.fill( pT );
      }

    // End of event loop.
    }

    // Statistics.
    pythia.readString("Stat:showPartonLevel = on");
    pythia.stat();

    // Get finishing time in seconds. Print used time.
    clock_t tEnd = clock();
    double tUsed = double(tEnd - tBegin) / double(CLOCKS_PER_SEC);
    cout << "\n This subrun took " << tUsed << " seconds \n" << endl;

  // End of loop over two runs.
  }

  // Form <pT>(n_charged) ratios.
  avgPTnChOld /= nChargedOld;
  avgPTnChNew /= nChargedNew;

  // Inline histogram printout.
  if (!usepyplot) {

    // Take ratios of new to old distributions.
    nChargedRat  = nChargedNew / nChargedOld;
    ySpecRat     = ySpecNew    / ySpecOld;
    pTSpecRat    = pTSpecNew    / pTSpecOld;
    avgPTnChRat  = avgPTnChNew / avgPTnChOld;
    xDistRat     = xDistNew    / xDistOld;
    pTDistRat    = pTDistNew   / pTDistOld;

    // Print histograms.
    cout << nChargedOld << nChargedNew << nChargedRat
         << ySpecOld    << ySpecNew    << ySpecRat
         << pTSpecOld   << pTSpecNew   << pTSpecRat
         << avgPTnChOld << avgPTnChNew << avgPTnChRat
         << xDistOld    << xDistNew    << xDistRat
         << pTDistOld   << pTDistNew   << pTDistRat;

  // Write Python code that can generate a PDF file with the distributions.
  } else {
    HistPlot hpl("plot202");
    hpl.frame( "fig202", "Charged multiplicity", "$n_{\\mathrm{charged}}$",
      "Rate");
    hpl.add( nChargedOld, "-,blue", "default");
    hpl.add( nChargedNew, "--,red", "PDF4LHC15_nlo");
    hpl.plot();
    hpl.frame( "", "Charged rapidity", "$y$", "Rate");
    hpl.add( ySpecOld, "-,blue", "default");
    hpl.add( ySpecNew, "--,red", "PDF4LHC15_nlo");
    hpl.plot();
    hpl.frame( "", "Charged transverse momentum", "$p_{\\perp}$", "Rate");
    hpl.add( pTSpecOld, "-,blue", "default");
    hpl.add( pTSpecNew, "--,red", "PDF4LHC15_nlo");
    hpl.plot( 0., 20., 0.1, 1e6,  true, false);
    hpl.frame( "", "Average charged transverse momentum",
      "$n_{\\mathrm{charged}}$", "$\\langle p_{\\perp} \\rangle$");
    hpl.add( avgPTnChOld, "-,blue", "default");
    hpl.add( avgPTnChNew, "--,red", "PDF4LHC15_nlo");
    hpl.plot();
  }

  // Done.
  return 0;
}
