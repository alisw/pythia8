// main94.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Dag Gillberg <dag.gillberg@cern.ch>

// Keywords: root; event display

// This is a program that uses ROOT to visualize the particles produced by
// Pythia. Particles are drawn in (y, phi)-space to depict the E/p flow.
// A pdf file is produced with multiple pages showing WH->qqbb events.

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"

using namespace Pythia8;

//==========================================================================

// Draws text to the canvas using non-direct coordinates:
// (x, y)=(0, 0) -> lower-left, (1, 1) -> upper-right.

void drawText(double x, double y, TString txt, int col= kBlack,
  double tsize = 0.032, int align = 11) {
  static auto tex = new TLatex();
  tex->SetTextColor(col);
  tex->SetTextSize(tsize);
  tex->SetTextFont(42);
  tex->SetNDC();
  tex->SetTextAlign(align);
  tex->DrawLatex(x, y, txt);
}

//==========================================================================

// Draws right-justified text, as above.

void drawTextR(double x, double y, TString txt, int col = kBlack,
  double tsize = 0.032) {
  drawText(x, y, txt, col, tsize, 31);
}

//==========================================================================

// Draw a particle in (y, phi)-space with a symbol and label.

void drawParticle(double x, double y, TString txt, int col = kGray + 1,
  double tsize = 0.03) {
  static auto m = new TMarker();
  m->SetMarkerColor(col);
  m->SetMarkerStyle(20);
  m->SetMarkerSize(1.2);
  m->DrawMarker(x, y);

  // Modify Pythia's default particle label to suit ROOT's TLatex format.
  if (txt.Contains("bar")) txt = "#bar{" + txt.ReplaceAll("bar", "") + "}";
  txt.ReplaceAll("pi", "#pi");
  txt.ReplaceAll("+", "^{+}");
  txt.ReplaceAll("-", "^{-}");
  txt.ReplaceAll("0", "^{0}");
  txt.ReplaceAll("_L", "_{L}");
  txt.ReplaceAll("gamma", "#gamma");
  txt.ReplaceAll("Lambda", "#Lambda");
  txt.ReplaceAll("Delta", "#Delta");
  txt.ReplaceAll("Sigma", "#Sigma");
  txt.ReplaceAll("rho", "#rho");
  txt.ReplaceAll("omega", "#omega");
  txt.ReplaceAll("eta", "#eta");
  txt.ReplaceAll("_S", "_{S}");
  static auto tex = new TLatex();
  tex->SetTextColor(col);
  tex->SetTextSize(tsize);
  tex->SetTextFont(42);
  tex->SetTextAlign(22);
  tex->DrawLatex(x, y - 0.2, txt);
}

//==========================================================================

// Draw a single particle.

void drawParticle(const Particle &p) {
  int col =
    // System.
    p.idAbs() == 90 ? kGray    :
    // Light quarks.
    p.idAbs() <=  4 ? kBlue    :
    // b quarks.
    p.idAbs() ==  5 ? kGreen + 2 :
    p.idAbs() == 21 ? kGray + 1  :
    p.isHadron()    ? kBlue    : kRed;
  drawParticle(p.y(), p.phi(), p.name().c_str(), col);
}

//==========================================================================

// Draw a line.

void drawLine(double x1, double y1, double x2, double y2, int col, int style) {
  static auto line = new TLine();
  line->SetLineColor(col);
  line->SetLineStyle(style);
  line->DrawLine(x1, y1, x2, y2);
}

//==========================================================================

// Example main program to draw an event display.

int main() {

  // Adjust ROOTs display options.
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  // Tick marks on top and RHS.
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "x");
  gStyle->SetTickLength(0.02, "y");

  // Define the canvas.
  auto can = new TCanvas();
  double x = 0.06, y = 0.96;
  // left-right-bottom-top.
  can->SetMargin(x, 0.02, 0.08, 0.05);

  // To draw axis frame, I find it easiest to create an empty histogram.
  double yMax = 5, phiMax = 3.4;
  // Adding a bit of margin around pi
  auto axis = new TH2D( "", ";Rapidity #it{y};Azimuth #it{#phi}", 1, -yMax,
    yMax, 1, -phiMax, phiMax);
  axis->GetYaxis()->SetTitleOffset(0.8);

  // Name of output pdf file.
  TString pdf = "main94.pdf";
  can->Print(pdf + "[");
  // Turn off loads of TCanvas::Print Info messages.
  // Current canvas added to pdf file.
  gErrorIgnoreLevel = 4000;

  // Configure and describe the process.
  // Description uses ROOT's TLatex notation
  TString desc = "#it{pp} #rightarrow #it{WH} "
    "#rightarrow #it{q#bar{q}b#bar{b}},  #sqrt{#it{s}} = 13.6 TeV";

  // Generator.
  Pythia pythia;

  // PYTHIA setup. Process selection. LHC initialization.
  pythia.readString("Beams:eCM = 13600.");
  pythia.readString("HiggsSM:ffbar2HW = on");
  // Force H->bb decays and hadronic W decays.
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfAny = 5");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = 1 2 3 4 5");
  pythia.init();

  // List of status codes to draw with associated descriptions.
  // Descriptions from "Particle Properties" in the Pythia8 manual.
  std::map<int, TString> statusCodeMap;
  //statusCodeMap[11] = "Beam particles";
  statusCodeMap[21] = "Hardest subprocess";
  statusCodeMap[31] = "Particles of subsequent subprocesses";
  statusCodeMap[41] = "Initial-state radiation";
  statusCodeMap[51] = "Final-state radiation";
  statusCodeMap[61] = "Beam remnant treatment";
  statusCodeMap[71] = "Preparation of hadronization";
  statusCodeMap[81] = "Primary hadrons";
  statusCodeMap[91] = "Decay products";

  // Event loop.
  for (int iEvent = 0; iEvent < 5; ++iEvent) {

    // Generate event. (Skip to next if pythia.next() returns false = error.)
    if (!pythia.next()) continue;

    // Clear and draw an empty canvas to paint on.
    axis->Reset();
    axis->Draw();
    drawLine( -yMax, TMath::Pi(), yMax, TMath::Pi(), kGray, 7);
    drawLine( -yMax, -TMath::Pi(), yMax, -TMath::Pi(), kGray, 7);

    // 1. Draw the hard process.
    for (int i = 0; i < pythia.process.size(); ++i) {
      auto &p = pythia.process[i];
      drawParticle(p);
      // printf("(id,pT,y,phi,status) = (%3i,%5.1f,%5.2f,%5.2f,%i)\n", p.id(),
      //   p.pT(), p.y(), p.phi(), p.status());
    }
    // Text.
    drawText( x, y, desc);
    drawTextR( 0.98, y, "Hard process");

    // Redraw the axis to make sure they are not covered by points.
    axis->Draw("axis same");
    can->Print(pdf);

    // 2. Draw intermediate particles.
    for (auto statusCode:statusCodeMap) {
      // Clear and re-draw the axis. Draw dashed-lines at +/- pi.
      axis->Reset();
      axis->Draw();
      drawLine( -yMax, TMath::Pi(), yMax, TMath::Pi(), kGray, 7);
      drawLine( -yMax, -TMath::Pi(), yMax, -TMath::Pi(), kGray, 7);
      int max = -statusCode.first, min = max-8;
      // printf("Range [%i,%i] %s\n", min, max, statusCode.second.Data());
      for (int i = 0; i < pythia.event.size(); ++i) {
        auto &p = pythia.event[i];
        // Impose desired status code(s) and rapidity intervals.
        if ( p.status() < min || p.status() > max) continue;
        if ( std::abs(p.y()) > yMax ) continue;
        drawParticle(p);
      }
      // Text.
      drawText( x, y, desc);
      drawTextR( 0.98, y, statusCode.second + Form(", status #in  [%.i,%.i]",
          min, max));
      // Redraw the axis to make sure they are not covered by points.
      axis->Draw("axis same");
      can->Print(pdf);
    }

    // 3. Draw the final hadrons.
    axis->Reset();
    axis->Draw();
    drawLine( -yMax, TMath::Pi(), yMax, TMath::Pi(), kGray, 7);
    drawLine( -yMax, -TMath::Pi(), yMax, -TMath::Pi(), kGray, 7);
    for (int i = 0; i < pythia.event.size(); ++i) {
      auto &p = pythia.event[i];
      if (not p.isFinal()) continue;
      if (std::abs(p.y())>yMax) continue;
      drawParticle(p);
    }
    drawText( x, y, desc);
    drawTextR( 0.98, y, "Stable particles");
    axis->Draw("axis same");
    can->Print(pdf);
  }

  // Close the pdf
  can->Print(pdf + "]");
  printf( "\nProduced %s\n\n", pdf.Data());

  // Done.
  return 0;
}
