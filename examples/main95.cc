// main95.cc is a part of the PYTHIA event generator.
// Copyright (C) 2024 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Dag Gillberg <dag.gillberg@cern.ch>

// Keywords: root; jets; event display

// This is a program to use ROOT to visualize different jet algoritms.
// The produced figure is used in the article "50 years of Quantum
// Chromodynamics" in celebration of the 50th anniversary of QCD (EPJC).

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TPave.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

//==========================================================================

// Hard-coded settings

// Jet and hadron pT thresholds.
// Will only show particles with pT > pTmin and |y| < yMax.
double pTmin_jet = 25;
double pTmin_hadron = 1;
double yMax = 4;

// Amount of pileup. Average number of inelastic pp collisions per event
// (=bunch-crossing). Set to zero to turn off pileup.
double mu = 60;

// Style format. Colours used by various drawn markers.
int colHS = kBlack, colPos = kRed, colNeg = kBlue;
int colNeut = kGreen + 3, colPU = kGray + 1;

using namespace Pythia8;

//==========================================================================

// Method to print descriptive text to the canvas.
// (x, y) are relative coordinates (NDC).

void drawText(double x, double y, TString txt, int align= 11,
  double tsize= 0.032) {
  static auto tex = new TLatex();
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
  tex->SetTextFont(42);
  tex->SetNDC();
  tex->DrawLatex(x, y, txt);
}

//==========================================================================

// Text to draw a marker at the (y, phi) coordinates of a particle.
// Absolute coordinates.

void drawParticleMarker(const Particle &p, int style, int col,
  double size= 1.0) {
  static auto m = new TMarker();
  m->SetMarkerStyle(style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->DrawMarker(p.y(), p.phi());
}

//==========================================================================

// Method to draw a marker+text of a particle.

void drawParticleText(const Particle &p) {
  // Draws a marker at (y, phi) of particle. Circle for parton, star
  // for boson.
  bool isParton =  (std::abs(p.id()) <= 5 || p.id() == 21);
  int col = colHS;
  drawParticleMarker( p, isParton?20:29, col, isParton?0.8:1.2);

  // Format the name-string of the particle according to ROOT's TLatex.
  // Print the text right under the marker.
  TString name = p.name();
  if (name.Contains("bar")) name = "#bar{" + name.ReplaceAll("bar", "") + "}";
  name.ReplaceAll("+", "^{+}").ReplaceAll("-", "^{-}").ReplaceAll("h0", "H");
  static auto tex = new TLatex();
  tex->SetTextSize(0.03);
  tex->SetTextFont(42);
  tex->SetTextAlign(11);
  tex->SetTextColor(col);
  tex->DrawLatex(p.y() + 0.1, p.phi() - 0.1, "#it{" + name + "}");
}

//==========================================================================

// Draws a box for text to appear.

void drawLegendBox(double x1, double y1, double x2, double y2) {
  static auto *box = new TPave(x1, y1, x2, y2, 1, "ndc");
  box->SetFillColor(kWhite);
  box->Draw();
}

//==========================================================================

// Draw a marker for legend.

void drawMarker(double x, double y, int style, int col, double size= 1.0) {
  auto m = new TMarker(x, y, style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->SetNDC(true);
  m->Draw();
}

//==========================================================================

// Example main program to vizualize jet algorithms.

int main() {

  // Adjust ROOTs default style.
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  // Tick marks on top and RHS.
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02, "x");
  gStyle->SetTickLength(0.015, "y");
  // Good with SetMax higher. 57, 91 and 104 also OK.
  gStyle->SetPalette(55);

  // Define the canvas.
  auto can = new TCanvas();
  double x = 0.06, y = 0.96;
  // Left-right-bottom-top
  can->SetMargin(x, 0.02, 0.08, 0.06);

  // Define the energy-flow histogram.
  int NyBins = 400/2, NphiBins = 314/2;
  double yMax = 4, phiMax = TMath::Pi();
  auto pTflow = new TH2D("",
    ";Rapidity #it{y};Azimuth #it{#phi};Jet #it{p}_{T} [GeV]",
    NyBins, -yMax, yMax, NphiBins, -phiMax, phiMax);
  pTflow->GetYaxis()->SetTitleOffset(0.8);
  pTflow->GetZaxis()->SetTitleOffset(1.1);

  // Name of output pdf file + open canvas for printing pages to it.
  TString pdf = "main95.pdf";
  can->Print(pdf + "[");

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  // Description of the process (using ROOT's TLatex notation).
  TString desc = "#it{pp} #rightarrow #it{WH} #rightarrow"
    " #it{q#bar{q}b#bar{b}},  #sqrt{#it{s}} = 13.6 TeV";
  pythia.readString("Beams:eCM = 13600.");
  pythia.readString("HiggsSM:ffbar2HW = on");
  // Force H->bb decays and hadronic W decays.
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfAny = 5");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = 1 2 3 4 5");
  pythia.init();

  // Pileup particles
  Pythia pythiaPU;
  pythiaPU.readString("Beams:eCM = 13600.");
  pythiaPU.readString("SoftQCD:inelastic = on");
  if (mu > 0) pythiaPU.init();

  // Setup fasjet. Create map with (key, value) = (descriptive text, jetDef).
  std::map<TString, fastjet::JetDefinition> jetDefs;
  jetDefs["Anti-#it{k_{t}} jets, #it{R} = 0.4"] = fastjet::JetDefinition(
    fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
  jetDefs["#it{k_{t}} jets, #it{R} = 0.4"] = fastjet::JetDefinition(
    fastjet::kt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
  jetDefs["Cambridge-Aachen jets,  #it{R} = 0.4"] = fastjet::JetDefinition(
    fastjet::cambridge_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
  jetDefs["Anti-#it{k_{t}} jets, #it{R} = 1.0"] = fastjet::JetDefinition(
    fastjet::antikt_algorithm, 1.0, fastjet::E_scheme, fastjet::Best);

  auto &event = pythia.event;
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;

    // Identify particles. Jets are built from all stable particles after
    // hadronization (particle-level jets).
    std::vector<Particle> VH, ptcls_hs, ptcls_pu;
    std::vector<fastjet::PseudoJet> stbl_ptcls;
    for (int i = 0; i < event.size(); ++i) {
      auto &p = event[i];
      if (p.isResonance() && p.status() == -62) VH.push_back(p);
      if (not p.isFinal()) continue;
      stbl_ptcls.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
      ptcls_hs.push_back(p);
    }

    // Should not happen!
    if (VH.size()!= 2) continue;

    // Want to show an example where one of the boson is boosted and hence
    // contained within a R=1.0 jet, and one is not.
    // The diboson system should also be fairly central.
    auto pVH = VH[0].p() + VH[1].p();
    double DR1 = event.RRapPhi(VH[0].daughter1(), VH[0].daughter2());
    double DR2 = event.RRapPhi(VH[1].daughter1(), VH[1].daughter2());
    // Central system.
    if ( std::abs(pVH.rap())>0.5 || std::abs(VH[0].phi())>2.5 ||
      std::abs(VH[1].phi())>2.5 ) continue;
    // One contained, one resolved.
    if ( (DR1<1.0 && DR2<1.0) || (DR1>1.0 && DR2>1.0) ) continue;

    // Add in ghost particles on the grid defined by the pTflow histogram.
    fastjet::PseudoJet ghost;
    double pTghost = 1e-100;
    for (int iy= 1;iy<= NyBins; ++iy) {
      for (int iphi= 1;iphi<= NphiBins; ++iphi) {
        double y   = pTflow->GetXaxis()->GetBinCenter(iy);
        double phi = pTflow->GetYaxis()->GetBinCenter(iphi);
        ghost.reset_momentum_PtYPhiM(pTghost, y, phi, 0);
        stbl_ptcls.push_back(ghost);
      }
    }

    // Add in pileup!
    int n_inel = gRandom->Poisson(mu);
    printf("Overlaying particles from %i pileup interactions!\n", n_inel);
    for (int i_pu= 0; i_pu<n_inel; ++i_pu) {
      if (!pythiaPU.next()) continue;
      for (int i = 0; i < pythiaPU.event.size(); ++i) {
        auto &p = pythiaPU.event[i];
        if (not p.isFinal()) continue;
        stbl_ptcls.push_back(
          fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
        ptcls_pu.push_back(p);
      }
    }

    can->SetLogz();
    can->SetRightMargin(0.13);
    bool first = true;
    for (auto jetDef:jetDefs) {
      fastjet::ClusterSequence clustSeq(stbl_ptcls, jetDef.second);
      auto jets = sorted_by_pt( clustSeq.inclusive_jets(pTmin_jet) );
      // Fill the pT flow.
      pTflow->Reset();
      // For each jet:
      for (auto jet:jets) {
        // For each particle:
        for (auto c:jet.constituents()) {
          if (c.pt() > 1e-50) continue;
          pTflow->Fill(c.rap(), c.phi_std(), jet.pt());
        }
      }
      pTflow->GetZaxis()->SetRangeUser(pTmin_jet/4,
        pTflow->GetBinContent(pTflow->GetMaximumBin())*4);
      // pTflow->GetZaxis()->SetRangeUser(8, 1100);
      // pTflow->GetZaxis()->SetMoreLogLabels();
      pTflow->Draw("colz");

      for (auto &p:ptcls_pu) {
        if ( std::abs(p.y()) < yMax && p.pT() > pTmin_hadron ) {
          drawParticleMarker(p, p.charge()?24:25, colPU, 0.4);
        }
      }

      // Draw the stable particles.
      for (auto &p:ptcls_hs) {
        if ( std::abs(p.y()) < yMax && p.pT() > pTmin_hadron ) {
          if (p.charge()>0) {
            drawParticleMarker(p, 5, colPos, 0.8);
          } else if (p.charge()<0) {
            drawParticleMarker(p, 5, colNeg, 0.8);
          } else {
            drawParticleMarker(p, 21, colNeut, 0.4);
            drawParticleMarker(p, 5, colNeut, 0.8);
          }
        }
      }

      // Draw the W and H.
      for (auto p:VH) {
        auto d1 = pythia.event[p.daughter1()];
        auto d2 = pythia.event[p.daughter2()];
        drawParticleText(p);
        drawParticleText(d1);
        drawParticleText(d2);
      }

      drawText(x, y, desc);
      drawText(0.87, y, jetDef.first +
        Form(", #it{p}_{T} > %.0f GeV", pTmin_jet), 31);

      // 'Hand-made' legend used to specific plot in the
      // '50 years of Quantum Chromodynamics', EPJC.
      if (first) {
        drawLegendBox(0.66, 0.67, 0.85, 0.925);
        drawText(0.715, 0.90, "Hard scatter", 12);
        drawMarker(0.68, 0.90, 20, colHS, 0.8);
        drawMarker(0.7, 0.90, 29, colHS, 1.2);

        drawText(0.675, 0.85, "Stable particles", 12);
        drawText(0.675, 0.824, "   +    #bf{#minus}    #scale[0.9]{neutral}",
          12);
        drawMarker(0.683, 0.82, 5, colPos, 0.8);
        drawMarker(0.717, 0.82, 5, colNeg, 0.8);
        drawMarker(0.75, 0.82, 21, colNeut, 0.4);
        drawMarker(0.75, 0.82, 5, colNeut, 0.8);

        drawText(0.675, 0.775, Form("Pileup  #it{#mu} = %.0f", mu), 12);
        drawText(0.675, 0.745, "   #pm    #scale[0.9]{neutral}", 12);
        drawMarker(0.683, 0.74, 24, colPU, 0.4);
        drawMarker(0.717, 0.74, 25, colPU, 0.4);
        drawText(0.70, 0.70, Form("#scale[0.8]{#it{p}_{T}^{ptcl} > %.1f GeV}",
            pTmin_hadron), 12);
      }
      first = false;
      can->Print(pdf);
    }
    break; // remove this to draw several events
  }

  // Close the pdf
  can->Print(pdf + "]");
  printf("Produced %s\n\n", pdf.Data());

  // Done.
  return 0;
}
