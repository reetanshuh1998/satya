#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TMath.h"

void plot_fig2_3() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kTemperatureMap);

    TFile *f = TFile::Open("TaggedNeutron-DIS-EicC.root");
    if (!f) return;
    TTree *tree = (TTree*)f->Get("tree");

    double xB, Q2, t, xpi, d4sigma, W2, MX2;
    TLorentzVector *elec_out = 0;
    TLorentzVector *neut_out = 0;

    tree->SetBranchAddress("xB", &xB);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("xpi", &xpi);
    tree->SetBranchAddress("d4sigma", &d4sigma);
    tree->SetBranchAddress("W2", &W2);
    tree->SetBranchAddress("MX2", &MX2);
    tree->SetBranchAddress("elec_out", &elec_out);
    tree->SetBranchAddress("neut_out", &neut_out);

    // Normalization factor: Rate (Hz)
    double V = 1.0 * 49.0 * 0.998 * 0.99; // Δ|t|=0.99, ΔQ²=49, ΔxL=0.998 (0.001→0.999), ΔxB≈0.99
    double N_gen = tree->GetEntries();
    double L = 1.0; // nb^-1 s^-1 (10^33 cm^-2 s^-1)
    double norm = L * V / N_gen;

    TH2D *hFig2a = new TH2D("hFig2a", "Event Rate;x_{B};Q^{2} (GeV^{2})", 100, 0, 0.3, 100, 1, 50);
    TH2D *hFig2b = new TH2D("hFig2b", "Event Rate;x_{#pi};|t| (GeV^{2})", 100, 0, 1, 100, 0, 1);

    TH2D *hFig3a = new TH2D("hFig3a", ";#eta_{e};p_{e} (GeV/c)", 100, -3.5, 3.5, 100, 0, 5);
    TH2D *hFig3b = new TH2D("hFig3b", ";#theta_{n} (Deg);p_{n} (GeV/c)", 100, 0, 6, 100, 12, 20);

    for (Long64_t i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        
        // PHYSICAL CUT: x_pi must be < 1
        if (xpi >= 1.0) continue;
        
        double weight = d4sigma * norm;
        hFig2a->Fill(xB, Q2, weight);
        hFig2b->Fill(xpi, -t, weight);
        
        // Figure 3
        hFig3a->Fill(elec_out->Eta(), elec_out->P(), weight);
        
        // Neutron polar angle relative to PROTON beam
        TVector3 pBeam_v3(sin(0.05)*20, 0, cos(TMath::Pi()-0.05)*20); // Approximation
        // Actually, better to get it from the boost vector or beam vector
        // But for EicC parameters:
        double theta_rel = neut_out->Vect().Angle(pBeam_v3) * 180.0 / TMath::Pi();
        
        hFig3b->Fill(theta_rel, neut_out->P(), weight);
    }

    TCanvas *c2 = new TCanvas("c2", "Figure 2", 1400, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig2a->GetZaxis()->SetTitle("Rate (Hz)");
    hFig2a->GetZaxis()->SetTitleOffset(1.2);
    hFig2a->GetZaxis()->SetLabelSize(0.04);
    hFig2a->GetZaxis()->SetRangeUser(1e-10, 1.0);
    hFig2a->Draw("COLZ");
    
    c2->cd(2);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig2b->GetZaxis()->SetTitle("Rate (Hz)");
    hFig2b->GetZaxis()->SetTitleOffset(1.2);
    hFig2b->GetZaxis()->SetLabelSize(0.04);
    hFig2b->GetZaxis()->SetRangeUser(1e-10, 1.0);
    hFig2b->Draw("COLZ");
    c2->SaveAs("Figure2.png");

    TCanvas *c3 = new TCanvas("c3", "Figure 3", 1400, 600);
    c3->Divide(2, 1);
    c3->cd(1);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig3a->GetZaxis()->SetTitle("Rate (Hz)");
    hFig3a->GetZaxis()->SetTitleOffset(1.2);
    hFig3a->GetZaxis()->SetLabelSize(0.04);
    hFig3a->Draw("COLZ");
    
    c3->cd(2);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig3b->GetZaxis()->SetTitle("Rate (Hz)");
    hFig3b->GetZaxis()->SetTitleOffset(1.2);
    hFig3b->GetZaxis()->SetLabelSize(0.04);
    hFig3b->Draw("COLZ");
    
    // Add 50 mrad line
    TLine *l50 = new TLine(2.865, 12, 2.865, 20); // 50 mrad in degrees
    l50->SetLineColor(kRed);
    l50->SetLineStyle(2);
    l50->SetLineWidth(2);
    l50->Draw();
    
    c3->SaveAs("Figure3.png");
}
