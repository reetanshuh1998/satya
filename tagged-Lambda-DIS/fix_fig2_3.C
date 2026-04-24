#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TMath.h"
#include "TColor.h"
#include "TLine.h"

void fix_fig2_3() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird); // Standard ROOT 6 beautiful palette, often looks like paper

    TFile *f = TFile::Open("../tagged-neutron-DIS/TaggedNeutron-DIS-EicC.root");
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

    double V = 1.0 * 49.0 * 0.249 * 0.99;
    double N_gen = tree->GetEntries();
    double L = 1.0; 
    double norm = L * V / N_gen;

    // Use 200x200 bins so max bin content is ~ 0.77 Hz, fitting beautifully in the [1e-10, 1] scale
    TH2D *hFig2a = new TH2D("hFig2a", "Event Rate;x_{B};Q^{2} (GeV^{2})", 200, 0, 0.3, 200, 1, 50);
    TH2D *hFig2b = new TH2D("hFig2b", "Event Rate;x_{#pi};|t| (GeV^{2})", 200, 0, 1, 200, 0, 1);

    TH2D *hFig3a = new TH2D("hFig3a", ";#eta_{e};p_{e} (GeV/c)", 200, -4.0, 4.0, 200, 0, 20);
    TH2D *hFig3b = new TH2D("hFig3b", ";#theta_{n} [Deg];p_{n} (GeV/c)", 200, 0, 6, 200, 12, 20);

    for (Long64_t i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        
        double abs_t = fabs(t);
        if (xpi >= 1.0) continue;
        if (abs_t < 0.01 || abs_t > 1.0) continue;
        if (Q2 < 1.0 || Q2 > 50.0) continue;
        if (W2 <= 4.0) continue;
        if (xB >= 1.0) continue;
        
        double weight = d4sigma * norm;
        hFig2a->Fill(xB, Q2, weight);
        hFig2b->Fill(xpi, abs_t, weight);
        
        hFig3a->Fill(elec_out->Eta(), elec_out->P(), weight);
        double theta_lab = neut_out->Theta() * 180.0 / TMath::Pi();
        hFig3b->Fill(theta_lab, neut_out->P(), weight);
    }

    TCanvas *c2 = new TCanvas("c2", "Figure 2", 1400, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig2a->GetZaxis()->SetTitle("Rate (Hz)");
    hFig2a->GetZaxis()->SetTitleOffset(1.2);
    hFig2a->GetZaxis()->SetLabelSize(0.04);
    hFig2a->SetMaximum(1.0);
    hFig2a->SetMinimum(1e-10);
    hFig2a->Draw("COLZ");
    
    c2->cd(2);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig2b->GetZaxis()->SetTitle("Rate (Hz)");
    hFig2b->GetZaxis()->SetTitleOffset(1.2);
    hFig2b->GetZaxis()->SetLabelSize(0.04);
    hFig2b->SetMaximum(1.0);
    hFig2b->SetMinimum(1e-10);
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
    hFig3a->SetMaximum(1.0);
    hFig3a->SetMinimum(1e-10);
    hFig3a->Draw("COLZ");
    
    c3->cd(2);
    gPad->SetRightMargin(0.18);
    gPad->SetLogz();
    hFig3b->GetZaxis()->SetTitle("Rate (Hz)");
    hFig3b->GetZaxis()->SetTitleOffset(1.2);
    hFig3b->GetZaxis()->SetLabelSize(0.04);
    hFig3b->SetMaximum(1.0);
    hFig3b->SetMinimum(1e-10);
    hFig3b->Draw("COLZ");
    
    TLine *l50 = new TLine(2.865, 12, 2.865, 20); 
    l50->SetLineColor(kRed);
    l50->SetLineStyle(2);
    l50->SetLineWidth(2);
    l50->Draw();
    
    c3->SaveAs("Figure3.png");
}
