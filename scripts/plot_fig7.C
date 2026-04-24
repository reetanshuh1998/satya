#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"

void plot_fig7() {
    gStyle->SetOptStat(0);

    // Load Signal (Pion Study - Neutron tagged)
    TFile *fSig = TFile::Open("../tagged-neutron-DIS/TaggedNeutron-DIS-EicC.root");
    if(!fSig) { cout<<"Signal file not found!"<<endl; return; }
    TTree *tSig = (TTree*)fSig->Get("tree");

    // Load Background (DIS - Neutron sample)
    TFile *fBkg = TFile::Open("../DIS_Background_EicC.root");
    if(!fBkg) { cout<<"Background file not found!"<<endl; return; }
    TTree *tBkg = (TTree*)fBkg->Get("tree");

    double xL, d4sigma, W2, MX2;
    TLorentzVector *neut_out = 0;
    tSig->SetBranchAddress("xL", &xL);
    tSig->SetBranchAddress("d4sigma", &d4sigma);
    tSig->SetBranchAddress("W2", &W2);
    tSig->SetBranchAddress("MX2", &MX2);
    tSig->SetBranchAddress("neut_out", &neut_out);

    double xL_bkg, eta_bkg, pT_bkg, W2_bkg, MX2_bkg;
    int id_bkg;
    tBkg->SetBranchAddress("xL", &xL_bkg);
    tBkg->SetBranchAddress("eta", &eta_bkg);
    tBkg->SetBranchAddress("pT", &pT_bkg);
    tBkg->SetBranchAddress("W2", &W2_bkg);
    tBkg->SetBranchAddress("MX2", &MX2_bkg);
    tBkg->SetBranchAddress("id", &id_bkg);

    TH1D *hEtaSig = new TH1D("hEtaSig", ";Neutron #eta;Events", 50, 2, 10);
    TH1D *hXLSig_pT03 = new TH1D("hXLSig_pT03", ";x_{L};Events", 50, 0.0, 1);
    TH1D *hXLSig_pT05 = new TH1D("hXLSig_pT05", ";x_{L};Events", 50, 0.0, 1);
    TH1D *hEtaBkg = new TH1D("hEtaBkg", ";Neutron #eta;Events", 50, 2, 10);
    TH1D *hXLBkg_pT03 = new TH1D("hXLBkg_pT03", ";x_{L};Events", 50, 0.0, 1);
    TH1D *hXLBkg_pT05 = new TH1D("hXLBkg_pT05", ";x_{L};Events", 50, 0.0, 1);

    // Normalization
    double V_sig = 1.0 * 49.0 * 0.998 * 0.99; // Δ|t|=0.99, ΔQ²=49, ΔxL=0.998 (0.001→0.999), ΔxB≈0.99
    double N_gen_sig = tSig->GetEntries();
    double L = 1.0; 
    double norm_sig = L * V_sig / N_gen_sig;

    double sigma_dis = 1000.0; 
    double norm_bkg = L * sigma_dis / 300000.0;

    // Fill Signal
    for (Long64_t i=0; i<tSig->GetEntries(); i++) {
        tSig->GetEntry(i);
        double weight = d4sigma * norm_sig;
        double eta = neut_out->Eta();
        double pT = neut_out->Pt();

        // (a) Cuts for eta plot
        if (xL > 0.75 && MX2 > 0.5*0.5 && W2 > 4.0) {
            hEtaSig->Fill(eta, weight);
        }

        // (b) Cuts for xL plot: two cases to compare pT threshold
        if (eta > 5.0 && pT < 0.3 && MX2 > 0.5*0.5 && W2 > 4.0) {
            hXLSig_pT03->Fill(xL, weight);
        }
        if (eta > 5.0 && pT < 0.5 && MX2 > 0.5*0.5 && W2 > 4.0) {
            hXLSig_pT05->Fill(xL, weight);
        }
    }

    // Fill Background
    for (Long64_t i=0; i<tBkg->GetEntries(); i++) {
        tBkg->GetEntry(i);
        if (abs(id_bkg) != 2112) continue; // Only neutrons
        
        // (a) Cuts for eta plot
        if (xL_bkg > 0.75 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hEtaBkg->Fill(eta_bkg, norm_bkg);
        }

        // (b) Cuts for xL plot: two cases to compare pT threshold
        if (eta_bkg > 5.0 && pT_bkg < 0.3 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hXLBkg_pT03->Fill(xL_bkg, norm_bkg);
        }
        if (eta_bkg > 5.0 && pT_bkg < 0.5 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hXLBkg_pT05->Fill(xL_bkg, norm_bkg);
        }
    }

    // Styling
    hEtaSig->SetLineColor(kOrange+7);
    hEtaSig->SetLineWidth(3);
    hEtaBkg->SetLineColor(kBlue+1);
    hEtaBkg->SetLineWidth(3);

    hXLSig_pT03->SetLineColor(kOrange+7);
    hXLSig_pT03->SetLineWidth(3);
    hXLBkg_pT03->SetLineColor(kBlue+1);
    hXLBkg_pT03->SetLineWidth(3);

    hXLSig_pT05->SetLineColor(kOrange+7);
    hXLSig_pT05->SetLineWidth(3);
    hXLBkg_pT05->SetLineColor(kBlue+1);
    hXLBkg_pT05->SetLineWidth(3);

    // ---- Figure 7 with pT < 0.3 GeV (matches paper body text and caption) ----
    TCanvas *c7 = new TCanvas("c7", "Figure 7 (pT < 0.3 GeV)", 1400, 600);
    c7->Divide(2, 1);
    
    c7->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hEtaSig->GetXaxis()->SetTitleSize(0.05);
    hEtaSig->GetYaxis()->SetTitleSize(0.05);
    double maxEta = max(hEtaSig->GetMaximum(), hEtaBkg->GetMaximum()) * 1.2;
    hEtaSig->SetMaximum(maxEta);
    hEtaSig->Draw("HIST");
    hEtaBkg->Draw("HIST SAME");
    TLine *lEta = new TLine(5, 0, 5, hEtaSig->GetMaximum()*1.1);
    lEta->SetLineStyle(2);
    lEta->SetLineWidth(2);
    lEta->Draw();
    TLegend *leg1 = new TLegend(0.18, 0.75, 0.45, 0.88);
    leg1->AddEntry(hEtaBkg, "DIS", "l");
    leg1->AddEntry(hEtaSig, "Sullivan", "l");
    leg1->SetBorderSize(1);
    leg1->Draw();

    c7->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hXLSig_pT03->GetXaxis()->SetTitleSize(0.05);
    hXLSig_pT03->GetYaxis()->SetTitleSize(0.05);
    double maxXL_03 = max(hXLSig_pT03->GetMaximum(), hXLBkg_pT03->GetMaximum()) * 1.2;
    hXLSig_pT03->SetMaximum(maxXL_03);
    hXLSig_pT03->Draw("HIST");
    hXLBkg_pT03->Draw("HIST SAME");
    TLine *lXL_03 = new TLine(0.75, 0, 0.75, hXLSig_pT03->GetMaximum()*1.1);
    lXL_03->SetLineStyle(2);
    lXL_03->SetLineWidth(2);
    lXL_03->Draw();
    TLegend *leg2 = new TLegend(0.18, 0.75, 0.45, 0.88);
    leg2->AddEntry(hXLBkg_pT03, "DIS", "l");
    leg2->AddEntry(hXLSig_pT03, "Sullivan", "l");
    leg2->SetBorderSize(1);
    leg2->Draw();

    c7->SaveAs("Figure7_pT03.png");

    // ---- Figure 7 with pT < 0.5 GeV (matches in-figure label in paper) ----
    TCanvas *c7b = new TCanvas("c7b", "Figure 7 (pT < 0.5 GeV)", 1400, 600);
    c7b->Divide(2, 1);

    c7b->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hEtaSig->GetXaxis()->SetTitleSize(0.05);
    hEtaSig->GetYaxis()->SetTitleSize(0.05);
    hEtaSig->SetMaximum(maxEta);
    hEtaSig->Draw("HIST");
    hEtaBkg->Draw("HIST SAME");
    TLine *lEta2 = new TLine(5, 0, 5, hEtaSig->GetMaximum()*1.1);
    lEta2->SetLineStyle(2);
    lEta2->SetLineWidth(2);
    lEta2->Draw();
    TLegend *leg3 = new TLegend(0.18, 0.75, 0.45, 0.88);
    leg3->AddEntry(hEtaBkg, "DIS", "l");
    leg3->AddEntry(hEtaSig, "Sullivan", "l");
    leg3->SetBorderSize(1);
    leg3->Draw();

    c7b->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hXLSig_pT05->GetXaxis()->SetTitleSize(0.05);
    hXLSig_pT05->GetYaxis()->SetTitleSize(0.05);
    double maxXL_05 = max(hXLSig_pT05->GetMaximum(), hXLBkg_pT05->GetMaximum()) * 1.2;
    hXLSig_pT05->SetMaximum(maxXL_05);
    hXLSig_pT05->Draw("HIST");
    hXLBkg_pT05->Draw("HIST SAME");
    TLine *lXL_05 = new TLine(0.75, 0, 0.75, hXLSig_pT05->GetMaximum()*1.1);
    lXL_05->SetLineStyle(2);
    lXL_05->SetLineWidth(2);
    lXL_05->Draw();
    TLegend *leg4 = new TLegend(0.18, 0.75, 0.45, 0.88);
    leg4->AddEntry(hXLBkg_pT05, "DIS", "l");
    leg4->AddEntry(hXLSig_pT05, "Sullivan", "l");
    leg4->SetBorderSize(1);
    leg4->Draw();

    c7b->SaveAs("Figure7_pT05.png");
}
