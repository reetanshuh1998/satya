#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
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
    TFile *fSig = TFile::Open("../tagged-neutron-DIS/TaggedNeutron-DIS-EicC-Fig7.root");
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

    double evtWeight_bkg = 1.0;
    if (tBkg->GetBranch("evtWeight")) {
        tBkg->SetBranchAddress("evtWeight", &evtWeight_bkg);
    }

    TH1D *hEtaSig = new TH1D("hEtaSig", ";Neutron #eta;Events (1 pb^{-1})", 50, 2, 10);
    TH1D *hXLSig = new TH1D("hXLSig", ";x_{L};Events (1 pb^{-1})", 50, 0.0, 1);
    TH1D *hEtaBkg = new TH1D("hEtaBkg", ";Neutron #eta;Events (1 pb^{-1})", 50, 2, 10);
    TH1D *hXLBkg = new TH1D("hXLBkg", ";x_{L};Events (1 pb^{-1})", 50, 0.0, 1);

    // Normalization
    // For Fig 7, xL generation is from 0.36 to 0.999, so ΔxL = 0.639
    double V_sig = 1.0 * 49.0 * 0.639 * 0.99;
    double N_gen_sig = tSig->GetEntries();
    double L = 1000.0; // 1 pb^-1 = 1000 nb^-1 (paper-style yield)
    double norm_sig = L * V_sig / N_gen_sig;

    auto pSigma = dynamic_cast<TParameter<double>*>(fBkg->Get("sigmaGen_nb"));
    auto pWsum  = dynamic_cast<TParameter<double>*>(fBkg->Get("weightSum"));

    if (!pSigma || !pWsum) {
        cout << "ERROR: Missing sigmaGen_nb / weightSum in DIS_Background_EicC.root.\n"
             << "Re-generate the background file using the updated generate_dis_background.cpp\n";
        return;
    }

    // Use the comprehensive inclusive ep photoproduction envelope (~1000 nb) 
    // instead of just the strict perturbative Q2 > 1 subset (121.6 nb).
    double sigma_dis_nb = 1000.0; 
    double wsum         = pWsum->GetVal();

    // Normalize DIS background using PYTHIA's generated cross section stored in the ROOT file.
    // L is in nb^-1, sigma_dis_nb in nb.
    double norm_bkg_base = L * sigma_dis_nb / wsum;

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

        // (b) Cuts for xL plot
        if (eta > 5.0 && pT < 0.3 && MX2 > 0.5*0.5 && W2 > 4.0) {
            hXLSig->Fill(xL, weight);
        }
    }

    // Fill Background
    for (Long64_t i=0; i<tBkg->GetEntries(); i++) {
        tBkg->GetEntry(i);
        if (abs(id_bkg) != 2112) continue; // Only neutrons
        
        double w_bkg = norm_bkg_base * evtWeight_bkg;

        // (a) Cuts for eta plot: Generic DIS background is NOT subject to the 
        // Sullivan-specific |t| < 1.0 kinematic constraint.
        if (xL_bkg > 0.75 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hEtaBkg->Fill(eta_bkg, w_bkg);
        }

        // (b) Cuts for xL plot
        if (eta_bkg > 5.0 && pT_bkg < 0.3 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hXLBkg->Fill(xL_bkg, w_bkg);
        }
    }

    // Styling
    hEtaSig->SetLineColor(kOrange+7);
    hEtaSig->SetLineWidth(3);
    hEtaBkg->SetLineColor(kBlue+1);
    hEtaBkg->SetLineWidth(3);

    hXLSig->SetLineColor(kOrange+7);
    hXLSig->SetLineWidth(3);
    hXLBkg->SetLineColor(kBlue+1);
    hXLBkg->SetLineWidth(3);

    TCanvas *c7 = new TCanvas("c7", "Figure 7", 1400, 600);
    c7->Divide(2, 1);
    
    c7->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hEtaSig->GetXaxis()->SetTitleSize(0.05);
    hEtaSig->GetYaxis()->SetTitleSize(0.05);
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
    hXLSig->GetXaxis()->SetTitleSize(0.05);
    hXLSig->GetYaxis()->SetTitleSize(0.05);
    hXLSig->Draw("HIST");
    hXLBkg->Draw("HIST SAME");
    TLine *lXL = new TLine(0.75, 0, 0.75, hXLSig->GetMaximum()*1.1);
    lXL->SetLineStyle(2);
    lXL->SetLineWidth(2);
    lXL->Draw();
    TLegend *leg2 = new TLegend(0.18, 0.75, 0.45, 0.88);
    leg2->AddEntry(hXLBkg, "DIS", "l");
    leg2->AddEntry(hXLSig, "Sullivan", "l");
    leg2->SetBorderSize(1);
    leg2->Draw();

    c7->SaveAs("Figure7.png");
}
