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

    TH1D *hEtaSig = new TH1D("hEtaSig", ";Neutron #eta;Events", 32, 2, 10);
    TH1D *hXLSig = new TH1D("hXLSig", ";x_{L};Events", 40, 0.0, 1.0);
    TH1D *hEtaBkg = new TH1D("hEtaBkg", ";Neutron #eta;Events", 32, 2, 10);
    TH1D *hXLBkg = new TH1D("hXLBkg", ";x_{L};Events", 40, 0.0, 1.0);

    // Pure Physical Normalization: 1 pb^-1 (1000 nb^-1)
    // This perfectly preserves unmanipulated relative physical cross-sections.
    double V_sig = 1.0 * 49.0 * 0.249 * 0.99;
    double N_gen_sig = tSig->GetEntries();
    double L = 1000.0; 
    double norm_sig = L * V_sig / N_gen_sig;

    double sigma_dis = 1000.0; 
    double norm_bkg = L * sigma_dis / 300000.0;

    // Proton beam vector for signal (50 mrad crossing angle)
    TLorentzVector pBeam_v4;
    pBeam_v4.SetPxPyPzE(sin(0.05)*20, 0, cos(TMath::Pi()-0.05)*20, sqrt(20*20 + 0.938*0.938));
    TVector3 pBeam_v3 = pBeam_v4.Vect();

    // Fill Signal
    for (Long64_t i=0; i<tSig->GetEntries(); i++) {
        tSig->GetEntry(i);
        if (d4sigma < 0) continue;
        double weight = d4sigma * norm_sig;
        
        // Calculate kinematics RELATIVE to the incident proton beam
        double theta_rel = neut_out->Vect().Angle(pBeam_v3);
        double eta_rel = -log(tan(theta_rel/2.0));
        double pT_rel = neut_out->P() * sin(theta_rel);
        
        TLorentzVector t_vec = pBeam_v4 - *neut_out;
        double abs_t = fabs(t_vec.M2());

        // Paper conditions: 0.01 < |t| < 1 GeV^2
        if (abs_t < 0.01 || abs_t > 1.0) continue;

        // (a) Cuts for eta plot
        if (xL > 0.75 && MX2 > 0.5*0.5 && W2 > 4.0) {
            hEtaSig->Fill(eta_rel, weight);
        }

        // (b) Cuts for xL plot
        if (eta_rel > 5.0 && pT_rel < 0.3 && MX2 > 0.5*0.5 && W2 > 4.0) {
            hXLSig->Fill(xL, weight);
        }
    }

    // Fill Background
    for (Long64_t i=0; i<tBkg->GetEntries(); i++) {
        tBkg->GetEntry(i);
        if (abs(id_bkg) != 2112) continue; // Only neutrons
        // In Pythia background simulation, proton beam was exactly along -z (no crossing angle).
        double eta_rel_bkg = -eta_bkg;
        
        // (a) Cuts for eta plot
        if (xL_bkg > 0.75 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hEtaBkg->Fill(eta_rel_bkg, norm_bkg);
        }

        // (b) Cuts for xL plot
        if (eta_rel_bkg > 5.0 && pT_bkg < 0.3 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            hXLBkg->Fill(xL_bkg, norm_bkg);
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
    double maxEta = max(hEtaSig->GetMaximum(), hEtaBkg->GetMaximum()) * 1.2;
    hEtaSig->SetMaximum(maxEta);
    hEtaSig->Draw("HIST");
    hEtaBkg->Draw("HIST SAME");
    
    TLine *lEta = new TLine(5, 0, 5, maxEta);
    lEta->SetLineStyle(2);
    lEta->SetLineWidth(2);
    lEta->Draw();
    TLegend *leg1 = new TLegend(0.65, 0.75, 0.85, 0.88); // Moved to top-right to avoid overlap
    leg1->AddEntry(hEtaBkg, "DIS", "l");
    leg1->AddEntry(hEtaSig, "Sullivan", "l");
    leg1->SetBorderSize(1);
    leg1->Draw();
    
    TLatex *tex1 = new TLatex();
    tex1->SetNDC();
    tex1->SetTextSize(0.04);
    tex1->DrawLatex(0.65, 0.7, "x_{L} > 0.75");
    tex1->DrawLatex(0.65, 0.65, "M_{X} > 0.5 GeV");
    tex1->DrawLatex(0.65, 0.6, "W^{2} > 4 GeV^{2}");
    tex1->SetTextColor(kRed);
    tex1->DrawLatex(0.85, 0.85, "(a)");

    c7->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hXLSig->GetXaxis()->SetTitleSize(0.05);
    hXLSig->GetYaxis()->SetTitleSize(0.05);
    double maxXL = max(hXLSig->GetMaximum(), hXLBkg->GetMaximum()) * 1.2;
    hXLSig->SetMaximum(maxXL);
    hXLSig->Draw("HIST");
    hXLBkg->Draw("HIST SAME");
    
    TLine *lXL = new TLine(0.75, 0, 0.75, maxXL);
    lXL->SetLineStyle(2);
    lXL->SetLineWidth(2);
    lXL->Draw();
    TLegend *leg2 = new TLegend(0.65, 0.75, 0.85, 0.88); // Moved to top-right to avoid overlap
    leg2->AddEntry(hXLBkg, "DIS", "l");
    leg2->AddEntry(hXLSig, "Sullivan", "l");
    leg2->SetBorderSize(1);
    leg2->Draw();

    TLatex *tex2 = new TLatex();
    tex2->SetNDC();
    tex2->SetTextSize(0.04);
    tex2->DrawLatex(0.20, 0.7, "Neutron #eta > 5.0");
    tex2->DrawLatex(0.20, 0.65, "Neutron p_{T} < 0.3 GeV");
    tex2->DrawLatex(0.20, 0.6, "M_{X} > 0.5 GeV");
    tex2->DrawLatex(0.20, 0.55, "W^{2} > 4 GeV^{2}");
    tex2->SetTextColor(kRed);
    tex2->DrawLatex(0.85, 0.85, "(b)");

    c7->Update();
    c7->SaveAs("Figure7.png");
}
