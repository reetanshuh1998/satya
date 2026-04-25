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
#include "TLatex.h"

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

    TH1D *hEtaSig = new TH1D("hEtaSig", ";Neutron #eta;Events", 20, 2, 10);
    TH1D *hXLSig = new TH1D("hXLSig", ";x_{L};Events", 35, 0.0, 1.0);
    TH1D *hXLSig_pT03 = new TH1D("hXLSig_pT03", ";x_{L};Events", 35, 0.0, 1.0);
    TH1D *hXLSig_pT05 = new TH1D("hXLSig_pT05", ";x_{L};Events", 35, 0.0, 1.0);
    TH1D *hEtaBkg = new TH1D("hEtaBkg", ";Neutron #eta;Events", 20, 2, 10);
    TH1D *hXLBkg = new TH1D("hXLBkg", ";x_{L};Events", 35, 0.0, 1.0);
    TH1D *hXLBkg_pT03 = new TH1D("hXLBkg_pT03", ";x_{L};Events", 35, 0.0, 1.0);
    TH1D *hXLBkg_pT05 = new TH1D("hXLBkg_pT05", ";x_{L};Events", 35, 0.0, 1.0);

    // Pure Physical Normalization: 1 pb^-1 (1000 nb^-1)
    // For Fig 7, xL generation is from 0.36 to 0.999, so ΔxL = 0.639
    double V_sig = 1.0 * 49.0 * 0.639 * 0.99;
    double N_gen_sig = tSig->GetEntries();
    double L = 1000.0; // 1 pb^-1
    double norm_sig = L * V_sig / N_gen_sig;


    auto pSigma = dynamic_cast<TParameter<double>*>(fBkg->Get("sigmaGen_nb"));
    auto pWsum  = dynamic_cast<TParameter<double>*>(fBkg->Get("weightSum"));

    if (!pSigma || !pWsum) {
        cout << "ERROR: Missing sigmaGen_nb / weightSum in DIS_Background_EicC.root.\n"
             << "Re-generate the background file using the updated generate_dis_background.cpp\n";
        return;
    }

    double sigma_dis_nb = pSigma->GetVal();
    double wsum         = pWsum->GetVal();

    // Normalize background to integrated luminosity L (nb^-1):
    double norm_bkg_base = L * sigma_dis_nb / wsum;

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

        // (b) Cuts for xL plot: two cases to compare pT threshold
        if (eta_rel > 5.0 && MX2 > 0.5*0.5 && W2 > 4.0) {
            if (pT_rel < 0.3) {
                hXLSig->Fill(xL, weight);
                hXLSig_pT03->Fill(xL, weight);
            }
            if (pT_rel < 0.5) {
                hXLSig_pT05->Fill(xL, weight);
            }
        }
    }

    // Fill Background
    for (Long64_t i=0; i<tBkg->GetEntries(); i++) {
        tBkg->GetEntry(i);
        if (abs(id_bkg) != 2112) continue; // Only neutrons
        // In Pythia background simulation, proton beam was exactly along -z (no crossing angle).
        double eta_rel_bkg = -eta_bkg;

        double w_bkg = norm_bkg_base * evtWeight_bkg;

        // Reconstruct neutron 4-vector to calculate |t|
        double pz_bkg = pT_bkg * sinh(eta_bkg);
        double E_bkg = sqrt(pT_bkg*pT_bkg + pz_bkg*pz_bkg + 0.939565*0.939565);
        TLorentzVector neut_bkg_v4(pT_bkg, 0, pz_bkg, E_bkg); // Assume phi=0 for |t| approx
        
        // Background was generated with proton beam along -z
        TLorentzVector pBeam_bkg_v4(0, 0, -20.0, sqrt(20.0*20.0 + 0.938*0.938));
        TLorentzVector t_vec_bkg = pBeam_bkg_v4 - neut_bkg_v4;
        double abs_t_bkg = fabs(t_vec_bkg.M2());

        // (a) Cuts for eta plot: apply consistent |t| < 1.0 cut as done for signal
        if (xL_bkg > 0.75 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0 && abs_t_bkg > 0.01 && abs_t_bkg < 1.0) {
            hEtaBkg->Fill(eta_rel_bkg, w_bkg);
        }

        // (b) Cuts for xL plot: two cases to compare pT threshold
        if (eta_rel_bkg > 5.0 && MX2_bkg > 0.5*0.5 && W2_bkg > 4.0) {
            if (pT_bkg < 0.3) {
                hXLBkg->Fill(xL_bkg, w_bkg);
                hXLBkg_pT03->Fill(xL_bkg, w_bkg);
            }
            if (pT_bkg < 0.5) {
                hXLBkg_pT05->Fill(xL_bkg, w_bkg);
            }
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

    hXLSig_pT03->SetLineColor(kOrange+7);
    hXLSig_pT03->SetLineWidth(3);
    hXLBkg_pT03->SetLineColor(kBlue+1);
    hXLBkg_pT03->SetLineWidth(3);
    
    hXLSig_pT05->SetLineColor(kOrange+7);
    hXLSig_pT05->SetLineWidth(3);
    hXLBkg_pT05->SetLineColor(kBlue+1);
    hXLBkg_pT05->SetLineWidth(3);

    std::cout << "Sullivan integral (7a): " << hEtaSig->Integral() << " events" << std::endl;
    std::cout << "DIS integral (7a): " << hEtaBkg->Integral() << " events" << std::endl;
    std::cout << "Sullivan integral (7b, pT<0.3): " << hXLSig->Integral() << " events" << std::endl;
    std::cout << "DIS integral (7b, pT<0.3): " << hXLBkg->Integral() << " events" << std::endl;

    auto draw_fig7 = [&](TH1D* sig, TH1D* bkg, const char* name, const char* title, const char* pT_text) {
        TCanvas *c = new TCanvas(name, title, 1400, 600);
        c->Divide(2, 1);
        
        c->cd(1);
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
        TLegend *leg1 = new TLegend(0.65, 0.75, 0.85, 0.88);
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

        c->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        sig->GetXaxis()->SetTitleSize(0.05);
        sig->GetYaxis()->SetTitleSize(0.05);
        double maxXL = max(sig->GetMaximum(), bkg->GetMaximum()) * 1.2;
        sig->SetMaximum(maxXL);
        sig->Draw("HIST");
        bkg->Draw("HIST SAME");
        
        TLine *lXL = new TLine(0.75, 0, 0.75, maxXL);
        lXL->SetLineStyle(2);
        lXL->SetLineWidth(2);
        lXL->Draw();
        TLegend *leg2 = new TLegend(0.20, 0.75, 0.40, 0.88);
        leg2->AddEntry(bkg, "DIS", "l");
        leg2->AddEntry(sig, "Sullivan", "l");
        leg2->SetBorderSize(1);
        leg2->Draw();

        TLatex *tex2 = new TLatex();
        tex2->SetNDC();
        tex2->SetTextSize(0.04);
        tex2->DrawLatex(0.20, 0.7, "Neutron #eta > 5.0");
        tex2->DrawLatex(0.20, 0.65, pT_text);
        tex2->DrawLatex(0.20, 0.6, "M_{X} > 0.5 GeV");
        tex2->DrawLatex(0.20, 0.55, "W^{2} > 4 GeV^{2}");
        tex2->SetTextColor(kRed);
        tex2->DrawLatex(0.85, 0.85, "(b)");

        c->Update();
        c->Print(Form("%s.png", name));
    };

    draw_fig7(hXLSig_pT03, hXLBkg_pT03, "Figure7_pT03", "Figure 7 (pT < 0.3)", "Neutron p_{T} < 0.3 GeV");
    draw_fig7(hXLSig_pT05, hXLBkg_pT05, "Figure7_pT05", "Figure 7 (pT < 0.5)", "Neutron p_{T} < 0.5 GeV");
    
    TCanvas *c7 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("Figure7_pT03");
    if(c7) c7->Print("Figure7.png");
}
