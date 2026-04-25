#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>

void debug_fig3a() {
    TFile *f = TFile::Open("tagged-neutron-DIS/TaggedNeutron-DIS-EicC.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    double d4sigma, xpi;
    TLorentzVector *elec_out = 0;
    TLorentzVector *neut_out = 0;
    tree->SetBranchAddress("d4sigma", &d4sigma);
    tree->SetBranchAddress("xpi", &xpi);
    tree->SetBranchAddress("elec_out", &elec_out);
    tree->SetBranchAddress("neut_out", &neut_out);
    
    // Check a few events
    int cnt = 0;
    double min_eta = 999, max_eta = -999;
    double min_pe = 999, max_pe = 0;
    double min_theta_n = 999, max_theta_n = 0;
    double min_pn = 999, max_pn = 0;
    
    for (int i = 0; i < tree->GetEntries() && i < 100000; i++) {
        tree->GetEntry(i);
        if (d4sigma < 0) continue;
        if (xpi >= 1.0) continue;
        
        double eta_e = elec_out->Eta();
        double pe = elec_out->P();
        double theta_n = 180.0 - neut_out->Theta() * 180.0 / TMath::Pi();
        double pn = neut_out->P();
        
        if (eta_e < min_eta) min_eta = eta_e;
        if (eta_e > max_eta) max_eta = eta_e;
        if (pe < min_pe) min_pe = pe;
        if (pe > max_pe) max_pe = pe;
        if (theta_n < min_theta_n) min_theta_n = theta_n;
        if (theta_n > max_theta_n) max_theta_n = theta_n;
        if (pn < min_pn) min_pn = pn;
        if (pn > max_pn) max_pn = pn;
        
        if (cnt < 5) {
            cout << "Event " << i << ": eta_e=" << eta_e 
                 << " pe=" << pe << " d4sigma=" << d4sigma 
                 << " theta_n=" << theta_n << " pn=" << pn << endl;
        }
        cnt++;
    }
    cout << "\nTotal valid events (first 100k): " << cnt << endl;
    cout << "eta_e range: [" << min_eta << ", " << max_eta << "]" << endl;
    cout << "p_e range: [" << min_pe << ", " << max_pe << "]" << endl;
    cout << "theta_n range: [" << min_theta_n << ", " << max_theta_n << "]" << endl;
    cout << "p_n range: [" << min_pn << ", " << max_pn << "]" << endl;
}
