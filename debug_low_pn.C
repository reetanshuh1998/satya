#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>

void debug_low_pn() {
    TFile *f = TFile::Open("tagged-neutron-DIS/TaggedNeutron-DIS-EicC.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    double xL, d4sigma, W2, MX2, t, Q2, xB, xpi;
    TLorentzVector *neut_out = 0;
    tree->SetBranchAddress("xL", &xL);
    tree->SetBranchAddress("d4sigma", &d4sigma);
    tree->SetBranchAddress("W2", &W2);
    tree->SetBranchAddress("MX2", &MX2);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("xB", &xB);
    tree->SetBranchAddress("xpi", &xpi);
    tree->SetBranchAddress("neut_out", &neut_out);
    
    int low_pn_count = 0;
    int printed = 0;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (d4sigma < 0) continue;
        if (xL <= 0.75) continue;
        
        double theta_lab = neut_out->Theta() * 180.0 / TMath::Pi();
        double theta_plot = 180.0 - theta_lab;
        double pn = neut_out->P();
        
        if (pn < 12.0) {
            low_pn_count++;
            if (printed < 10) {
                cout << "LOW pn=" << pn << " xL=" << xL << " |t|=" << fabs(t) 
                     << " Q2=" << Q2 << " xB=" << xB << " xpi=" << xpi
                     << " W2=" << W2 << " MX2=" << MX2 
                     << " d4sigma=" << d4sigma 
                     << " theta_n=" << theta_plot << endl;
                printed++;
            }
        }
    }
    cout << "\nTotal events with pn < 12 and xL > 0.75: " << low_pn_count << endl;
    cout << "Total events with xL > 0.75 and d4sigma > 0: ";
    int total = 0;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (d4sigma > 0 && xL > 0.75) total++;
    }
    cout << total << endl;
    
    // Check W2 distribution for low pn events
    int low_pn_lowW2 = 0;
    int low_pn_lowMX2 = 0;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (d4sigma < 0 || xL <= 0.75) continue;
        double pn = neut_out->P();
        if (pn < 12.0) {
            if (W2 < 4.0) low_pn_lowW2++;
            if (MX2 < 0.25) low_pn_lowMX2++;
        }
    }
    cout << "Of those, W2 < 4: " << low_pn_lowW2 << endl;
    cout << "Of those, MX2 < 0.25: " << low_pn_lowMX2 << endl;
}
