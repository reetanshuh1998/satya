#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace Pythia8;

int main() {
    // Generator.
    Pythia pythia;

    // Beams: asymmetric 3.5 GeV e- on 20 GeV p.
    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eA = 3.5");
    pythia.readString("Beams:eB = 20.0");

    // Process selection.
    pythia.readString("WeakBosonExchange:all = on");
    pythia.readString("PhaseSpace:Q2Min = 1.0");

    // Initialize.
    pythia.init();

    // ROOT setup.
    TFile *f = new TFile("DIS_Background_EicC.root", "RECREATE");
    TTree *tree = new TTree("tree", "Pythia8 DIS Background");
    double xL, eta, pT, Q2, xB, W2, MX2;
    int id;
    tree->Branch("xL", &xL, "xL/D");
    tree->Branch("eta", &eta, "eta/D");
    tree->Branch("pT", &pT, "pT/D");
    tree->Branch("Q2", &Q2, "Q2/D");
    tree->Branch("xB", &xB, "xB/D");
    tree->Branch("W2", &W2, "W2/D");
    tree->Branch("MX2", &MX2, "MX2/D");
    tree->Branch("id", &id, "id/I");

    // Number of events.
    int nEvents = 100; 

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        // Manual DIS kinematic calculation
        Vec4 pLeptonIn = pythia.event[1].p();
        Vec4 pHadronIn = pythia.event[2].p();
        Vec4 pLeptonOut;
        bool foundLepton = false;
        
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal() && pythia.event[i].id() == 11) {
                pLeptonOut = pythia.event[i].p();
                foundLepton = true;
                break;
            }
        }
        
        if (!foundLepton) continue;

        Vec4 q = pLeptonIn - pLeptonOut;
        Q2 = -q.m2Calc();
        double pDotQ = pHadronIn * q;
        xB = Q2 / (2.0 * pDotQ);
        W2 = (pHadronIn + q).m2Calc();

        // Only keep DIS events in the kinematic region
        if (W2 < 4.0) continue;

        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal() && (abs(pythia.event[i].id()) == 3122 || abs(pythia.event[i].id()) == 2112)) {
                Vec4 pNeutron = pythia.event[i].p();
                id = pythia.event[i].id();
                eta = pythia.event[i].eta();
                pT = pythia.event[i].pT();
                xL = pNeutron.e() / 20.0; 
                MX2 = (pHadronIn + q - pNeutron).m2Calc();
                
                tree->Fill();
            }
        }
        if (iEvent % 10000 == 0) std::cout << iEvent << " events processed" << std::endl;
    }

    pythia.stat();
    tree->Write();
    f->Close();
    
    // Print the generated cross section
    double sigma_gen_mb = pythia.info.sigmaGen();
    std::cout << "Pythia 8 DIS cross-section (mb): " << sigma_gen_mb << std::endl;
    std::cout << "Pythia 8 DIS cross-section (nb): " << sigma_gen_mb * 1e6 << std::endl;

    return 0;
}
