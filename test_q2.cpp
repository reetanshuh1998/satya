#include "Pythia8/Pythia.h"
#include <iostream>

using namespace Pythia8;

int main() {
    Pythia pythia;
    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eA = 3.5");
    pythia.readString("Beams:eB = 20.0");
    pythia.readString("WeakBosonExchange:all = on");
    pythia.readString("PhaseSpace:Q2Min = 0.01"); // Test Q2 > 0.1
    pythia.init();
    
    // Just run 10 events to get the cross section
    for (int iEvent = 0; iEvent < 10; ++iEvent) pythia.next();
    
    pythia.stat();
    std::cout << "SigmaGen: " << pythia.info.sigmaGen() * 1e6 << " nb" << std::endl;
    return 0;
}
