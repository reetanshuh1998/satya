#include "../tagged-neutron-DIS/KineCal.h"
#include "../tagged-neutron-DIS/TaggedN_DIS.h"
#include "../tagged-neutron-DIS/piIMParton.h"
#include "../tagged-neutron-DIS/KineCal.cpp"
#include "../tagged-neutron-DIS/TaggedN_DIS.cpp"
#include "../tagged-neutron-DIS/piIMParton.cpp"
#include "TRandom3.h"
#include <iostream>

void estimate_sigma(){
    TaggedN_DIS dis;
    dis.SetQ2max(50);
    dis.SetQ2min(1);
    dis.SetTmax(1.0);
    dis.SetTmin(0.01);
    dis.SetxLmax(0.999);
    dis.SetxLmin(0.37);
    dis.SetxBmax(1.0);
    dis.SetxBmin(0.0001);
    dis.SetElecBeamEnergy(3.5);
    dis.SetProtBeamEnergy(20);
    dis.SetBeamCrossAngle(0.05);
    dis.SetSamplingMode(0); 
    dis.SetQuiet(1);

    TRandom3 rnd(0);
    int n = 100000;
    double sum_sigma = 0;
    
    double s = (3.5 + 20.0) * (3.5 + 20.0) - (20.0 - 3.5) * (20.0 - 3.5); // Very rough s
    // Actually use the class s
    // dis.s is private... damn.
    
    for(int i=0; i<n; i++){
        double xB = rnd.Uniform(0.0001, 1.0);
        double Q2 = rnd.Uniform(1.0, 50.0);
        double xL = rnd.Uniform(0.37, 0.999);
        double t = rnd.Uniform(-1.0, -0.01);
        
        double d4s = dis.d4sigma_dQ2dxBdxLdt_GRV(Q2, xB, xL, t);
        if(d4s > 0) {
            sum_sigma += d4s;
        }
    }
    
    double avg_sigma = sum_sigma / n;
    double V = (50.0-1.0) * (1.0-0.0001) * (0.999-0.37) * (1.0-0.01);
    double total_sigma = avg_sigma * V;
    
    std::cout << "Average d4sigma: " << avg_sigma << " nb/GeV^4" << std::endl;
    std::cout << "Phase space volume V: " << V << " GeV^4" << std::endl;
    std::cout << "Total Cross Section sigma: " << total_sigma << " nb" << std::endl;
}
