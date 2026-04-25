#include"KineCal.h"
#include"TaggedN_DIS.h"
#include"piIMParton.h"

#include"KineCal.cpp"
#include"TaggedN_DIS.cpp"
#include"piIMParton.cpp"

void run_sim_pion_scan(int nEvents = 500000){
	TaggedN_DIS dis;
	dis.SetQ2max(50);
	dis.SetQ2min(1);
	dis.SetTmax(1.0);
	dis.SetTmin(0.01);
	dis.SetxLmax(0.999);
	dis.SetxLmin(0.10); // Ultra-low to see where physics cuts off

	dis.SetxBmax(1.0);
	dis.SetxBmin(0.0001);

	dis.SetOutputFileName("scan_xL.root"); 

	dis.SetElecBeamEnergy(3.5);
	dis.SetProtBeamEnergy(20);
	dis.SetBeamCrossAngle(0.05);

	dis.SetSamplingMode(0); 
	dis.Generate(nEvents);
}
