#include"KineCal.h"
#include"TaggedN_DIS.h"
#include"piIMParton.h"

#include"KineCal.cpp"
#include"TaggedN_DIS.cpp"
#include"piIMParton.cpp"

void run_sim_pion(int nEvents = 1000000){
	TaggedN_DIS dis;
	dis.SetQ2max(50);
	dis.SetQ2min(1);
	dis.SetTmax(1.0);
	dis.SetTmin(0.01);
	dis.SetxLmax(0.999);
	dis.SetxLmin(0.001); // Full physical range [0,1]; paper cuts (xL>0.75 etc.) are applied at analysis level

	dis.SetxBmax(1.0);
	dis.SetxBmin(0.0001);

	dis.SetOutputFileName("TaggedNeutron-DIS-EicC.root"); 

	dis.SetElecBeamEnergy(3.5);
	dis.SetProtBeamEnergy(20);
	dis.SetBeamCrossAngle(0.05);  //// 50 mrad

	dis.SetSamplingMode(0); // Save all events and use d4sigma for weighting
	dis.Generate(nEvents);
}
