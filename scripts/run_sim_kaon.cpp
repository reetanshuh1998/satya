#include"KineCal.h"
#include"TaggedLambda_DIS.h"
#include"piIMParton.h"

#include"KineCal.cpp"
#include"TaggedLambda_DIS.cpp"
#include"piIMParton.cpp"

void run_sim_kaon(int nEvents = 1000000){
	TaggedLambda_DIS dis;
	dis.SetQ2max(50);
	dis.SetQ2min(1);
	dis.SetTmax(1.0);
	dis.SetTmin(0.01);
	dis.SetxLmax(0.999);
	dis.SetxLmin(0.6);

	dis.SetOutputFileName("TaggedLambda-DIS-EicC.root"); 

	dis.SetElecBeamEnergy(3.5);
	dis.SetProtBeamEnergy(20);
	dis.SetBeamCrossAngle(0.05);  //// 50 mrad

	dis.SetSamplingMode(0); 
	dis.Generate(nEvents);
}
