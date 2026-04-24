#include"KineCal.h"
#include"TaggedN_DIS.h"
#include"piIMParton.h"

#include"KineCal.cpp"
#include"TaggedN_DIS.cpp"
#include"piIMParton.cpp"


int test(){

	TaggedN_DIS dis;

	dis.SetQ2max(50);
	dis.SetQ2min(1);
	dis.SetTmax(1);
	dis.SetTmin(0.01);
	dis.SetxLmax(0.995);
	dis.SetxLmin(0.5);

	//char filename[50] = "TaggedNeutron-DIS-EicC.root";
	TString filename = "TaggedNeutron-DIS-EicC.root";
	dis.SetOutputFileName("TaggedNeutron-DIS-EicC.root"); 


	dis.SetElecBeamEnergy(3.5);
	dis.SetProtBeamEnergy(20);
	dis.SetBeamCrossAngle(0.05);  //// 50 mrad
	dis.SetBeamCrossAngle(0);  //// 0 mrad

	//dis.SetQuiet(1);
	dis.SetSamplingMode(1);
	dis.Generate(100000);

	return 0;
}


