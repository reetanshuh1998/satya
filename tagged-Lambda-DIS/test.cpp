#include"KineCal.h"
#include"TaggedLambda_DIS.h"
#include"piIMParton.h"

#include"KineCal.cpp"
#include"TaggedLambda_DIS.cpp"
#include"piIMParton.cpp"


int test(){

	TaggedLambda_DIS dis;

	dis.SetQ2max(60);
	dis.SetQ2min(1);
	dis.SetTmax(2);
	dis.SetTmin(0.001);
	dis.SetxLmax(0.995);
	dis.SetxLmin(0.2);

	//char filename[50] = "TaggedLambda-DIS-EicC.root";
	//TString filename = "TaggedLambda-DIS-EicC.root";
	dis.SetOutputFileName("TaggedLambda-DIS-EicC.root"); 


	dis.SetElecBeamEnergy(3.5);
	dis.SetProtBeamEnergy(20);
	dis.SetBeamCrossAngle(0.05);  //// 50 mrad
	dis.SetBeamCrossAngle(0);   //// 0 mrad

	//dis.SetQuiet(1);
	dis.SetSamplingMode(1);
	dis.Generate(50000);

	return 0;
}


