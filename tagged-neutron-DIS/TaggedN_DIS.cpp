#include"TaggedN_DIS.h"

#include<iostream>
#include<string.h>
#include<cmath>
#include"TMath.h"

using namespace std;


TaggedN_DIS::TaggedN_DIS(){
	cout<<"****EicC Meson Structure Project"<<endl;
	cout<<"****Coding issues, contact rwang@impcas.ac.cn"<<endl;
	cout<<endl<<endl;
	cout<<"    Simulation starting..."<<endl;
	cout<<"    Process: e p --> e n X"<<endl;

	///// the kinematical ranges for MC sampling
	xBmin = 0.0001;
	xBmax = 0.95;
	Q2min = 1;
	Q2max = 100;
	xLmin = 0.1;
	xLmax = 0.95;
	t0 = -0.001;
	t1 = -100;
	ymin = 0.01;
	ymax = 0.99;
	Tmin = 0.01;
	Tmax = 50;

	sampling_flag = 0;
	quiet_flag = 0;
	max_d4sigma = 2200;

	//// nucleon mass and electron mass
	mN = 0.938272;
	mpi = 0.14;
	me = 0.000511;

	PI = 3.141592653;

	strFileName = new char[100];
	strcpy(strFileName, "DIS_tagN.root");

	//// EicC optimal collision energy
	eBeamE = 3.5;
	pBeamE = 20;
	beam_cross_angle = 0.05; //// 50 mrad
	eBeam = new TLorentzVector(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam = new TLorentzVector(p_px, 0, p_pz, pBeamE);
	BoostToEIC = new TVector3(0,0,0);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();

	elec_out = new TLorentzVector(0, 0, 0, 0);
	neut_out = new TLorentzVector(0, 0, 0, 0);

	/// random seed
	random = new TRandom3(0);
	/// kinematic calculator
	kine = new KineCal();
	///// pion PDFs by IMParton Collaboration.
	pionPDFIMP = new piIMParton();

}
TaggedN_DIS::~TaggedN_DIS(){
	//tree->Write();
	//fout->Write();
	//fout->Close();
	//cout<<"    Data file saved and closed~"<<endl<<endl;
	//
	//delete random;
	//delete tree;
	//delete fout;
	//delete kine;
	//delete eBeam;
	//delete pBeam;
	//delete elec_out;
	//delete neut_out;
}

int TaggedN_DIS::Generate(int N = 20000){

	MakeROOTFile(strFileName);

	cout<<"    To generate "<<N<<" events..."<<endl;
	TVector3 zdirection = eBeam->Vect();
	TVector3 ydirection(0, 1, 0);
	TVector3 xdirection(zdirection.Z(), 0, -zdirection.X());

	for(int i=0; i<N; ){
		xB = random->Uniform(xBmin, xBmax);
		Q2 = random->Uniform(Q2min, Q2max);
		if(Q2>(xB*ymax*(s-mN*mN)))continue;
		xL = random->Uniform(xLmin, xLmax);
		W2 = (1.0/xB-1)*Q2 + mN*mN;
		//if(W2<2.6)continue;
		if(W2<3.5)continue;
		MX2 = (1-xL)*(W2+Q2) - Q2;
		if(MX2<mpi*mpi)continue;
		//if(MX2<3.5)continue;
		t0 = kine->calTMin(W2, -Q2, mN*mN, MX2, mN*mN);
		if(!(t0==t0))continue;
		if(t0<(-Tmax))continue;
		if(t0>(-Tmin))t0 = -Tmin;
		t1 = kine->calTMax(W2, -Q2, mN*mN, MX2, mN*mN);
		if(!(t1==t1))continue;
		if(t1>(-Tmin))continue;
		if(t1<(-Tmax))t1 = -Tmax;
		//cout<<"t0 "<<t0<<endl<<t1<<endl;  //test code.
		t = random->Uniform(t1, t0);		
		xpi = xB/(1-xL);
		y = Q2 / xB / (s-mN*mN);

		double nv = Q2/2.0/mN/xB;
		double eOutE = eBeamENRest - nv;
		//cout<<eBeamENRest<<"\t"<<nv<<endl;   // test code.
		double cosetheta = 1 - Q2/2.0/eBeamENRest/eOutE;
		double sinetheta = sqrt(1 - cosetheta*cosetheta);
		double ephi = random->Uniform(0, 2*PI);
		double emom = sqrt(eOutE*eOutE - me*me);
		TVector3 emom_v3 = emom*sinetheta*cos(ephi) * xdirection.Unit();
		emom_v3 += emom*sinetheta*sin(ephi) * ydirection.Unit();
		emom_v3 += emom*cosetheta * zdirection.Unit();
		elec_out->SetXYZT(emom_v3.X(), emom_v3.Y(), emom_v3.Z(), eOutE);

		double neutE = (2*mN*mN-t) / 2.0 / mN;
		double nphi = random->Uniform(0, 2*PI);
		double neutmom = sqrt(neutE*neutE-mN*mN);
		double qx = - elec_out->Px();
		double qy = - elec_out->Py();
		double qz = sqrt(eBeamENRest*eBeamENRest-me*me) - elec_out->Pz();
		double k1 = qx*neutmom*cos(nphi) + qy*neutmom*sin(nphi);
		double k2 = qz*neutmom;
		double k3 = nv*neutE - xL*mN*nv;
		double A = k1*k1 + k2*k2;
		double B = -2.0*k2*k3;
		double C = k3*k3 - k1*k1;
		if(B*B < 4*A*C)continue;
		double cosntheta = (-B + sqrt(B*B-4*A*C)) / 2.0 / A;
		if(cosntheta<-1)continue;
		if(cosntheta>1)cosntheta = (-B - sqrt(B*B-4*A*C)) / 2.0 / A;
		if(cosntheta<-1 || cosntheta>1)continue;
		double sinntheta = sqrt(1-cosntheta*cosntheta);
		neut_out->SetXYZT(neutmom*sinntheta*cos(nphi), neutmom*sinntheta*sin(nphi), neutmom*cosntheta, neutE);

		d4sigma = d4sigma_dQ2dxBdxLdt_GRV(Q2, xB, xL, t);

		//// Boost from proton beam rest frame to collider frame!
		elec_out->Boost(*BoostToEIC);
		neut_out->Boost(*BoostToEIC);

		if(sampling_flag)
			if(d4sigma < random->Uniform(0,max_d4sigma))continue;
		tree->Fill();
		i++;
		if(!quiet_flag) if(i%1000==0) cout<<i<<" events"<<endl;
	}


	cout<<"    Event generation done! "<<endl;


	tree->Write();
	//fout->Write();
	fout->Close();
	cout<<"    Data file saved and closed~"<<endl<<endl;

	return N;
}

//// Create a ROOT file and a TTree.
void TaggedN_DIS::MakeROOTFile(char *filename){
	//// create the output file and the output TTree
	cout<<"    Creating the output file: "<<filename<<endl;
	fout = new TFile(filename,"recreate");
	tree = new TTree("tree","TaggedDIS");
	tree->Branch("xB", &xB, "xB/D");
	tree->Branch("Q2", &Q2, "Q2/D");
	tree->Branch("xL", &xL, "xL/D");
	tree->Branch("t", &t, "t/D");
	tree->Branch("xpi", &xpi, "xpi/D");
	tree->Branch("y", &y, "y/D");
	tree->Branch("W2", &W2, "W2/D");
	tree->Branch("MX2", &MX2, "MX2/D");
	tree->Branch("s", &s, "s/D");
	tree->Branch("d4sigma", &d4sigma, "d4sigma/D");
	tree->Branch("elec_out", "TLorentzVector", elec_out);
	tree->Branch("neut_out", "TLorentzVector", neut_out);
}
void TaggedN_DIS::SetOutputFileName(char *filename){
	strcpy(strFileName, filename);
}
void TaggedN_DIS::SetOutputFileName(TString filename){
	strcpy(strFileName, filename.Data());
}

void TaggedN_DIS::SetElecBeamEnergy(double ebeamenergy){
	if(ebeamenergy<0.001){cout<<"Error: electron beam energy is too small!!!"<<endl; return;}
	if(ebeamenergy>1e6){cout<<"Error: electron beam energy is too high!!!"<<endl; return;}
	eBeamE = ebeamenergy;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
}

void TaggedN_DIS::SetProtBeamEnergy(double pbeamenergy){
	if(pbeamenergy<1){cout<<"Error: proton beam energy is too small!!!"<<endl; return;}
	if(pbeamenergy>1e6){cout<<"Error: proton beam energy is too high!!!"<<endl; return;}
	pBeamE = pbeamenergy;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
}
//// set beam crossing angle
void TaggedN_DIS::SetBeamCrossAngle(double _angle){
	beam_cross_angle = _angle;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
}
double TaggedN_DIS::GetBeamCrossAngle(){return beam_cross_angle;}



//// GRV's pionic PDF model, published in the year of 1992.
double TaggedN_DIS::xuv_pi_GRV(double xpi, double Q2){
	double N, a, A, D, s;
	if(xpi>0.999)return 0.0;
	s = log( log(Q2 / (0.232 * 0.232)) / log(0.25 / (0.232 * 0.232)) );
	N = 0.519 + 0.180 * s - 0.011 * s * s;
	a = 0.499 - 0.027 * s;
	A = 0.381 - 0.419 * s;
	D = 0.367 + 0.563 * s;
	return N * pow(xpi,a) * (1 + A * sqrt(xpi)) * pow((1 - xpi),D);
}
double TaggedN_DIS::xdv_pi_GRV(double xpi, double Q2){
	double N, a, A, D, s;
	if(xpi>0.999)return 0.0;
	s = log( log(Q2 / (0.232 * 0.232)) / log(0.25 / (0.232 * 0.232)) );
	N = 0.519 + 0.180 * s - 0.011 * s * s;
	a = 0.499 - 0.027 * s;
	A = 0.381 - 0.419 * s;
	D = 0.367 + 0.563 * s;
	return N * pow(xpi,a) * (1 + A * sqrt(xpi)) * pow((1 - xpi),D);
}
double TaggedN_DIS::xsbar_pi_GRV(double xpi, double Q2){
	double s, a, A, B, D, E, F;
	if(xpi>0.999)return 0.0;
	s = log( log(Q2 / (0.232 * 0.232)) / log(0.25 / (0.232 * 0.232)) );
	a = 2.538 -0.763 * s;
	A = -0.748;
	B = 0.313 + 0.935 * s;
	D = 3.359;
	E = 4.433 + 1.301 * s;
	F = 9.3 + 0.887 * s;	
	return pow(s,0.55) / pow(log(1 / xpi),a) * (1 + A * sqrt(xpi) + B * xpi) * pow((1 - xpi),D) * exp( -E + sqrt(F * pow(s,0.56) * log(1 / xpi)) );
}
double TaggedN_DIS::F2_pi_GRV(double xpi, double Q2){
	double y;
	y = 4.0/9.0*(xuv_pi_GRV(xpi,Q2) + 2*xsbar_pi_GRV(xpi,Q2)) + 1.0/9.0*(xdv_pi_GRV(xpi,Q2) + 2*xsbar_pi_GRV(xpi,Q2) + 2*xsbar_pi_GRV(xpi,Q2));
	return y;
}
double TaggedN_DIS::F2_pi_IMParton(double xpi, double Q2){
	double y;
	y =  4.0/9.0*(pionPDFIMP->getPDF(1,xpi,Q2) + pionPDFIMP->getPDF(-1,xpi,Q2))
	   + 1.0/9.0*(pionPDFIMP->getPDF(2,xpi,Q2) + pionPDFIMP->getPDF(-2,xpi,Q2) + pionPDFIMP->getPDF(3,xpi,Q2) + pionPDFIMP->getPDF(-3,xpi,Q2));
	return xpi * y;
}
//// return cross-section in the unit of nb/GeV^4.
double TaggedN_DIS::d4sigma_dQ2dxBdxLdt_GRV(double Q2, double xB, double xL, double t)
{
	double _xpi = xB / (1-xL);
	double _F2pi = F2_pi_IMParton(_xpi, Q2);
	double _y = Q2 / xB / (s-mN*mN);
        double alpha2 = pow(1/137.0, 2);
      	double f_pi_in_p = 1/2.0/PI * 13.6 * (1-xL) * (-t)/pow(mpi*mpi-t, 2) * exp(-0.93*0.93*(mpi*mpi-t)/(1-xL));
	//// return cross-section in nb.
	double sigma = 3.881e5 * 4*PI*alpha2 /xB /pow(Q2,2) * (1.0-_y+_y*_y/2.0) * _F2pi * f_pi_in_p;
	if(xpi>0.999)return -30;
	//// return 0.0 if the value is NaN.
	if(!(sigma==sigma))return -20;
	if(sigma<0)return -10;
	else return sigma;
}


//// set sampling ranges
void TaggedN_DIS::SetxBmin(double min){xBmin = min;}
void TaggedN_DIS::SetxBmax(double max){xBmax = max;}
void TaggedN_DIS::SetQ2min(double min){Q2min = min;}
void TaggedN_DIS::SetQ2max(double max){Q2max = max;}
void TaggedN_DIS::SetxLmin(double min){xLmin = min;}
void TaggedN_DIS::SetxLmax(double max){xLmax = max;}
void TaggedN_DIS::SetTmin(double min){Tmin = min;}
void TaggedN_DIS::SetTmax(double max){Tmax = max;}
void TaggedN_DIS::Setymin(double min){ymin = min;}
void TaggedN_DIS::Setymax(double max){ymax = max;}

int TaggedN_DIS::SetSamplingMode(int flag){
	sampling_flag = flag;
	return sampling_flag;
}
int TaggedN_DIS::SetQuiet(int flag){quiet_flag = flag; return quiet_flag;}
