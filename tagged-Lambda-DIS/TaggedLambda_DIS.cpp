#include"TaggedLambda_DIS.h"

#include<iostream>
#include<string.h>
#include<cmath>
#include"TMath.h"

using namespace std;


TaggedLambda_DIS::TaggedLambda_DIS(){
	cout<<"****EicC Meson Structure Project"<<endl;
	cout<<"****Coding issues, contact rwang@impcas.ac.cn"<<endl;
	cout<<endl<<endl;
	cout<<"    Simulation starting..."<<endl;
	cout<<"    DIS process: e p --> e Lambda X"<<endl;
	cout<<"    Decay process: Lambda --> pi^zero n"<<endl;

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
	max_d4sigma = 100;


	//// nucleon mass and electron mass
	mN = 0.938272;
	mLambda = 1.115;
	mpi =0.139;
	mkaon = 0.49368;
	me = 0.000511;
	pi0mass = 0.135;
	
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
	lambda_out = new TLorentzVector(0, 0, 0, 0);
	decay_neut = new TLorentzVector(0, 0, 0, 0);
	decay_neut_smeared = new TLorentzVector(0, 0, 0, 0);
	decay_pi0 = new TLorentzVector(0, 0, 0, 0);

	decay_gamma_1 = new TLorentzVector(0, 0, 0, 0);
	decay_gamma_2 = new TLorentzVector(0, 0, 0, 0);
	decay_gamma_1_smeared = new TLorentzVector(0, 0, 0, 0);
	decay_gamma_2_smeared = new TLorentzVector(0, 0, 0, 0);

	/// random seed
	random = new TRandom3(0);
	/// kinematic calculator
	kine = new KineCal();
	///// pion PDFs by IMParton Collaboration.
	Kaonpdf = new piIMParton();
	Kaonpdf->setDataSet(321);

}
TaggedLambda_DIS::~TaggedLambda_DIS(){
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
	//delete lambda_out;
}

int TaggedLambda_DIS::Generate(int N = 20000000){

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
		if(MX2<mLambda*mLambda)continue;
		//if(MX2<3.5)continue;
		t0 = kine->calTMin(W2, -Q2, mN*mN, MX2, mLambda*mLambda);
		if(!(t0==t0))continue;
		if(t0<(-Tmax))continue;
		if(t0>(-Tmin))t0 = -Tmin;
		t1 = kine->calTMax(W2, -Q2, mN*mN, MX2, mLambda*mLambda);
		if(!(t1==t1))continue;
		if(t1>(-Tmin))continue;
		if(t1<(-Tmax))t1 = -Tmax;
		//cout<<"t0 "<<t0<<endl<<t1<<endl;  //test code.
		t = random->Uniform(t1, t0);		
		xkaon = xB/(1-xL);
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

		double LambdaE = (mN*mN+mLambda*mLambda-t) / 2.0 / mN;
		double Lambdaphi = random->Uniform(0, 2*PI);
		double Lambdamom = sqrt(LambdaE*LambdaE-mLambda*mLambda);
		double qx = - elec_out->Px();
		double qy = - elec_out->Py();
		double qz = sqrt(eBeamENRest*eBeamENRest-me*me) - elec_out->Pz();
		double k1 = qx*Lambdamom*cos(Lambdaphi) + qy*Lambdamom*sin(Lambdaphi);
		double k2 = qz*Lambdamom;
		double k3 = nv*LambdaE - xL*mN*nv;
		double A = k1*k1 + k2*k2;
		double B = -2.0*k2*k3;
		double C = k3*k3 - k1*k1;
		if(B*B < 4*A*C)continue;
		double cosLambdatheta = (-B + sqrt(B*B-4*A*C)) / 2.0 / A;
		if(cosLambdatheta<-1)continue;
		if(cosLambdatheta>1)cosLambdatheta = (-B - sqrt(B*B-4*A*C)) / 2.0 / A;
		if(cosLambdatheta<-1 || cosLambdatheta>1)continue;
		double sinLambdatheta = sqrt(1-cosLambdatheta*cosLambdatheta);
		lambda_out->SetXYZT(Lambdamom*sinLambdatheta*cos(Lambdaphi), Lambdamom*sinLambdatheta*sin(Lambdaphi), Lambdamom*cosLambdatheta, LambdaE);

		d4sigma = d4sigma_dQ2dxBdxLdt_kaon(Q2, xB, xL, t);
		if(sampling_flag)
			if(d4sigma < random->Uniform(0,max_d4sigma))continue;


/*************************************************************************************/
		// Lambda decay
		double costheta_decay = random->Uniform(-1,1);
 		double phi_decay  =  random->Uniform(-3.141592653, 3.141592653);
		double sintheta_decay = sqrt(1 - costheta_decay*costheta_decay);
		decay_neut->SetXYZT(0.104*sintheta_decay*cos(phi_decay), 0.104*sintheta_decay*sin(phi_decay), 0.104*costheta_decay, sqrt(0.104*0.104 + mN*mN)  );
		decay_pi0->SetXYZT(-0.104*sintheta_decay*cos(phi_decay), -0.104*sintheta_decay*sin(phi_decay), -0.104*costheta_decay, sqrt(0.104*0.104 + pi0mass*pi0mass) );

/*************************************************************************************/
		// pi0 decay    
		double costheta_decay_pi0 = random->Uniform(-1,1);
 		double phi_decay_pi0  =  random->Uniform(-3.141592653, 3.141592653);
		double sintheta_decay_pi0 = sqrt(1 - costheta_decay_pi0*costheta_decay_pi0);
		decay_gamma_1->SetXYZT( 0.067*sintheta_decay_pi0*cos(phi_decay_pi0),  0.067*sintheta_decay_pi0*sin(phi_decay_pi0),  0.067*costheta_decay_pi0, 0.067 );
		decay_gamma_2->SetXYZT(-0.067*sintheta_decay_pi0*cos(phi_decay_pi0), -0.067*sintheta_decay_pi0*sin(phi_decay_pi0), -0.067*costheta_decay_pi0, 0.067 );




		//// Boost from Lambda rest frame to proton (target) rest frame
		decay_neut->Boost( lambda_out->BoostVector() );
		decay_pi0->Boost( lambda_out->BoostVector() );

		//// Boost from pi0 rest frame to proton (target) rest frame
		decay_gamma_1->Boost( decay_pi0->BoostVector() );
		decay_gamma_2->Boost( decay_pi0->BoostVector() );

		//// Boost from proton beam rest frame to collider frame!

		elec_out->Boost(*BoostToEIC);
		lambda_out->Boost(*BoostToEIC);
		decay_neut->Boost(*BoostToEIC);
		decay_pi0->Boost(*BoostToEIC);
		decay_gamma_1->Boost(*BoostToEIC);
		decay_gamma_2->Boost(*BoostToEIC);

		if(lambda_out->E()>4){
			double e_ = lambda_out->E();
			double theta_ = lambda_out->Theta();
			double phi_ = lambda_out->Phi();
			double delta_ = e_ * random->Gaus(0, 0.35/sqrt(e_));
			e_ = e_ + delta_;
			double p_ = sqrt(e_*e_ - mLambda*mLambda);
			delta_ = theta_ * random->Gaus(0, 0.001);
			theta_ = theta_ + delta_;
			if(theta_>=3.14159265) theta_ = 3.14159265;
			TLorentzVector lam_(p_*sin(theta_)*cos(phi_), p_*sin(theta_)*sin(phi_), p_*cos(theta_), e_);
			TLorentzVector qq_ = *eBeam - *elec_out;
			xL_smeared = (lam_*qq_) / ((*pBeam)*qq_);

			e_ = decay_neut->E();
			theta_ = decay_neut->Theta();
			phi_ = decay_neut->Phi();
			delta_ = e_ * random->Gaus(0, 0.35/sqrt(e_));
			e_ = e_ + delta_;
			p_ = sqrt(e_*e_ - mN*mN);
			delta_ = theta_ * random->Gaus(0, 0.001);
			theta_ = theta_ + delta_;
			if(theta_>=3.14159265) theta_ = 3.14159265;
			decay_neut_smeared->SetXYZT(p_*sin(theta_)*cos(phi_), p_*sin(theta_)*sin(phi_), p_*cos(theta_), e_);
		
			if(decay_gamma_1->E()>0.1){
				e_ = decay_gamma_1->E();
				theta_ = decay_gamma_1->Theta();
				phi_ = decay_gamma_1->Phi();
				delta_ = e_ * random->Gaus(0, 0.03/sqrt(e_));
				e_ = e_ + delta_;
				p_ = e_;
				delta_ = theta_ * random->Gaus(0, 0.0004);
				theta_ = theta_ + delta_;
				if(theta_>=3.14159265) theta_ = 3.14159265;
				decay_gamma_1_smeared->SetXYZT(p_*sin(theta_)*cos(phi_), p_*sin(theta_)*sin(phi_), p_*cos(theta_), e_);
			}
			else decay_gamma_1_smeared->SetXYZT(0, 0, 0, 0);

			if(decay_gamma_2->E()>0.1){
				e_ = decay_gamma_2->E();
				theta_ = decay_gamma_2->Theta();
				phi_ = decay_gamma_2->Phi();
				delta_ = e_ * random->Gaus(0, 0.03/sqrt(e_));
				e_ = e_ + delta_;
				p_ = e_;
				delta_ = theta_ * random->Gaus(0, 0.0004);
				theta_ = theta_ + delta_;
				if(theta_>=3.14159265) theta_ = 3.14159265;
				decay_gamma_2_smeared->SetXYZT(p_*sin(theta_)*cos(phi_), p_*sin(theta_)*sin(phi_), p_*cos(theta_), e_);
			}
			else decay_gamma_2_smeared->SetXYZT(0, 0, 0, 0);

			mpi0_rec = ((*decay_gamma_1_smeared) + (*decay_gamma_2_smeared)).M();
			mLambda_rec = ((*decay_gamma_1_smeared) + (*decay_gamma_2_smeared) + (*decay_neut_smeared)).M();

		}
		else xL_smeared = 0;



/*************************************************************************************/

		//// calculate the decay vertex of Lambda and the momentum of the decays		
		double ctau = 0.07896;
 		double length = lambda_out->Beta() * lambda_out->Gamma() * ctau;
		double lengthMC = random->Exp(length);

 		gen_lamb_SecVx = lengthMC * sin(lambda_out->Theta()  ) * cos(lambda_out->Phi()  );
 		gen_lamb_SecVy = lengthMC * sin(lambda_out->Theta()  ) * sin(lambda_out->Phi()  );
 		gen_lamb_SecVz = lengthMC * cos(lambda_out->Theta()  ) ;

		//// calculate the second decay vertex of pi0 and the momentum of the decays
		double ctau_pi0 = 2.52e-8;
 		double length_pi0 = decay_pi0->Beta() * decay_pi0->Gamma() * ctau_pi0;
		double lengthMC_pi0 = random->Exp(length_pi0);

 		gen_pi0_ThirdVx = lengthMC_pi0 * sin(decay_pi0->Theta()  ) * cos(decay_pi0->Phi()  );
 		gen_pi0_ThirdVy = lengthMC_pi0 * sin(decay_pi0->Theta()  ) * sin(decay_pi0->Phi()  );
 		gen_pi0_ThirdVz = lengthMC_pi0 * cos(decay_pi0->Theta()  ) ;
		
				
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
void TaggedLambda_DIS::MakeROOTFile(char *filename){
	//// create the output file and the output TTree
	cout<<"    Creating the output file: "<<filename<<endl;
	fout = new TFile(filename,"recreate");
	tree = new TTree("tree","TaggedDIS");
	tree->Branch("xB", &xB, "xB/D");
	tree->Branch("Q2", &Q2, "Q2/D");
	tree->Branch("xL", &xL, "xL/D");
	tree->Branch("t", &t, "t/D");
	tree->Branch("xkaon", &xkaon, "xkaon/D");
	tree->Branch("y", &y, "y/D");
	tree->Branch("W2", &W2, "W2/D");
	tree->Branch("MX2", &MX2, "MX2/D");
	tree->Branch("s", &s, "s/D");
	tree->Branch("gen_lamb_SecVx", &gen_lamb_SecVx, "gen_lamb_SecVx/D");
	tree->Branch("gen_lamb_SecVy", &gen_lamb_SecVy, "gen_lamb_SecVy/D");
	tree->Branch("gen_lamb_SecVz", &gen_lamb_SecVz, "gen_lamb_SecVz/D");
	tree->Branch("gen_pi0_ThirdVx", &gen_pi0_ThirdVx, "gen_pi0_ThirdVx/D");
	tree->Branch("gen_pi0_ThirdVy", &gen_pi0_ThirdVy, "gen_pi0_ThirdVy/D");
	tree->Branch("gen_pi0_ThirdVz", &gen_pi0_ThirdVz, "gen_pi0_ThirdVz/D");
	tree->Branch("d4sigma", &d4sigma, "d4sigma/D");
	tree->Branch("elec_out", "TLorentzVector", elec_out);
	tree->Branch("lambda_out", "TLorentzVector", lambda_out);
	tree->Branch("decay_neut", "TLorentzVector", decay_neut);
	tree->Branch("decay_pi0", "TLorentzVector", decay_pi0);
	tree->Branch("decay_gamma_1", "TLorentzVector", decay_gamma_1);
	tree->Branch("decay_gamma_2", "TLorentzVector", decay_gamma_2);	
	tree->Branch("decay_neut_smeared", "TLorentzVector", decay_neut_smeared);
	tree->Branch("decay_gamma_1_smeared", "TLorentzVector", decay_gamma_1_smeared);
	tree->Branch("decay_gamma_2_smeared", "TLorentzVector", decay_gamma_2_smeared);	

	tree->Branch("xL_smeared", &xL_smeared, "xL_smeared/D");
	tree->Branch("mpi0_rec", &mpi0_rec, "mpi0_rec/D");
	tree->Branch("mLambda_rec", &mLambda_rec, "mLambda_rec/D");
}
void TaggedLambda_DIS::SetOutputFileName(char *filename){
	strcpy(strFileName, filename);
}
void TaggedLambda_DIS::SetOutputFileName(TString filename){
	strcpy(strFileName, filename.Data());
}

void TaggedLambda_DIS::SetElecBeamEnergy(double ebeamenergy){
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

void TaggedLambda_DIS::SetProtBeamEnergy(double pbeamenergy){
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
void TaggedLambda_DIS::SetBeamCrossAngle(double _angle){
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
double TaggedLambda_DIS::GetBeamCrossAngle(){return beam_cross_angle;}



// kaon structure function 
// there we use pion SF to replace kaon SF
double TaggedLambda_DIS::F2_kaon_IMParton(double xkaon, double Q2){
	double y;
	y =  4.0/9.0*(Kaonpdf->getPDF(2,xkaon,Q2) + Kaonpdf->getPDF(-2,xkaon,Q2))
	   + 1.0/9.0*(Kaonpdf->getPDF(1,xkaon,Q2) + Kaonpdf->getPDF(-1,xkaon,Q2) + Kaonpdf->getPDF(3,xkaon,Q2) + Kaonpdf->getPDF(-3,xkaon,Q2));
	return xkaon * y;
}

//// return cross-section in the unit of nb/GeV^4.
double TaggedLambda_DIS::d4sigma_dQ2dxBdxLdt_kaon(double Q2, double xB, double xL, double t)
{
	double _xkaon = xB / (1-xL);
	double _F2kaon = F2_kaon_IMParton(_xkaon, Q2);
	double _y = Q2 / xB / (s-mN*mN);
        double alpha2 = pow(1/137.0, 2);
      	double f_kaon_in_p = 1/2.0/PI * 14.688 * (1-xL) * (-t)/pow(mkaon*mkaon-t, 2) * exp(-1*1*(mkaon*mkaon-t)/(1-xL));
	//// return cross-section in nb.
	double sigma = 3.881e5 * 4*PI*alpha2 /xB /pow(Q2,2) * (1.0-_y+_y*_y/2.0) * _F2kaon * f_kaon_in_p;
	if(xkaon>0.999)return -30;
	//// return 0.0 if the value is NaN.
	if(!(sigma==sigma))return -20;
	if(sigma<0)return -10;
	else return sigma;
}


//// set sampling ranges
void TaggedLambda_DIS::SetxBmin(double min){xBmin = min;}
void TaggedLambda_DIS::SetxBmax(double max){xBmax = max;}
void TaggedLambda_DIS::SetQ2min(double min){Q2min = min;}
void TaggedLambda_DIS::SetQ2max(double max){Q2max = max;}
void TaggedLambda_DIS::SetxLmin(double min){xLmin = min;}
void TaggedLambda_DIS::SetxLmax(double max){xLmax = max;}
void TaggedLambda_DIS::SetTmin(double min){Tmin = min;}
void TaggedLambda_DIS::SetTmax(double max){Tmax = max;}
void TaggedLambda_DIS::Setymin(double min){ymin = min;}
void TaggedLambda_DIS::Setymax(double max){ymax = max;}

int TaggedLambda_DIS::SetSamplingMode(int flag){
	sampling_flag = flag;
	return sampling_flag;
}
int TaggedLambda_DIS::SetQuiet(int flag){quiet_flag = flag; return quiet_flag;}
