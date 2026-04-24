#ifndef _TAGGEDLAMBDADIS_H
#define _TAGGEDLAMBDADIS_H 1

#include"TRandom3.h"
#include"TFile.h"
#include"TTree.h"
#include"TLorentzVector.h"
#include"TVector3.h"
#include"TString.h"

#include"KineCal.h"

#include"piIMParton.h"


class TaggedLambda_DIS{

	public:
		TaggedLambda_DIS();
		~TaggedLambda_DIS();

		/// the function to generate and dump N events into root file
		int Generate(int N);

		void SetOutputFileName(char *filename);
		void SetOutputFileName(TString filename);

	
		
	
		double F2_kaon_IMParton(double xkaon, double Q2);
		double d4sigma_dQ2dxBdxLdt_kaon(double Q2, double xB, double xL, double t);


		//// set sampling ranges
		void SetxBmin(double min);
		void SetxBmax(double max);
		void SetQ2min(double min);
		void SetQ2max(double max);
		void SetxLmin(double min);
		void SetxLmax(double max);
		void SetTmin(double min);
		void SetTmax(double max);
		void Setymin(double min);
		void Setymax(double max);
		int SetSamplingMode(int flag);
		int SetQuiet(int flag);

		//// set beam energies and crossing angle
		void SetElecBeamEnergy(double ebeamenergy);		
		void SetProtBeamEnergy(double pbeamenergy);		
		void SetBeamCrossAngle(double _angle);
		double GetBeamCrossAngle();


	private:
		int sampling_flag;
		int quiet_flag;
		double max_d4sigma;

		double me;
		double mkaon;
		double mN;
		double mLambda;
		double pi0mass;
		double mpi;
		double PI;

		double xB;
		double xL;
		double Q2;
		double t;
		double xkaon;
		double y;
		double W2;
		double MX2;
		double s;

		double d4sigma;

		double xBmin;
		double xBmax;
		double Q2min;
		double Q2max;
		double xLmin;
		double xLmax;
		double t0;
		double t1;
		double ymin;
		double ymax;
		double Tmin;
		double Tmax;
		
		double gen_lamb_SecVx;
		double gen_lamb_SecVy;
		double gen_lamb_SecVz;

		double gen_pi0_ThirdVx;
 		double gen_pi0_ThirdVy;
 		double gen_pi0_ThirdVz;



		double xL_smeared;
		double mpi0_rec;
		double mLambda_rec;

		TLorentzVector *elec_out;
		TLorentzVector *lambda_out;
		
		// Option A and Option B	
		// A
	/*	TLorentzVector *decay_prot;
		TLorentzVector *decay_pim;*/
		
		// B
		TLorentzVector *decay_neut;
		TLorentzVector *decay_neut_smeared;
		TLorentzVector *decay_pi0;
		
		TLorentzVector *decay_gamma_1;
		TLorentzVector *decay_gamma_2;
		TLorentzVector *decay_gamma_1_smeared;
		TLorentzVector *decay_gamma_2_smeared;

		double beam_cross_angle;
		double eBeamE;
		double pBeamE;
		double eBeamENRest;
		TLorentzVector *eBeam;
		TLorentzVector *pBeam;

		TVector3 *BoostToEIC;

		TRandom3 *random;
		TFile *fout;
		TTree *tree;
		KineCal *kine;
		piIMParton *Kaonpdf;
		char *strFileName;

		void MakeROOTFile(char *filename);
};


#endif
