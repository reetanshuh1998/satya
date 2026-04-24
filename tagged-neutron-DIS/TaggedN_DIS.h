#ifndef _TAGGEDNDIS_H
#define _TAGGEDNDIS_H 1

#include"TRandom3.h"
#include"TFile.h"
#include"TTree.h"
#include"TLorentzVector.h"
#include"TVector3.h"
#include"TString.h"

#include"KineCal.h"

#include"piIMParton.h"


class TaggedN_DIS{

	public:
		TaggedN_DIS();
		~TaggedN_DIS();

		/// the function to generate and dump N events into root file
		int Generate(int N);

		void SetOutputFileName(char *filename);
		void SetOutputFileName(TString filename);

		//// GRV's pionic PDF model, published in the year of 1992.
		double xuv_pi_GRV(double xpi, double Q2);
		double xdv_pi_GRV(double xpi, double Q2);
		double xsbar_pi_GRV(double xpi, double Q2);
		double F2_pi_GRV(double xpi, double Q2);
		double F2_pi_IMParton(double xpi, double Q2);
		double d4sigma_dQ2dxBdxLdt_GRV(double Q2, double xB, double xL, double t);


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
		double mpi;
		double mN;

		double PI;

		double xB;
		double xL;
		double Q2;
		double t;
		double xpi;
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

		TLorentzVector *elec_out;
		TLorentzVector *neut_out;

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
		piIMParton *pionPDFIMP;

		char *strFileName;

		void MakeROOTFile(char *filename);
};


#endif
