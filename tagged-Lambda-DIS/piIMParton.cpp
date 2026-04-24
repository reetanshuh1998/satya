#include<iostream>
#include<fstream>
#include<string>
extern "C"{
#include<math.h>
}
#include "piIMParton.h"
using namespace std;


//a method used to choose a data set
void piIMParton::setDataSet(int dataset)
{
        if(dataset==1)
        {
                grid = gridA;
		hadronID = 211;
		dataSetNO = dataset;
                cout<<"    Using data set A."<<endl;
        }
        else if(dataset==2)
        {
                grid = gridB;
		hadronID = 211;
		dataSetNO = dataset;
                cout<<"    Using data set B."<<endl;
        }
        else if(dataset==211)
        {
                grid = gridPion;
		hadronID = 211;
		dataSetNO = dataset;
                cout<<"    Using data set of the pion from the global fit."<<endl;
        }
        else if(dataset==321)
        {
                grid = gridKaon;
		hadronID = 321;
		dataSetNO = dataset;
                cout<<"    Using data set of the kaon from the global fit."<<endl;
        }
        else
        {
                cout<<"!!->Unknown data set "<<dataset<<"."<<endl;
                cout<<"!!->Data set should be 211, or 321, or 1, or 2"<<endl;
        }
}

//return the parton distributions of different kinds at x and Q^2
//Here we use the numbering scheme of Monte Carlo suggested by particle data group
double piIMParton::getPDF(int Iparton, double x, double Q2) const
{
	if(hadronID==211){
		if(Iparton==1)return getPDFraw(-2, x, Q2);
		else if(Iparton==-1)return getPDFraw(2, x, Q2);
		else if(Iparton==2)return getPDFraw(1, x, Q2);
		else if(Iparton==-2)return getPDFraw(-1, x, Q2);
		else if(Iparton==3 || Iparton==-3)return getPDFraw(Iparton, x, Q2);
		else if(Iparton==4 || Iparton==-4)return getPDFraw(Iparton, x, Q2);
		else if(Iparton==21)return getPDFraw(0, x, Q2);
		else 
        	{
                	cout<<"!!->Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                	cout<<"!!->Iparton should be one of these: [-4,-3, ... 3, 4, 21]"<<endl;
                	cout<<"!!->Using PDG MC numbering scheme."<<endl;
                	return 0;
        	}
	}
	else if(hadronID==321){
		if(Iparton==1)return getPDFraw(-2, x, Q2);
		else if(Iparton==-1)return getPDFraw(-2, x, Q2);
		else if(Iparton==2)return getPDFraw(1, x, Q2);
		else if(Iparton==-2)return getPDFraw(-1, x, Q2);
		else if(Iparton==3)return getPDFraw(-3, x, Q2);
		else if(Iparton==-3)return getPDFraw(2, x, Q2) - getPDFraw(-2, x, Q2) + getPDFraw(-3, x, Q2);
		else if(Iparton==4 || Iparton==-4)return getPDFraw(Iparton, x, Q2);
		else if(Iparton==21)return getPDFraw(0, x, Q2);
		else 
        	{
                	cout<<"!!->Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                	cout<<"!!->Iparton should be one of these: [-4,-3, ... 3, 4, 21]"<<endl;
                	cout<<"!!->Using PDG MC numbering scheme."<<endl;
                	return 0;
        	}
	}
	else{
		cout<<"!!Unknown hadron."<<endl;
		return 0;
	}
}

//the constructor and initialization
piIMParton::piIMParton()
{
	cout<<"    piIMParton version - 1.0, for pion & kaon parton distribution functions"<<endl;
	char filename[50];
	ifstream datain;
	double x, Q2;
	unsigned int i, j;
	xMax=61;
	Q2Max=32;
	flavorMax=7;
	lnxstep=log(10)/((xMax-1)/6);
	lnQ2step=log(2.0);
	gridA = new double[Q2Max*xMax*flavorMax];
	gridB = new double[Q2Max*xMax*flavorMax];
	gridPion = new double[Q2Max*xMax*flavorMax];
	gridKaon = new double[Q2Max*xMax*flavorMax];
	//read grid data for interpolation
	//reading data set of the pion from the global fit
	sprintf(filename,"grid_pion_global_fit.dat");
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridPion+(xMax*i+j)*7))>>(*(gridPion+(xMax*i+j)*7+1))>>(*(gridPion+(xMax*i+j)*7+2))>>(*(gridPion+(xMax*i+j)*7+3))>>(*(gridPion+(xMax*i+j)*7+4))>>(*(gridPion+(xMax*i+j)*7+5))>>(*(gridPion+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	//reading data set of the kaon from the global fit
	sprintf(filename,"grid_kaon_global_fit.dat");
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridKaon+(xMax*i+j)*7))>>(*(gridKaon+(xMax*i+j)*7+1))>>(*(gridKaon+(xMax*i+j)*7+2))>>(*(gridKaon+(xMax*i+j)*7+3))>>(*(gridKaon+(xMax*i+j)*7+4))>>(*(gridKaon+(xMax*i+j)*7+5))>>(*(gridKaon+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	//reading data set A
	sprintf(filename,"grid_pion_SetA.dat");
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridA+(xMax*i+j)*7))>>(*(gridA+(xMax*i+j)*7+1))>>(*(gridA+(xMax*i+j)*7+2))>>(*(gridA+(xMax*i+j)*7+3))>>(*(gridA+(xMax*i+j)*7+4))>>(*(gridA+(xMax*i+j)*7+5))>>(*(gridA+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	//reading data set B
        sprintf(filename,"grid_pion_SetB.dat");
        cout<<"    Loading "<<filename<<endl;
        datain.open(filename);
        if(!datain.good())cout<<"!!->Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
        else
        for(i=0;i<Q2Max;i++)
        {
                for(j=0;j<xMax;j++)
                {
                        datain>>Q2>>x>>(*(gridB+(xMax*i+j)*7+0))>>(*(gridB+(xMax*i+j)*7+1))>>(*(gridB+(xMax*i+j)*7+2))>>(*(gridB+(xMax*i+j)*7+3))>>(*(gridB+(xMax*i+j)*7+4))>>(*(gridB+(xMax*i+j)*7+5))>>(*(gridB+(xMax*i+j)*7+6));
                }
        }
        datain.close();
	//the default is set of the pion from the global fit
	grid = gridPion;

	hadronID = 211;
	dataSetNO = 211;

}

//the deconstructor
piIMParton::~piIMParton(void)
{
	delete[] gridA;
        delete[] gridB;
        delete[] gridPion;
        delete[] gridKaon;

}

//a method which returns xuv
double piIMParton::getXUV(double x, double Q2) const
{
	return getPDFType(1,x,Q2); 
}

//a method which returns xdv
double piIMParton::getXDV(double x, double Q2) const
{
	return getPDFType(2,x,Q2);
}

//a method which returns xusea
double piIMParton::getXUSea(double x, double Q2) const
{
	return getPDFType(3,x,Q2);
}

//a method which returns xdsea
double piIMParton::getXDSea(double x, double Q2) const
{
	return getPDFType(4,x,Q2);
}

//a method which returns xssea
double piIMParton::getXSSea(double x, double Q2) const
{
        return getPDFType(5,x,Q2);
}

//a method which returns xcsea
double piIMParton::getXCSea(double x, double Q2) const
{
        return getPDFType(6,x,Q2);
}

//a method which returns xgluon
double piIMParton::getXGluon(double x, double Q2) const
{
        return getPDFType(0,x,Q2);
}

//a method which returns different types of distributions
double piIMParton::getPDFType(int Iparton, double x, double Q2) const
{
	if(Iparton<0 || Iparton>6)
	{
		cout<<"!!->Wrong Iparton input for getPDFType(int Iparton, double x, double Q2)."<<endl;
		return 0;
	}
	else
	{
		double lnx, lnQ2;
		int i=(int)(lnx=log(x*1e6)/lnxstep);
		int j=(int)(lnQ2=log(Q2*8)/lnQ2step);
		double g0[3], g1[3], g2[3], g[3]={0};
		if(i<0)i=0;
		if(i>(int)(xMax-3))i=xMax-3;
		if(j<0)j=0;
		if(j>29)j=29;
		//avoid log(1-x) calculation in below algorithm
		if(x>0.9999)return 0.0;
		//if x>0.5, we use A(1-x)^B to do the interpolation
		else if(x>0.5 && Iparton!=6)
		{
			double vln1_x[2]={log(1-exp(log(1e-6)+i*lnxstep)),log(1-exp(log(1e-6)+(i+1)*lnxstep))};
			g0[0]=log(grid[(xMax*j+i)*7+Iparton]);
			g0[1]=log(grid[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(grid[(xMax*j+i)*7+Iparton]);
			g1[1]=log(grid[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(grid[(xMax*j+i)*7+Iparton]);
			g2[1]=log(grid[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(log(1-x),vln1_x,g0));
			g[1]=exp(fitLinear(log(1-x),vln1_x,g1));
			g[2]=exp(fitLinear(log(1-x),vln1_x,g2));
		}
		//if x<1e-5, we use A*x^B to do the interpolation
		//for valance quark, B>0; for gluon and sea quark, B<0
		else if(x<1e-4)
		{
			double vlnx[2]={(double)i,(double)(i+1)};
			g0[0]=log(grid[(xMax*j+i)*7+Iparton]);
			g0[1]=log(grid[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(grid[(xMax*j+i)*7+Iparton]);
			g1[1]=log(grid[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(grid[(xMax*j+i)*7+Iparton]);
			g2[1]=log(grid[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(lnx,vlnx,g0));
			g[1]=exp(fitLinear(lnx,vlnx,g1));
			g[2]=exp(fitLinear(lnx,vlnx,g2));
		}
		//we use quadratic interpolation method for other situations
		else
		{
			double vlnx[3]={(double)i,(double)(i+1),(double)(i+2)};
			g0[0]=grid[(xMax*j+i)*7+Iparton];
			g0[1]=grid[(xMax*j+i+1)*7+Iparton];
			g0[2]=grid[(xMax*j+i+2)*7+Iparton];
			j++;
			g1[0]=grid[(xMax*j+i)*7+Iparton];
			g1[1]=grid[(xMax*j+i+1)*7+Iparton];
			g1[2]=grid[(xMax*j+i+2)*7+Iparton];
			j++;
			g2[0]=grid[(xMax*j+i)*7+Iparton];
			g2[1]=grid[(xMax*j+i+1)*7+Iparton];
			g2[2]=grid[(xMax*j+i+2)*7+Iparton];
			g[0]=fitQuadratic(lnx,vlnx,g0);
			g[1]=fitQuadratic(lnx,vlnx,g1);
			g[2]=fitQuadratic(lnx,vlnx,g2);
		}
		//if Q2>1, we do the interpolation to the variable ln(Q^2)
		if(Q2>1)
		{
			double vlnQ2[3]={(double)(j-2),(double)(j-1),(double)j};
			return fitQuadratic(lnQ2,vlnQ2,g);
		}
		//if Q2<1, we do the interpolation to the variable Q^2
		else 
		{
			double vQ2[3]={0.125*pow(2,j-2),0.125*pow(2,j-1),0.125*pow(2,j)};
			return fitQuadratic(Q2,vQ2,g);
		}	
	}
}


//quadratic interpolation method
double piIMParton::fitQuadratic(double x, double* px, double* pf) const
{
	double f01=(pf[1]-pf[0])/(px[1]-px[0]);
	double f12=(pf[2]-pf[1])/(px[2]-px[1]);
	double f012=(f12-f01)/(px[2]-px[0]);
	return pf[0]+f01*(x-px[0])+f012*(x-px[0])*(x-px[1]);
}

//linear interpolation method
double piIMParton::fitLinear(double x, double* px, double* pf) const
{
	double f01=(pf[1]-pf[0])/(px[1]-px[0]);
	return pf[0]+f01*(x-px[0]);
}


//return the parton distributions of different kinds at x and Q^2
double piIMParton::getPDFraw(int Iparton, double x, double Q2) const
{
        if(Iparton==-4 || Iparton==4)return getXCSea(x,Q2)/2.0/x;
        else if(Iparton==-3 || Iparton==3)return getXSSea(x,Q2)/2.0/x;
        else if(Iparton==-2)return getXDSea(x,Q2)/2.0/x;
        else if(Iparton==2)return getXDSea(x,Q2)/2.0/x+getXDV(x,Q2)/x;
        else if(Iparton==-1)return getXUSea(x,Q2)/2.0/x;
        else if(Iparton==1)return getXUSea(x,Q2)/2.0/x+getXUV(x,Q2)/x;
        else if(Iparton==0)return getXGluon(x,Q2)/x;
        else
        {
                cout<<"!!->Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                cout<<"!!->Iparton should be one of these: [-4,-3, ... 3, 4]"<<endl;
                return 0;
        }

}


