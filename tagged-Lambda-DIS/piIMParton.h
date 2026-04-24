#ifndef IMPARTON_H_
#define IMPARTON_H_

class piIMParton
//this class is for users
{
public:
	piIMParton();   
	virtual void setDataSet(int);                            //Choose a data set, 1 is for set A and 2 is for set B
	virtual double getPDF(int, double x, double Q2) const;   //user function to get parton distribution functions, see ReadMe.txt for details       
        virtual double getXUV(double x, double Q2) const;        //return x(u -ubaar)
        virtual double getXDV(double x, double Q2) const;        //return x(d - dbar)
        virtual double getXUSea(double x, double Q2) const;      //return 2x*ubar
        virtual double getXDSea(double x, double Q2) const;      //return 2x*dbar
        virtual double getXSSea(double x, double Q2) const;      //return 2x*sbar
        virtual double getXCSea(double x, double Q2) const;      //return 2x*cbar
        virtual double getXGluon(double x, double Q2) const;     //return x*gluon
	virtual ~piIMParton(void);                                 //deconstructor function

private:
	unsigned int xMax;        //grid points number for variable x
	unsigned int Q2Max;       //grid points number for variable Q^2
	unsigned int flavorMax;   //flavor number in the grid data
	double lnxstep;           //step length of ln(x) for the grid data
	double lnQ2step;          //step length of ln(Q^2) for the grid data

	int hadronID;
	int dataSetNO;

        double * grid;            //grid data array
        double * gridA;           //data array storing data set A
        double * gridB;           //data array storing data set B
        double * gridPion;           //data array storing data set of the pion from the global fit
        double * gridKaon;           //data array storing data set of the kaon from the global fit

	double fitLinear(double x, double* px, double* pf) const;       //linear interpolation function
	double fitQuadratic(double x, double* px, double* pf) const;    //quadratic interpolation function
        double getPDFType(int, double x, double Q2) const;
        double getPDFraw(int, double x, double Q2) const;
  
};

#endif
