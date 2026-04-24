#include"KineCal.h"
#include<cmath>

KineCal::KineCal(){}
KineCal::~KineCal(){}

double KineCal::calTMin(double s, double m1_sq, double m2_sq, double m3_sq, double m4_sq){
	double part1 = pow((m1_sq-m3_sq-m2_sq+m4_sq)/2.0/sqrt(s), 2);
	double E1cm = (s+m1_sq-m2_sq)/2.0/sqrt(s);
	double E3cm = (s+m3_sq-m4_sq)/2.0/sqrt(s);
	double p1cm = sqrt(E1cm*E1cm-m1_sq);
	double p3cm = sqrt(E3cm*E3cm-m3_sq);
	return part1 - pow(p1cm-p3cm, 2);
}

double KineCal::calTMax(double s, double m1_sq, double m2_sq, double m3_sq, double m4_sq){
	double part1 = pow((m1_sq-m3_sq-m2_sq+m4_sq)/2.0/sqrt(s), 2);
	double E1cm = (s+m1_sq-m2_sq)/2.0/sqrt(s);
	double E3cm = (s+m3_sq-m4_sq)/2.0/sqrt(s);
	double p1cm = sqrt(E1cm*E1cm-m1_sq);
	double p3cm = sqrt(E3cm*E3cm-m3_sq);
	return part1 - pow(p1cm+p3cm, 2);
}

double KineCal::calW2(double xB, double Q2){
	double mN = 0.938272;
	return (1.0/xB-1)*Q2 + mN*mN;
}

double KineCal::calXB(double W2, double Q2){
	double mN = 0.938272;
	return Q2/(Q2+W2-mN*mN);
}

double KineCal::calXi(double xB, double m_M, double Q2){
	return xB/(2-xB)*(1+m_M*m_M/Q2);
}

double KineCal::calY(double s, double xB, double Q2){
	double mN = 0.938272;
	/// fractional energy loss
	return Q2/xB/(s-mN*mN);
}


double KineCal::calEpsilon(double y){
	return 2.0*(1-y)/(1+(1-y)*(1-y)); 
}
double KineCal::calEpsilon(double y, double Q2, double s){
	double mN = 0.938272;
	double fourE2 = (s-mN*mN) * (s-mN*mN) /mN/mN;
	return (1 - y - Q2/fourE2) / (1 - y + y*y/2.0 + Q2/fourE2);
}
double KineCal::calEpsilonW2Q2(double W2, double Q2, double s){
	double xb = calXB(W2, Q2);
	double y = calY(s, xb, Q2);  
	return calEpsilon(y, Q2, s);
}


