#ifndef _KINECAL_H
#define _KINECAL_H 1

class KineCal
{
	public:
		KineCal();
		~KineCal();

		double calTMin(double s, double m1_sq, double m2_sq, double m3_sq, double m4_sq);
		double calTMax(double s, double m1_sq, double m2_sq, double m3_sq, double m4_sq);
		double calW2(double xB, double Q2);
		double calXB(double W2, double Q2);
		double calXi(double xB, double m_M, double Q2);

		double calY(double s, double xB, double Q2); /// fractional energy loss
		double calEpsilon(double y);
		double calEpsilon(double y, double Q2, double s);
		double calEpsilonW2Q2(double W2, double Q2, double s);

	private:


};

#endif
