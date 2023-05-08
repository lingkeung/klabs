#ifndef PV_H
#define PV_H

#include "matrix.h"

class PV : public Matrix
{

public:
	PV() : Matrix(3, 1) {}
	PV(Matrix A) : Matrix(3, 1, vector<double>{A.getVdbl()[0], A.getVdbl()[1], A.getVdbl()[2]}){};
	PV(vector<double> vdbl) : Matrix(3, 1, vdbl) {}
	PV(double a1, double a2, double a3) : Matrix(3, 1, vector<double>{a1, a2, a3}) {}
	PV(double m, double a, double b, double g);

	void set_vdbl(vector<double> vdbl) { this->vdbl = vdbl; }
	void set(double a1, double a2, double a3);
	void mdc(double &mag, double &cosa, double &cosb, double &cosg);
	void mabg(double &m, double &a, double &b, double &g);
	PV unit();

	friend PV operator+(const PV &A, const PV &B);
	friend PV operator-(const PV &A, const PV &B);
	friend PV operator*(const double &scalar, const PV &A);
	friend PV operator/(const double &scalar, const PV &A);
	friend double operator%(const PV &u, const PV &v); // dot product
	friend PV operator&(const PV &u, const PV &v);	   // vector product

private:
	double r2d(double rad)
	{
		const double pi = 3.14159265358979323846;
		return rad / pi * 180;
	}

	double d2r(double deg)
	{
		const double pi = 3.14159265358979323846;
		return deg / 180 * pi;
	}
};

#endif
