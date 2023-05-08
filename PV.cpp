#include "pv.h"
#include "cmath"

PV::PV(double mag, double a, double b, double g)
{
	m = 3;
	n = 1;
	vdbl = vector<double>{mag * cos(d2r(a)), mag * cos(d2r(b)), mag * cos(d2r(g))};
}

void PV::set(double a1, double a2, double a3)
{
	vector<double> vdbl{a1, a2, a3};
	set_vdbl(vdbl);
}

PV operator+(const PV &A, const PV &B)
{
	PV result;
	for (int i = 1; i <= 3; i++)
	{
		result(i) = A(i) + B(i);
	}
	return result;
}

PV operator-(const PV &A, const PV &B)
{
	PV result;
	for (int i = 1; i <= 3; i++)
	{
		result(i) = A(i) - B(i);
	}
	return result;
}

PV operator*(const double &scalar, const PV &A)
{
	PV result;
	for (int i = 1; i <= 3; i++)
	{
		result(i) = scalar * A(i);
	}
	return result;
}

PV operator/(const double &scalar, const PV &A)
{
	PV result;
	for (int i = 1; i <= 3; i++)
	{
		result(i) = A(i) / scalar;
	}
	return result;
}

double operator%(const PV &u, const PV &v)
{
	double result = 0;
	for (int i = 1; i <= 3; i++)
	{
		result += u(i) * v(i);
	}
	return result;
}

PV operator&(const PV &u, const PV &v)
{
	PV r;
	r(1) = u(2) * v(3) - u(3) * v(2);
	r(2) = u(3) * v(1) - u(1) * v(3);
	r(3) = u(1) * v(2) - u(2) * v(1);
	return r;
}

void PV::mdc(double &mag, double &cosa, double &cosb, double &cosg)
{
	mag = vNorm2(*this);
	cosa = (*this)(1) / mag;
	cosb = (*this)(2) / mag;
	cosg = (*this)(3) / mag;
}

void PV::mabg(double &m, double &a, double &b, double &g)
{
	mdc(m, a, b, g);
	a = r2d(acos(a));
	b = r2d(acos(b));
	g = r2d(acos(g));
}

PV PV::unit()
{
	double m, ca, cb, cg;
	mdc(m, ca, cb, cg);
	return PV(ca, cb, cg);
}