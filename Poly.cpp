#include "poly.h"
#include "matrix.h"
#include "cmatrix.h"
#include "eig.h"
#include "eigen.h"

Poly::Poly(Matrix C)
{
	int m = C.getM(), n = C.getN();
	if (m == 1 || n == 1)
	{
		coef = C;
		deg = (m >= n) ? m - 1 : n - 1;
		isCol = (m >= n) ? true : false;
	}
	else
	{
		coef = Matrix(m, 1);
		deg = 0;
		isCol = true;
	}
}

void Poly::print()
{
	if (isCol)
	{
		coef.transpose().print();
	}
	else
	{
		coef.print();
	}
}

Poly operator*(double scalar, Poly P)
{
	return Poly(scalar * P.get_coef());
}

Matrix compan(Poly P)
{
	int n = P.get_deg();
	Matrix result(n, n);
	Matrix v = (P.get_isCol()) ? (P.get_coef()).transpose() : P.get_coef();
	v = (-1 / v(1)) * v;
	result.setblk(1, 1, v(1, 1, 2, n + 1));
	result.setblk(2, 1, identity(n - 1));
	return result;
}

Matrix compane(Poly P)
{
	int n = P.get_deg();
	Matrix result(n, n);
	Matrix v = (P.get_isCol()) ? (P.get_coef()).transpose() : P.get_coef();
	v = (-1 / v(1)) * v;
	Matrix u = v;
	for (int i = 1; i <= n + 1; i++)
	{
		v(i) = u(n + 2 - i);
	}
	result.setblk(n, 1, v(1, 1, 1, n));
	result.setblk(1, 2, identity(n - 1));
	return result;
}

void Poly::pzero()
{
	Matrix M = compane(*this);
	Eigen m(M);
	m.pAllVal();
}

cMatrix Poly::gzero(int k1, int k2)
{
	Matrix M = compane(*this);
	Eigen m(M);
	return m.allEigval(k1, k2);
}

Poly Poly::derivative()
{
	Matrix C(deg, 1, vector<double>(deg, 0));
	for (int i = 1; i <= deg; i++)
	{
		C(i) = (deg + 1 - i) * coef(i);
	}
	return Poly(C);
}

double Poly::value(double x)
{
	double s = coef(1);
	for (int i = 1; i <= deg; i++)
	{
		s = x * s + coef(i + 1);
	}
	return s;
};

double Poly::rzero(double x0, double tol)
{
	Poly dp = derivative();
	double x = x0;
	double x1 = 0;
	for (int i = 0; i < 100; i++)
	{
		x1 = x - value(x) / dp.value(x);
		if (abs(x1 - x) < tol)
		{
			return x1;
		}
		x = x1;
	}
	return x1;
}