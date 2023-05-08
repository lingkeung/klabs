#include "simplex.h"
#include "matrix.h"
#include "iostream"
#include "qr.h"
#include "qrc.h"
#include <cmath>
#include "plu.h"

using namespace std;

void phase0(Matrix &A, Matrix &b, Matrix &c, Matrix s)
{ // convert to standardized LP (all ='s)

	int m = A.getM(), n = A.getN();
	Matrix Ac = A, cc = c;
	int k = m;
	for (int i = 1; i <= m; i++)
	{ // s is the sign matrix; -1/0/1 means >=/=/<= constraint
		if (s(i) == 0)
			k--; // get number of slack variables required (k <= m)
	}

	A = Matrix(m, n + k); // set new A
	A.setblk(1, 1, Ac);
	int offset = 0;
	for (int i = 1; i <= m; i++)
	{ // set sign of slack variables
		if (s(i) == 0)
			offset++;
		else
			A(i, n + i - offset) = s(i);
	}

	c = Matrix(1, n + k); // set new c
	c.setblk(1, 1, cc);
	c.setblk(1, n + 1, Matrix(1, k));

	b = b; // no change for b
}

int phase2(Matrix A, Matrix b, Matrix c, vector<int> &Bi, vector<int> &Ni, Matrix &xB, Matrix &sP, int &tout)
{

	int m = A.getM(), n = A.getN();
	// prepare various elements of the simplex method
	int output = 1;
	for (int t = 1; t <= 100; t++)
	{
		Matrix cB(1, m), cN(1, n - m);
		Matrix B(m, m);
		for (int i = 1; i <= m; i++)
		{
			B.setblk(1, i, A(1, m, Bi[i - 1], Bi[i - 1]));
			cB.setblk(1, i, c(1, 1, Bi[i - 1], Bi[i - 1]));
		}
		Matrix N(m, n - m);
		for (int i = 1; i <= n - m; i++)
		{
			N.setblk(1, i, A(1, m, Ni[i - 1], Ni[i - 1]));
			cN.setblk(1, i, c(1, 1, Ni[i - 1], Ni[i - 1]));
		}

		QR q(B);			// QR used for its stability;
		auto Q = q.getQ(); // solve lamda = cB * (B inverse)
		auto R = q.getR();
		auto P = fsub(R.transpose(), cB.transpose());
		auto lamda = (Q * P).transpose();
		Matrix r = cN - lamda * N;
		double Nmin = 0;
		int Nidx = 0;
		bool rIsPos = true;
		for (int i = 1; i <= r.getN(); i++)
		{ // pick varibale from Ni to enter Bi
			if (r(i) < Nmin)
			{
				Nmin = r(i);
				Nidx = i - 1;
				rIsPos = false;
			}
		}

		xB = q.axb(b); // solve xB = (B inverse) * b
		if (rIsPos)
		{ // check if no more cost reduction is possible
			output = 0;
			sP = lamda; // sP stands for shadow price; lamda stands for Lagrange multipliers
			tout = t;
			for (int i = 1; i <= xB.getM(); i++)
			{
				if (xB(i) != xB(i))
				{
					return 1;
				}
				if (xB(i) < -1e-10)
				{
					return 1;
				}
				else
					xB(i) = abs(xB(i));
			}
			break;
		}

		Matrix u = N(1, m, Nidx + 1, Nidx + 1);
		Matrix v = q.axb(u); // solve v = (B inverse) * u;
		double Bmin = 1e100; // a large number
		int Bidx = 0;
		bool isNInf = true;
		for (int i = 1; i <= v.getM(); i++)
		{ // pick variable to leave Bi
			if (v(i) > 0)
			{
				isNInf = false;
				if (((xB)(i) / v(i)) < Bmin)
				{
					Bmin = (xB)(i) / v(i);
					Bidx = i - 1;
				}
			}
		}
		if (isNInf)
		{ // check if min is -Inf
			output = 1;
			cout << " isNInf = 1 " << endl;
			break;
		}
		int temp = Ni[Nidx]; // exchange entering and leaving variable
		Ni[Nidx] = Bi[Bidx];
		Bi[Bidx] = temp;
	}
	return output; // signals execution status
}

void phase1(Matrix A, Matrix b, Matrix c, Matrix s, vector<int> &Bi, vector<int> &Ni)
{

	int m = A.getM(), n = A.getN();

	for (int i = 1; i <= m; i++)
	{ // get equivalent LP with nonnegative b elements
		if (b(i) < 0)
		{
			A.setblk(i, 1, -1 * A(i, i, 1, n));
			b(i) = -1 * b(i);
		}
	}

	Matrix Ac = A; // construct alternate LP : min sum(w) s.t. Ax + w = b and x,w >= 0
	A = Matrix(m, n + m);
	A.setblk(1, 1, Ac);
	A.setblk(1, n + 1, identity(m));
	Matrix cc = c;
	c = Matrix(1, n + m);
	c.setblk(1, n + 1, Matrix(1, m, vector<double>(m, 1)));

	for (int i = 1; i <= n + m; i++)
	{ // construct Bi and Ni of a feasible solution of alternate LP
		if (i < n + 1)
		{
			Ni.push_back(i);
		}
		else
		{
			Bi.push_back(i);
		}
	}

	Matrix xB, sP;
	int t;
	phase2(A, b, c, Bi, Ni, xB, sP, t); // get optimal solution of alternate LP hence Bi of original LP

	Ni.clear(); // construct Ni of original LP

	for (int i = 1; i <= n; i++)
	{
		bool isInBi = false;
		for (int j = 0; j < Bi.size(); j++)
		{
			if (i == Bi[j])
			{
				isInBi = true;
				break;
			}
		}
		if (!isInBi)
			Ni.push_back(i);
	}
	// Bi and Ni of one feasible solution for orginal LP are now ready
}

int simplex(Matrix A, Matrix b, Matrix c, Matrix s, vector<int> &Bi, vector<int> &Ni, Matrix &xB, Matrix &sP, int &t)
{

	phase0(A, b, c, s); // get standardized LP

	phase1(A, b, c, s, Bi, Ni); // get Bi and Ni of one basic and feasible solution for standardized LP

	int output = phase2(A, b, c, Bi, Ni, xB, sP, t); // get Bi, Ni and xB for optimal solution of original LP

	return output; // status of execution : 1 means abnormal; 0 means normal
}

void lps(Matrix A, Matrix b, Matrix c, Matrix s, Matrix &x, double &minCost, int &status, Matrix &sP, int &t)
{

	int m = A.getM();
	int n = A.getN();
	Matrix xB;
	vector<int> Bi, Ni;
	status = simplex(A, b, c, s, Bi, Ni, xB, sP, t); // see comments in simplex(.) code
	x = Matrix(n, 1);								 // matrix holding results
	for (int i = 1; i <= m; i++)
	{ // get x
		if (Bi[i - 1] <= n)
		{
			x(Bi[i - 1]) = (xB(i));
		}
	}
	minCost = (c * x)(1); // get minCost

	Matrix br = A * x; // equality of br and b within reasonal numerical accuray implies feasible solution

	bool status1 = 0; // get feasibility status (in addition to execution status)
	for (int i = 1; i <= m; i++)
	{
		if (s(i) == 1 && !((abs(br(i) - b(i)) < 1e-6) || br(i) < b(i)))
		{
			status1 = 1;
			break;
		}
		else if (s(i) == -1 && !((abs(br(i) - b(i)) < 1e-6) || br(i) > b(i)))
		{
			status1 = 1;
			break;
		}
		else if (s(i) == 0 && !(abs(br(i) - b(i)) < 1e-6))
		{
			status1 = 1;
			break;
		}
	}

	status = (status || status1); // 0 means usable result; 1 means result is not usable
}
