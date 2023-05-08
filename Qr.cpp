#include "qr.h"
#include "plu.h"

void house(Matrix x, Matrix &v, double &beta)
{
	int n = x.getM();
	if (vNorm2(x) < 1e-15)
	{ // handles small or zero x
		v = Matrix(n, 1);
		beta = 0;
		return;
	}

	double eta = vNormInf(x);
	x = (1 / eta) * x;
	double sigma = (x(2, n, 1, 1).transpose() * x(2, n, 1, 1))(1);
	v = Matrix(n, 1);
	v.setblk(2, 1, x(2, n, 1, 1));
	if (sigma <= 1e-10)
		beta = 0;
	else
	{
		double alpha = sqrt(x(1) * x(1) + sigma);
		if (abs(x(1)) <= 1e-15)
			v(1) = x(1) - alpha;
		else
			v(1) = -1 * sigma / (x(1) + alpha);

		beta = 2 * v(1) * v(1) / (sigma + v(1) * v(1));
		v = (1 / v(1)) * v;
	}
}

void qr(Matrix A, Matrix &Q, Matrix &R, Matrix &Q1, Matrix &R1)
{
	Matrix Ac = A;
	int m = A.getM();
	int n = A.getN();
	Matrix d(n, 1);
	for (int j = 1; j <= n; j++)
	{
		if (j < m)
		{
			Matrix v;
			double beta = 0;
			house(A(j, m, j, j), v, beta);
			A.setblk(j, j, (identity(m - j + 1) - beta * v * v.transpose()) * A(j, m, j, n));
			d(j) = beta;
			A.setblk(j + 1, j, v(2, m - j + 1, 1, 1));
		}
	}
	R = Matrix(m, n);
	for (int i = 1; i <= m; i++)
	{
		for (int j = i; j <= n; j++)
		{
			R(i, j) = A(i, j);
		}
	}
	Matrix V = A;
	for (int i = 1; i <= n; i++)
	{
		for (int j = i; j <= n; j++)
		{
			V(i, j) = (i == j) ? 1 : 0;
		}
	}
	Q = identity(m);
	Matrix Im = identity(m);
	for (int j = 1; j <= n; j++)
	{
		Q = (Im - d(j) * V(1, m, j, j) * V(1, m, j, j).transpose()) * Q;
	}
	Q = Q.transpose();

	Q1 = Q(1, m, 1, n);

	R1 = R(1, n, 1, n);
}

Matrix qrAxb(Matrix A, Matrix b)
{
	Matrix Q, R, Q1, R1;
	qr(A, Q, R, Q1, R1);
	Matrix c1 = Q1.transpose() * b;
	return bsub(R1, c1);
}

Matrix qrBasis(Matrix A)
{
	Matrix Q, R, Q1, R1;
	qr(A, Q, R, Q1, R1);
	return Q;
}

Givens::Givens(double c, double s)
{
	this->c = c;
	this->s = s;
	rot = Matrix(2, 2, vector<double>({this->c, this->s, -1 * this->s, this->c}));
}

Givens givens(double a, double b)
{
	double c, s;

	if (abs(b) <= 1e-10)
	{
		c = 1;
		s = 0;
	}
	else
	{
		if (abs(b) > abs(a))
		{
			double tau = 1 * a / b;
			s = 1 / (sqrt(1 + tau * tau));
			c = s * tau;
		}
		else
		{
			double tau = 1 * b / a;
			c = 1 / (sqrt(1 + tau * tau));
			s = c * tau;
		}
	}
	Givens result(c, s);
	return result;
}

void grm(Givens G, int i, int k, Matrix &A)
{
	int m = A.getM();
	for (int j = 1; j <= m; j++)
	{
		double tau1 = A(j, i);
		double tau2 = A(j, k);
		A(j, i) = G.c * tau1 + G.s * tau2;
		A(j, k) = -1 * G.s * tau1 + G.c * tau2;
	}
}

void glm(Givens G, int i, int k, Matrix &A)
{
	int n = A.getN();
	for (int j = 1; j <= n; j++)
	{
		double tau1 = A(i, j);
		double tau2 = A(k, j);
		A(i, j) = G.c * tau1 + G.s * tau2;
		A(k, j) = -1 * G.s * tau1 + G.c * tau2;
	}
}

void qrGivens(Matrix A, Matrix &Q, Matrix &R, Matrix &Q1, Matrix &R1)

{
	int m = A.getM(), n = A.getN();
	Givens g;
	Matrix QT = identity(m);
	for (int j = 1; j <= n; j++)
	{
		for (int i = m; i >= (j + 1); i--)
		{
			g = givens(A(i - 1, j), A(i, j));
			Matrix mat = A(i - 1, i, j, n);
			glm(g, 1, 2, mat);
			A.setblk(i - 1, j, mat);
			glm(g, i - 1, i, QT);
		}
	}
	Q = QT.transpose();
	R = A;
	Q1 = Q(1, m, 1, n);
	R1 = R(1, n, 1, n);
}