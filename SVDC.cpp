#include "svdc.h"
#include "svd.h"

SVD::SVD(Matrix B)
{
	A = B;
	svd(B, U, Sigma, VT);
}

int SVD::rank(double tol)
{
	int m = Sigma.getM(), n = Sigma.getN();
	int rank = 0;
	int limit = (m < n) ? m : n;
	for (int i = 1; i <= limit; i++)
	{
		if (Sigma(i, i) > tol)
		{
			rank++;
		}
	}
	return rank;
}

Matrix SVD::null(double tol)
{

	Matrix V = VT.transpose();
	int m = V.getM(), n = V.getN();
	int r = rank(tol);
	Matrix null(m, n - r);
	null.setblk(1, 1, V(1, m, r + 1, n));
	return null;
}

Matrix SVD::range(double tol)
{

	int m = U.getM(), n = U.getN();
	int r = rank(tol);
	Matrix range(m, r);
	range.setblk(1, 1, U(1, m, 1, r));
	return range;
}

Matrix SVD::sumRank1(int lower, int upper, double tol)
{
	int r = rank(tol);
	int m = A.getM(), n = A.getN();
	Matrix V = VT.transpose();
	Matrix E(m, n);
	upper = (r <= upper) ? r : upper;
	lower = (lower >= 1) ? lower : 1;
	if (upper >= lower)
	{
		for (int i = lower; i <= upper; i++)
		{
			E = E + Sigma(i, i) * U(1, m, i, i) * V(1, n, i, i).transpose();
		}
		return E;
	}
	return E;
}

Matrix SVD::pinv(double tol)
{
	int m = A.getM(), n = A.getN();
	int r = rank(tol);
	Matrix SigmaPinv(n, m);
	for (int i = 1; i <= r; i++)
	{
		SigmaPinv(i, i) = 1.0 / Sigma(i, i);
	}
	return VT.transpose() * SigmaPinv * U.transpose();
}

void SVD::polar(Matrix &Q, Matrix &S)
{
	Q = U * VT;
	S = VT.transpose() * Sigma * VT;
}

Matrix SVD::axb(Matrix b, bool &c, bool &u, Matrix &p, Matrix &Q, double tol)
{
	Matrix Ap = pinv();
	c = vNormInf(A * Ap * b - b) < tol;
	p = Ap * b;
	Q = (identity(A.getN()) - Ap * A);
	u = normInf(Q) < tol;
	return p;
}

double SVD::cond(double tol)
{
	int r = rank(tol);
	if (r == A.getM() && r == A.getN())
	{
		return Sigma(1, 1) / Sigma(r, r);
	}
	else
		return 1e99;
}

Matrix SVD::inv(double tol)
{
	Matrix SigmaInv(Sigma.getM(), Sigma.getN());
	int r = rank(tol);
	for (int i = 1; i <= r; i++)
	{
		SigmaInv(i, i) = 1.0 / Sigma(i, i);
	}
	if (r == A.getM() && r == A.getN())
	{
		return (VT.transpose() * SigmaInv * U.transpose());
	}
	else
		return identity(A.getM());
}