#include "plu.h"
#include "pluc.h"

PLU::PLU(Matrix B)
{
	A = B;
	isInv = plu(A, P, L, U, det);
	isChk = check();
}

Matrix PLU::axb(Matrix b)
{
	Matrix bc = b;
	if (isInv)
		return bsub(U, fsub(L, P * b));
	else
	{
		cout << "Matrix may be singular!" << endl;
		return bc;
	}
}

Matrix PLU::inv()
{
	int n = A.getN();
	Matrix I = identity(n);
	if (isInv)
	{
		for (int i = 1; i <= n; i++)
		{
			I.setblk(1, i, bsub(U, fsub(L, P * I(1, n, i, i))));
		}
		return I;
	}
	else
	{
		cout << "Matrix may be singular" << endl;
		return A;
	}
}

bool PLU::check(double tol)
{
	bool result = false;
	if (normInf(P * A - L * U) < tol)
		result = true;
	return result;
}
