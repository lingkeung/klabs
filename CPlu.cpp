#include "cplu.h"
#include <cmath>
#include <iostream>

using namespace std;

cMatrix fsub(cMatrix L, cMatrix b)
{
	auto m = L.getM();
	for (int j = 1; j <= m - 1; j++)
	{
		b(j) = b(j) / L(j, j);
		b.setblk(j + 1, 1, b(j + 1, m, 1, 1) - b(j) * L(j + 1, m, j, j));
	}
	b(m) = b(m) / L(m, m);
	return b;
}

cMatrix bsub(cMatrix U, cMatrix b)
{
	auto n = U.getN();
	for (auto j = n; j >= 2; j--)
	{
		b(j) = b(j) / U(j, j);
		b.setblk(1, 1, b(1, j - 1, 1, 1) - b(j) * U(1, j - 1, j, j));
	}
	b(1) = b(1) / U(1, 1);
	return b;
}

void rexch(int i, int j, cMatrix &A)
{
	int n = A.getN();
	cMatrix temp = A(i, i, 1, n);
	A.setblk(i, 1, A(j, j, 1, n));
	A.setblk(j, 1, temp);
}

/*double abs(C cn) {
	double a = cn.a;
	double b = cn.b;
	double r = sqrt(a*a + b*b);
	return r;
}*/

bool cPlu(cMatrix A, cMatrix &P, cMatrix &L, cMatrix &U, C &det)
{
	int n = A.getN();
	int p = 0;
	P = cIdentity(n);
	det = C(1, 0);

	for (int k = 1; k <= n - 1; k++)
	{
		double max = 0;
		for (int i = k; i <= n; i++)
		{
			if (abs(A(i, k)) > max)
			{
				max = abs(A(i, k));
				p = i;
			}
		}
		if (p != k)
		{
			rexch(p, k, A);
			rexch(p, k, P);
			det = C(-1, 0) * det;
		}

		if (abs(A(k, k)) >= 1e-6)
		{
			A.setblk(k + 1, k, (C(1, 0) / A(k, k)) * A(k + 1, n, k, k));
			A.setblk(k + 1, k + 1, (A(k + 1, n, k + 1, n) - A(k + 1, n, k, k) * A(k, k, k + 1, n)));
		}
		else
		{
			return false;
		}
	}

	L = cIdentity(n);
	for (int i = 2; i < n + 1; i++)
	{
		for (int j = 1; j < i; j++)
		{
			L(i, j) = A(i, j);
		}
	}
	U = cMatrix(n, n);
	for (int i = 1; i < n + 1; i++)
	{
		for (int j = n; j >= i; j--)
		{
			if (j >= i)
			{
				U(i, j) = A(i, j);
			}
		}
	}

	for (int i = 1; i <= n; i++)
	{
		det = det * U(i, i);
	}
	if (abs(det) == 0)
		return false;
	return true;
}

cMatrix cAxb(cMatrix A, cMatrix b)
{
	cMatrix bc = b;
	cMatrix P, L, U;
	C det;
	if (cPlu(A, P, L, U, det))
		return bsub(U, fsub(L, P * b));
	else
	{
		cout << "Matrix may be singular!" << endl;
		return bc;
	}
}

C det(cMatrix A)
{
	cMatrix P, L, U;
	C det;
	if (cPlu(A, P, L, U, det))
		return det;
	else
	{
		cout << "Matrix may be singular !" << endl;
		return C(0, 0);
	}
}

cMatrix cInverse(cMatrix A)
{
	int n = A.getN();
	cMatrix I = cIdentity(n);
	cMatrix P, L, U;
	C det;
	if (cPlu(A, P, L, U, det))
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

cMatrix cNull(cMatrix A)
{
	int m = A.getM(), n = A.getN();
	int size = (m <= n) ? m : n;
	int p = 0;
	for (int k = 1; k <= size - 1; k++)
	{
		double max = 0;
		for (int i = k; i <= size; i++)
		{
			if (abs(A(i, k)) > max)
			{
				max = abs(A(i, k));
				p = i;
			}
		}
		if (p != k)
		{
			rexch(p, k, A);
		}

		if (abs(A(k, k)) >= 1e-6)
		{
			A.setblk(k + 1, k, (C(1, 0) / A(k, k)) * A(k + 1, n, k, k));
			A.setblk(k + 1, k + 1, A(k + 1, n, k + 1, n) - A(k + 1, n, k, k) * A(k, k, k + 1, n));
		}
		for (int i = k + 1; i <= m; i++)
		{
			A(i, k) = C(0, 0);
		}
	}
	int rank = 0;
	Matrix u(size, 1);
	for (int i = 1; i <= size; i++)
	{
		if (abs(A(i, i)) >= 1e-6)
		{
			rank++;
			u(i) = 1;
		}
	}

	vector<C> vU;
	vector<int> eqN;
	for (int i = 1; i <= size; i++)
	{
		if (u(i) == 1)
		{
			eqN.push_back(i);
			for (int j = 1; j <= n; j++)
			{
				if (u(j) == 1)
				{
					vU.push_back(A(i, j));
				}
			}
		}
	}
	cMatrix U(rank, rank, vU);
	vector<C> vC;
	vector<cMatrix> vb;
	cMatrix soln(n, size - rank);
	int k = 1;
	for (int i = 1; i <= size; i++)
	{
		if (u(i) == 0)
		{
			soln(i, k) = C(1, 0);
			k++;
			for (int idx = 0; idx < eqN.size(); idx++)
			{
				vC.push_back(-1 * A(eqN[idx], i));
			}
			cMatrix mat(rank, 1, vC);
			vb.push_back(mat);
			vC.clear();
		}
	}
	for (int i = 0; i < vb.size(); i++)
	{
		for (int j = 0; j < eqN.size(); j++)
		{
			soln(eqN[j], i + 1) = bsub(U, vb[i])(j + 1);
		}
	}
	return soln;
}