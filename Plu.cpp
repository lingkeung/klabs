#include "plu.h"

Matrix fsub(Matrix L, Matrix b)
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

Matrix bsub(Matrix U, Matrix b)
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

void rexch(int i, int j, Matrix &A)
{
	int n = A.getN();
	Matrix temp = A(i, i, 1, n);
	A.setblk(i, 1, A(j, j, 1, n));
	A.setblk(j, 1, temp);
}

bool plu(Matrix A, Matrix &P, Matrix &L, Matrix &U, double &det)
{
	int n = A.getN();
	int p = 0;
	P = identity(n);
	det = 1;
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
			det *= -1;
		}
		if (abs(A(k, k)) >= 1e-6)
		{
			A.setblk(k + 1, k, (1 / A(k, k)) * A(k + 1, n, k, k));
			A.setblk(k + 1, k + 1, (A(k + 1, n, k + 1, n) - (A(k + 1, n, k, k) * A(k, k, k + 1, n))));
		}
		else
		{
			return false;
		}
	}

	L = identity(n);
	for (int i = 2; i < n + 1; i++)
	{
		for (int j = 1; j < i; j++)
		{
			L(i, j) = A(i, j);
		}
	}
	U = Matrix(n, n);
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
		det *= U(i, i);
	}
	if (det == 0)
		return false;
	return true;
}

Matrix axb(Matrix A, Matrix b)
{
	Matrix bc = b;
	Matrix P, L, U;
	double det;
	if (plu(A, P, L, U, det))
		return bsub(U, fsub(L, P * b));
	else
	{
		cout << "Matrix may be singular!" << endl;
		return bc;
	}
}

double det(Matrix A)
{
	Matrix P, L, U;
	double det;
	if (plu(A, P, L, U, det))
		return det;
	else
	{
		cout << "Matrix may be singular !" << endl;
		return 0;
	}
}

Matrix inverse(Matrix A)
{
	int n = A.getN();
	Matrix I = identity(n);
	Matrix P, L, U;
	double det;
	if (plu(A, P, L, U, det))
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

Matrix null(Matrix A)
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
			A.setblk(k + 1, k, (1 / A(k, k)) * A(k + 1, n, k, k));
			A.setblk(k + 1, k + 1, A(k + 1, n, k + 1, n) - A(k + 1, n, k, k) * A(k, k, k + 1, n));
		}
		for (int i = k + 1; i <= m; i++)
		{
			A(i, k) = 0;
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

	vector<double> vU;
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
	Matrix U(rank, rank, vU);
	vector<double> vdbl;
	vector<Matrix> vb;
	Matrix soln(n, size - rank);
	int k = 1;
	for (int i = 1; i <= size; i++)
	{
		if (u(i) == 0)
		{
			soln(i, k) = 1;
			k++;
			for (int idx = 0; idx < eqN.size(); idx++)
			{
				vdbl.push_back(-1 * A(eqN[idx], i));
			}
			Matrix mat(rank, 1, vdbl);
			vb.push_back(mat);
			vdbl.clear();
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

vector<int> icols(Matrix A)
{
	vector<int> result;
	int m = A.getM(), n = A.getN();

	int p = 0;
	for (int k = 1; k <= m - 1; k++)
	{
		double max = 0;
		for (int i = k; i <= m; i++)
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
			A.setblk(k + 1, k, (1 / A(k, k)) * A(k + 1, m, k, k)); // store multipliers (L)
			A.setblk(k + 1, k + 1, A(k + 1, m, k + 1, n) - A(k + 1, m, k, k) * A(k, k, k + 1, n));
			for (int i = k + 1; i <= m; i++)
			{
				A(i, k) = 0; // restore zeros (erase multipliers)
			}
			result.push_back(k);
		}
	}

	for (int i = 0; i <= n - m; i++)
	{
		if (abs(A(m, m + i)) >= 1e-6)
		{
			result.push_back(m + i); // check for last column pivot
			break;
		}
	}

	// A.print(); // for debugging
	return result;
}
