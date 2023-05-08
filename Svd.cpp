#include "svd.h"
#include "qr.h"
#include "plu.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

void bidiag(Matrix &A, Matrix &UT, Matrix &V)
{

	int n = A.getN();
	int m = A.getM();
	UT = identity(m);
	V = identity(n);
	for (int k = 1; k <= n; k++)
	{
		Matrix x = A(k, m, k, k);
		Matrix v;
		double beta;
		house(x, v, beta);
		A.setblk(k, k, (identity(m - k + 1) - beta * v * v.transpose()) * A(k, m, k, n));
		Matrix vm(m, 1);
		for (int i = 0; i < v.getM(); i++)
		{
			vm(m - i) = v(v.getM() - i);
		}
		UT = (identity(m) - beta * vm * vm.transpose()) * UT;
		if (k <= n - 2)
		{
			Matrix x2 = A(k, k, k + 1, n).transpose();
			Matrix v2;
			double beta2;
			house(x2, v2, beta2);
			A.setblk(k, k + 1, A(k, m, k + 1, n) * (identity(n - k) - beta2 * v2 * v2.transpose()));
			Matrix v2m(n, 1);
			for (int i = 0; i < v2.getM(); i++)
			{
				v2m(n - i) = v2(v2.getM() - i);
			}
			V = V * (identity(n) - beta2 * v2m * v2m.transpose());
		}
	}
}

void wilkinson(Matrix &B, Matrix &PT, Matrix &Q)
{
	int m = B.getM(), n = B.getN();

	double dn, gn_1, dn_1, gn_2;
	dn = B(n, n);
	gn_1 = B(n - 1, n - 2);
	dn_1 = B(n - 1, n - 1);
	gn_2 = B(n - 2, n - 1);

	double d1 = B(1, 1);
	double g1 = B(1, 2);

	double alpha = dn * dn + gn_1 * gn_1;
	double delta = (dn_1 * dn_1 + gn_2 * gn_2 - alpha) / 2;
	double beta = dn_1 * gn_1;
	double signdel;
	if (delta >= 0)
	{
		signdel = 1;
	}
	else
	{
		signdel = -1;
	}
	double mu = alpha - beta * beta / (delta + signdel * sqrt(delta * delta + beta * beta));
	double y = d1 * d1 - mu;
	double z = d1 * g1;

	PT = identity(m);
	Q = identity(n);

	for (int k = 1; k <= n - 1; k++)
	{
		Givens G;
		G = givens(y, z);
		grm(G, k, k + 1, B);
		y = B(k, k);
		z = B(k + 1, k);
		grm(G, k, k + 1, Q);
		Givens G1;
		G1 = givens(y, z);
		glm(G1, k, k + 1, B);
		glm(G1, k, k + 1, PT);
		if (k < n - 1)
		{
			y = B(k, k + 1);
			z = B(k, k + 2);
		}
	}
}

void pqn(Matrix &B, int &i1, int &i2, int &j1, int &j2, int &m)
{

	double tol = 1e-8;
	int r = B.getM();
	for (int i = 1; i <= r - 1; i++)
	{
		if (abs(B(i, i + 1)) <= tol * (abs(B(i, i)) + abs(B(i + 1, i + 1))))
		{
			B(i, i + 1) = 0;
		}
		if (abs(B(i, i)) <= tol * normInf(B))
		{
			B(i, i) = 0;
		}
	}
	m = 0;

	int row = r; // ensures initially m = r - row = 0
	for (int i = r - 1; i >= 1; i--)
	{
		if (abs(B(i, i + 1)) <= tol)
		{
			row = i;
		}
		else
		{
			break;
		}
	}
	m = r - row;
	if (m == r - 1)
	{
		m = r;
	}
	int nonzero = 0;
	int row1 = row;
	for (int j = row - 1; j >= 1; j--)
	{
		if (abs(B(j, j + 1)) >= tol)
		{
			row1 = j;
			nonzero++;
		}
		else
			break;
	}
	int n = nonzero + 1;

	i1 = row1;
	i2 = n + i1 - 1;
	j1 = row1;
	j2 = n + j1 - 1;
}

void svd42(Matrix &A, Matrix &U, Matrix &SIGMA, Matrix &VT)
{
	// handles 2x2 matrices only
	int m = A.getM(), n = A.getN();
	U = identity(m);
	VT = identity(n);
	Givens G; // user-defined Givens class
	for (int i = 0; i < 10; i++)
	{								  // direct reduction to diagonal matrix ...
		G = givens(A(2, 2), A(1, 2)); // ... using Givens rotations
		G.s = -1 * G.s;
		glm(G, 1, 2, A); // glm(.) and grm(.) used instead of matrix multiplication
		glm(G, 1, 2, U);
		G = givens(A(2, 2), A(2, 1));
		G.s = -1 * G.s;
		grm(G, 1, 2, A);
		grm(G, 1, 2, VT);
	}
	U = U.transpose();
	SIGMA = A;
	for (int j = 1; j <= n; j++)
	{ // ensures non-negative singular values
		if (SIGMA(j, j) < 0)
		{
			SIGMA(j, j) = -1 * SIGMA(j, j);
			for (int i = 1; i <= m; i++)
			{
				U(i, j) = -1 * U(i, j);
			}
		}
	}
	vector<double> sv; // sort singular values in decending order
	for (int i = 1; i <= SIGMA.getN(); i++)
	{
		sv.push_back(SIGMA(i, i));
	}
	vector<int> idx(SIGMA.getN());
	std::iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(), [&](int i, int j)
		 { return sv[i] > sv[j]; });
	Matrix sU = U;
	Matrix sV = VT;
	Matrix sSIGMA = SIGMA;
	for (int i = 0; i < SIGMA.getN(); i++)
	{														  // based on A*U = V*SIGMA using ...
		sSIGMA(i + 1, i + 1) = SIGMA(idx[i] + 1, idx[i] + 1); // ... the same permutation for SIGMA, U and V
		sU.setblk(1, i + 1, U(1, U.getM(), idx[i] + 1, idx[i] + 1));
		sV.setblk(1, i + 1, VT(1, VT.getM(), idx[i] + 1, idx[i] + 1));
	}
	U = sU; // final 'outputs'
	SIGMA = sSIGMA;
	VT = sV.transpose();
}

void svd(Matrix &A, Matrix &U, Matrix &SIGMA, Matrix &VT)
{
	// notation used : UT*A*V=SIGMA or A=U*SIGMA*VT
	int m = A.getM(), n = A.getN();
	if (m == 2 && n == 2)
	{ // svd42(.) handles 2x2 matrices
		svd42(A, U, SIGMA, VT);
		return;
	}
	if (n > m)
	{							   // recursive call handles n > m cases
		Matrix AT = A.transpose(); // for AT m >= n is met
		svd(AT, U, SIGMA, VT);	   // recursive call
		SIGMA = SIGMA.transpose(); // A=U*SIGMA*VT -> AT=V*SIGMAT*UT
		Matrix Uc = U;
		U = VT.transpose();
		VT = Uc.transpose();
		return;
	}
	// core code to handle m >= n cases except for 2x2
	bidiag(A, U, VT);  // transformation to bidiagonal matrix
	U = U.transpose(); // bidiag(.) 'returns' UT and V; svd(.) does U and VT
	int i1, i2, j1, j2, q;
	int itrn;
	int t = 50; // usually increasing t does not help
	for (itrn = 0; itrn < t; itrn++)
	{
		// cout << "iteration " << itrn << endl;//****
		// cout << "A = " << endl;//*****
		// A.print(10,15);// *****
		pqn(A, i1, i2, j1, j2, q); // pqn(.) analyses state of reduction of bidiagonal matrices
								   // cout << i1 << j1 << q << endl;
		if (q == m)
			break; // ****** signals reduction completed ******

		Matrix A22 = A(i1, i2, j1, j2); // diagonal block of A that needs furthur reduction
		int n22 = A22.getN();
		Matrix In22 = identity(n22);
		int zero = 0;
		Givens G;
		// finds zero diagonal elements and zero the rest of the corresponding rows
		for (int i = 1; i < n22 - 1; i++)
		{
			if (abs(A22(i + 1, i + 1)) <= 1e-12)
			{
				zero++;
				// cout << "zero = " << zero << endl;
				for (int k = 1; k < n22 - i; k++)
				{
					G = givens(A22(i + k + 1, i + k + 1), A22(i + 1, i + k + 1));
					G.s = -1 * G.s;
					glm(G, i + 1, i + k + 1, A22);	// modifies A22
					glm(G, i + 1, i + k + 1, In22); // captures matrix operations
				}
			}
		}
		if (zero != 0)
		{							// at least 1 zero diagonal element detected and row zeroing done
			A.setblk(i1, j1, A22);	// inserts modified A22 back into A
			Matrix I = identity(m); // constructs new U
			I.setblk(i1, j1, In22);
			U = U * I;
		}
		else
		{ // ****core of core code: reduces A; constructs U and VT with similar procedures
			Matrix PT22;
			Matrix Q22;
			// cout << "A22 in = " << endl;
			// A22.print();
			wilkinson(A22, PT22, Q22); // **** core modification using wilkinson(.) ****
									   // cout << "A22 out = " << endl;
									   // A22.print();
			A.setblk(i1, j1, A22);
			Matrix Im = identity(m);
			Im.setblk(i1, j1, PT22.transpose()); // captures corresponding matrix operations
			U = U * Im;
			Matrix In = identity(n);
			In.setblk(i1, j1, Q22);
			VT = VT * In; // generates new VT
		}
	}
	SIGMA = A; // *** At this stage A has been transformed to the final SIGMA ***
	for (int j = 1; j <= n; j++)
	{ // ensures SIGMA with non-negative values only
		if (SIGMA(j, j) < 0)
		{
			SIGMA(j, j) = -1 * SIGMA(j, j);
			for (int i = 1; i <= m; i++)
			{ // corresponding modifications to U (A=(U*SIGMA)*VT)
				U(i, j) = -1 * U(i, j);
			}
		}
	}
	vector<double> sv; // sorts results such that the singular values appear in decending order
	for (int i = 1; i <= SIGMA.getN(); i++)
	{
		sv.push_back(SIGMA(i, i));
	}
	// std::iota below requires <numeric>
	vector<int> idx(SIGMA.getN());		  // eventually contains indices of SIGMA elements ...
	std::iota(idx.begin(), idx.end(), 0); // ... according to decending element values;
	sort(idx.begin(), idx.end(), [&](int i, int j)
		 { return sv[i] > sv[j]; }); // <algorithm> req'd
	Matrix sU = U;
	Matrix sV = VT;
	Matrix sSIGMA = SIGMA;
	for (int i = 0; i < SIGMA.getN(); i++)
	{																 // based on the relation A*V = U*SIGMA
		sSIGMA(i + 1, i + 1) = SIGMA(idx[i] + 1, idx[i] + 1);		 // SIGMA elements and columns of V and U ...
		sU.setblk(1, i + 1, U(1, U.getM(), idx[i] + 1, idx[i] + 1)); //... permutated in the same way
		sV.setblk(1, i + 1, VT(1, VT.getM(), idx[i] + 1, idx[i] + 1));
	}
	U = sU; 
	SIGMA = sSIGMA;
	VT = sV.transpose();
}