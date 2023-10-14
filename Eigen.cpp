#include "eigen.h"
#include "eig.h"
#include "cmatrix.h"
#include <iostream>
#include "leqs.h"
#include "qr.h"
#include "qrc.h"
using namespace std;

Eigen::Eigen(Matrix B, double tol, int iter)
{
	A = B; // initialize Eigen::A

	keig(A, tol, allVal, iter); // keig(.) defined in Eig.cpp
	vSort(allVal);
	disVal = distinct(allVal); // distinct(.) defined in cMatrix.cpp
	eigIsReal = isReal(disVal);
	Lamda = cMatrix(A.getM(), A.getN()); // get complex Lamda
	for (int i = 1; i <= A.getM(); i++)
	{
		Lamda(i, i) = allVal(i);
	}

	if (!eigIsReal)
	{
		S = cMatrix(A.getM(), A.getN());
		cMatrix s = eigvecN(1, disVal.getM());
		S.setblk(1, 1, s);
		isDiag = (disVal.getN() == A.getM() || S.getN() == s.getN());
		Sinv = (isDiag) ? cInverse(S) : cIdentity(A.getM());
	}
	else
	{
		Matrix rDistVal = real(disVal);
		rS = Matrix(A.getM(),A.getN());
		Matrix rs = reigvecN(1, disVal.getM(), rDistVal);
		rS.setblk(1,1,rs);
		rS=normalize(rS);		
		if (isDiag)
		{
			QR q(rS);
			rSinv = q.inv();
		}
		else
		{
			rSinv = identity(A.getM());
		}

		rLamda = real(getLamda());
	}
}

void Eigen::pAllVal()
{
	for (int i = 1; i <= allVal.getM(); i++)
	{
		cout << i << "    ";
		allVal(i).print();
	}
	cout << endl;
}

cMatrix Eigen::allEigval(int k1, int k2)
{
	vector<C> cV;
	for (int i = k1; i <= k2 && i <= allVal.getM(); i++)
	{
		cV.push_back(allVal(i));
	}
	return cMatrix(cV.size(), 1, cV);
}

void Eigen::pDisVal()
{
	for (int i = 1; i <= disVal.getM(); i++)
	{
		cout << i << "    ";
		disVal(i).print();
	}
	cout << endl;
}

cMatrix Eigen::disEigval(int k1, int k2)
{
	vector<C> cV;
	for (int i = k1; i <= k2 && i <= disVal.getM(); i++)
	{
		cV.push_back(disVal(i));
	}
	return cMatrix(cV.size(), 1, cV);
}

cMatrix Eigen::eigvecN(int k1, int k2)
{
	cMatrix B = eigVal2VecN(A, disVal(k1));
	for (int i = k1 + 1; i <= k2 && i <= disVal.getM(); i++)
	{
		B = combine(B, eigVal2VecN(A, disVal(i)));
	}
	return (B);
}

cMatrix Eigen::eigvec(int k1, int k2)
{
	cMatrix B = eigVal2Vec(A, disVal(k1));
	for (int i = k1 + 1; i <= k2 && i <= disVal.getM(); i++)
	{
		B = combine(B, eigVal2Vec(A, disVal(i)));
	}
	return (B);
}

Matrix Eigen::reigvecN(int k1, int k2, Matrix rdistVal)
{
	Matrix B = reigVal2VecN(A, rdistVal(k1));
	for (int i = k1 + 1; i <= k2 && i <= rdistVal.getM(); i++)
	{
		B = combine(B, reigVal2VecN(A, rdistVal(i)));
	}
	return (B);
}

Matrix Eigen::reigVal2VecN(Matrix A, double eigval)
{
	int n = A.getN();
	A = A - (eigval * identity(n));
	return nullspace(A);
}

bool Eigen::realLamda(Matrix &B)
{ // get real Lamda if all eigen values are real
	if (isReal(Lamda))
	{
		B = real(Lamda);
		return true;
	}
	else
	{
		B = identity(A.getM());
		return false;
	}
}

bool Eigen::check(double tol)
{
	Matrix residue = real(S * Lamda * Sinv) - A;
	if (normInf(residue) < tol)
		return true;
	else
		return false;
}

double norm2(Matrix A)
{
	A = A.transpose() * A;
	Eigen a(A);
	return sqrt(a.allEigval(1, 1)(1).a);
}
