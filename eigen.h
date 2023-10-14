#ifndef EIGEN_H
#define EIGEN_H
#include "matrix.h"
#include "cmatrix.h"

class Eigen
{
private:
	Matrix A;
	cMatrix allVal;
	cMatrix disVal;
	bool isDiag;
	bool eigIsReal;
	cMatrix S;
	cMatrix Lamda;
	cMatrix Sinv;
	Matrix rS;
	Matrix rLamda;
	Matrix rSinv;

public:
	Eigen(){};
	Eigen(Matrix B, double tol = 1e-8, int iter = 50);
	void pAllVal();
	void pDisVal();
	Matrix getA() { return A; }
	cMatrix get_allVal() { return allVal; }
	cMatrix get_disVal() { return disVal; }
	bool get_isDiag() { return isDiag; }
	cMatrix getS() { return S; }
	cMatrix getSinv() { return Sinv; }
	cMatrix getLamda() { return Lamda; }
	Matrix getrS() { return rS; }
	Matrix getrSinv() { return rSinv; }
	Matrix getrLamda() { return rLamda; }
	bool getEigIsReal(){return eigIsReal;}
	bool realLamda(Matrix &B);
	cMatrix allEigval(int k1, int k2);
	cMatrix disEigval(int k1, int k2);
	cMatrix eigvecN(int k1, int k2);
	cMatrix eigvec(int k1, int k2);
	Matrix reigVal2VecN(Matrix A, double eigval);
	Matrix reigvecN(int k1, int k2, Matrix rdistVal);
	bool check(double tol = 1e-3);
};

double norm2(Matrix A);
#endif
