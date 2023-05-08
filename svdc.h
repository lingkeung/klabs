#ifndef SVDC_H
#define SVDC_H
#include "matrix.h"

class SVD
{
private:
	Matrix A;
	Matrix U;
	Matrix Sigma;
	Matrix VT;

public:
	SVD(){};
	SVD(Matrix A);
	Matrix getA() { return A; }
	Matrix getU() { return U; }
	Matrix getSigma() { return Sigma; }
	Matrix getVT() { return VT; }

	int rank(double tol = 1e-8);
	Matrix null(double tol = 1e-8);
	Matrix range(double tol = 1e-8);
	Matrix sumRank1(int k1, int k2, double tol = 1e-8);
	Matrix pinv(double tol = 1e-8);
	void polar(Matrix &Q, Matrix &S);
	Matrix axb(Matrix b, bool &c, bool &u, Matrix &P, Matrix &Q, double tol = 1e-8);
	double cond(double tol = 1e-8);
	Matrix inv(double tol = 1e-8);
};
#endif