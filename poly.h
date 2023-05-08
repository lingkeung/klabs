#ifndef POLY_H
#define POLY_H
#include "matrix.h"
#include "cmatrix.h"

class Poly
{
private:
	Matrix coef;
	bool isCol;
	int deg;

public:
	Poly() {}
	Poly(Matrix C);
	void print();
	Matrix get_coef() { return coef; }
	bool get_isCol() { return isCol; }
	int get_deg() { return deg; };
	void pzero();
	cMatrix gzero(int k1, int k2);
	double rzero(double x0, double tol);
	Poly derivative();
	double value(double x);
	friend Poly operator*(double scalar, Poly P);
};

Matrix compan(Poly P);
Matrix compane(Poly P);

#endif