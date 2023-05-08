#ifndef PLUC_H
#define PLUC_H

class PLU
{
private:
	Matrix A;
	Matrix P;
	Matrix L;
	Matrix U;
	double det;
	bool isInv;
	bool isChk;

public:
	PLU() {}
	PLU(Matrix B);
	Matrix getA() { return A; }
	Matrix getP() { return P; }
	Matrix getL() { return L; }
	Matrix getU() { return U; }
	double get_det() { return det; }
	bool get_isInv() { return isInv; }
	bool get_isChk() { return isChk; }
	Matrix axb(Matrix b);
	Matrix inv();
	bool check(double tol = 1e-8);
};

#endif