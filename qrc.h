#ifndef QRC_H
#define QRC_H
#include "matrix.h"

using namespace std;

class QR
{
	Matrix A;
	Matrix Q;
	Matrix Q1;
	Matrix R;
	Matrix R1;

public:
	QR() {}
	QR(Matrix A);
	QR(char q, Matrix A);
	Matrix getA() { return A; }
	Matrix getQ() { return Q; }
	Matrix getQ1() { return Q1; }
	Matrix getR() { return R; }
	Matrix getR1() { return R1; }
	Matrix axb(Matrix b);
	Matrix basis();
	Matrix inv();

private:
	Matrix Rinv();
};

#endif