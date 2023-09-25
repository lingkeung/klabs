#include "qrc.h"
#include "qr.h"
#include "plu.h"

/* If A is square and non-singular, then */
/* A = QR where Q is orthogonal and R is upper triangular */

QR::QR(Matrix A)/* use Householder reflection */
{
	qr(A, Q, R, Q1, R1);
}

QR::QR(char q, Matrix A) /* use Givens rotation */
{
	qrGivens(A, Q, R, Q1, R1);
}

Matrix QR::axb(Matrix b) /* A=QR & Ax=b -> QRx=b -> Rx=QTb -> x=bsub(R,QTb) */
{
	Matrix c1 = Q1.transpose() * b;
	return bsub(R1, c1);
}

Matrix QR::basis() /* A=QR->A's columns are linear combinations of Q's columns */
{
	return getQ();
}

Matrix QR::Rinv()/* A=QR -> Ainv=RinvQinv=RinvQT */
{
	int m = R.getM(), n = R.getN();
	Matrix result = identity(m);
	if (m != n)
		return result;
	for (int j = 1; j <= n; j++)/* RRinv=I -> Rinv=bsub(R,I) */
	{
		result.setblk(1, j, bsub(R, result(1, m, j, j)));
	}
	return result;
}

Matrix QR::inv() /* A=QR -> Ainv=RinvQinv=RinvQT */
{
	if (A.getM() == A.getN())
	{
		return (Rinv() * Q.transpose());
	}
	else
		return identity(A.getM());
}