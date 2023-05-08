#ifndef EIG_H
#define EIG_H

#include "matrix.h"
#include "cmatrix.h"
#include "cplu.h"

Matrix qrSchur(Matrix A, int iter);
Matrix hessenberg(Matrix A);
void deflate(Matrix &H, double tol);
void lmn(Matrix H, int &i1, int &i2, int &j1, int &j2, int &m);
void francis(Matrix &H);
void schur2eig(Matrix Schur, cMatrix &eig);
void keig(Matrix A, double tol, cMatrix &eig, int iter = 50);
cMatrix eigVal2Vec(Matrix A, C eigval);
cMatrix eigVal2VecN(Matrix A, C eigval);

#endif