#ifndef PLU_H
#define PLU_H

#include "matrix.h"

Matrix fsub(Matrix L, Matrix b);

Matrix bsub(Matrix U, Matrix b);

void rexch(int i, int j, Matrix &A);

bool plu(Matrix A, Matrix &P, Matrix &L, Matrix &U, double &det);

Matrix axb(Matrix A, Matrix b);

double det(Matrix A);

Matrix inverse(Matrix A);

Matrix null(Matrix A);

vector<int> icols(Matrix A); //find independent columns

#endif