#ifndef SVD_H
#define SVD_H

#include "matrix.h"

void bidiag(Matrix &A, Matrix &UT, Matrix &V);

void wilkinson(Matrix &B, Matrix &PT, Matrix &Q);

void pqn(Matrix &B, int &i1, int &i2, int &j1, int &j2, int &m);

void svd42(Matrix &A, Matrix &U, Matrix &SIGMA, Matrix &VT);

void svd(Matrix &A, Matrix &U, Matrix &SIGMA, Matrix &VT);

#endif