#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "matrix.h"

using namespace std;

void phase0(Matrix &A, Matrix &b, Matrix &c, Matrix s);

int phase2(Matrix A, Matrix b, Matrix c, vector<int> &Bi, vector<int> &Ni, Matrix &xB, Matrix &sP, int &t);

void phase1(Matrix A, Matrix b, Matrix c, Matrix s, vector<int> &Bi, vector<int> &Ni);

int simplex(Matrix A, Matrix b, Matrix c, Matrix s, vector<int> &Bi, vector<int> &Ni, Matrix &xB, Matrix &sP, int &t);

void lps(Matrix A, Matrix b, Matrix c, Matrix s, Matrix &x, double &minCost, int &status, Matrix &sP, int &t);

#endif
