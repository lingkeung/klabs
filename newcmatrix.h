#ifndef NEWCMATRIX_H
#define NEWCMATRIX_H

// #include "matrix.h"
#include <vector>
#include <complex>

using namespace std;

typedef complex<double> C;

class newcMatrix
{
public:
    vector<complex<double>> vC;
    int m, n;

public:
    newcMatrix(const int m = 0, const int n = 0);
    newcMatrix(const int m, const int n, vector<C> vC);
    void print();
    /*int getM() { return m; }
    int getN() { return n; }
    vector<C> getvC() { return vC; }
    cMatrix hermitian();*/
    C operator()(const int i, const int j) const;
    C &operator()(const int i, const int j);
    C operator()(const int k) const;
    C &operator()(const int k);
    
    newcMatrix getblk(int i1, int i2, int j1, int j2);
    void setblk(int r, int c, newcMatrix A);
    newcMatrix operator()(int i1, int i2, int j1, int j2);
    friend newcMatrix operator+(const newcMatrix &A, const newcMatrix &B);
    friend newcMatrix operator-(const newcMatrix &A, const newcMatrix &B);
    friend newcMatrix operator*(const C &scalar, const newcMatrix &A);
    friend newcMatrix operator*(const double &scalar, const newcMatrix &A);
    friend newcMatrix operator*(const newcMatrix &A, const newcMatrix &B);
};
/*
cMatrix cIdentity(int order);
cMatrix cplex(Matrix A);
double vNorm2(cMatrix x);
void vNormal(cMatrix &v);
double vNormInf(cMatrix x);
double normInf(cMatrix A);
void vSort(cMatrix &v);
cMatrix distinct(cMatrix v);
cMatrix combine(cMatrix A, cMatrix B);
Matrix real(cMatrix A);
bool isReal(cMatrix A);
cMatrix ndft(Matrix x); // naive discrete fourier transform
Matrix nidft(cMatrix X); // naive inverse discrete fourier transform
cMatrix fdft(Matrix x); // fast discrete fourier transform (radix-2)
cMatrix fidftHelp(cMatrix x); // helper fast discrete fourier transform, N = 2^n
Matrix fidft(cMatrix X);  */

#endif