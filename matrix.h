#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Matrix
{
protected:
	vector<double> vdbl;
	int m, n;

public:
	Matrix(const int m = 0, const int n = 0);/* Note 1*/
	Matrix(const int m, const int n, vector<double> vdbl);/* Note 2 */
	void print(int precision = 2, int width = 12);
	int getM() { return m; }
	int getN() { return n; }
	vector<double> getVdbl() { return vdbl; }
	Matrix transpose();
	double operator()(const int i, const int j) const;/* Note 3 */
	double &operator()(const int i, const int j);/* Note 4 */
	double operator()(const int k) const;
	double &operator()(const int k);
	Matrix getblk(int i1, int i2, int j1, int j2);
	void setblk(int r, int c, Matrix A);
	Matrix operator()(int i1, int i2, int j1, int j2);
	friend Matrix operator+(const Matrix &A, const Matrix &B);
	friend Matrix operator-(const Matrix &A, const Matrix &B);
	friend Matrix operator*(const double &scalar, const Matrix &A);
	friend Matrix operator*(const Matrix &A, const Matrix &B);
};

Matrix identity(int order);
double vNorm2(Matrix x);
double vNormInf(Matrix x);
double normInf(Matrix A);
double norm1(Matrix A);
double normF(Matrix A);
Matrix diag(Matrix v);
Matrix gDiag(Matrix A);
Matrix hilbert(int n);
Matrix hankel(Matrix C, Matrix R);
Matrix vander(Matrix C);
Matrix linspace(int n, double lower, double upper);
void save(const char *fileName, Matrix A);
Matrix load(const char *fileName);
Matrix normalize(Matrix A);
#endif
/* Note 1: All parameters of this function are given default values, therefore it is the default constructor. Thus declaration such as >Matrix A;< is valid. */

/* Note 2: Elements of >Matrix< are of type double and stored in an >vector<double><, a 1-dimensional container, in a row-major manner. Vectors are constructed and handled as >Matrices< with at least one of its dimensions being 1. */

/* Note 3: The function >double operator()(const int i, const int j) const< is explained below.
(1) Name : operator()(int, int) const.
(2) Use : Overloads the operator () to get the value of an element of a Matrix object.
(3) This is a const function which is valid for a member function only. The const qualifier ensuures that it is invlid to modify the calling object. Roughly speaking it makes the function read-only.
(4) As a member function it can be invoked as >double a = A(i,j)<, where A is a Matrix object.
(5) The parameters are const qualified to ensure it is invalid to modify them. */

/* Note 4: The function >double &operator()(const int i, const int j)< is explained below:
(1) Name : operator(int, int).
(2) Use : Overloads the operator () to assign value to an element of a Matrix object.
(3) As a member function it is invoked as A(i,j) similar to function explained in Note 3. The difference this time is this function returns a reference to the element and is used in an assignment statement such as A(i,j) = a.*/