#include "matrix.h"
#include <numeric>
#include <algorithm>
#include <fstream>

Matrix::Matrix(const int m, const int n) : m(m), n(n)
{
	if (m <= 0 || n <= 0)
		return;
	vdbl = vector<double>(m * n, 0);
}

Matrix::Matrix(const int m, const int n, vector<double> vdbl)
{
	if (m * n == int(vdbl.size()))
	{
		this->m = m;
		this->n = n;
		this->vdbl = vdbl;
	}
	else
	{
		this->m = m;
		this->n = n;
		this->vdbl = vector<double>(m * n, 0);
	}
}

void Matrix::print(int precision, int width)
{
	cout << endl
		 << fixed << setprecision(precision);
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << setw(width) << (*this)(i, j);
		cout << endl;
	}
}

double Matrix::operator()(const int i, const int j) const
{
	int k = (i - 1) * n + (j - 1);
	return vdbl[k];
}

double &Matrix::operator()(const int i, const int j)
{
	int k = (i - 1) * n + (j - 1);
	return vdbl[k];
}

double Matrix::operator()(const int k) const
{
	return vdbl[k - 1];
}

double &Matrix::operator()(const int k)
{
	return vdbl[k - 1];
}

Matrix Matrix::getblk(int i1, int i2, int j1, int j2)
{
	int m = abs(i2 - i1 + 1);
	int n = abs(j2 - j1 + 1);
	Matrix result(m, n);
	if (i1 > i2 || j1 > j2)
	{
		return result;
	}
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			result(i, j) = (*this)(i + i1 - 1, j + j1 - 1);
		}
	}
	return result;
}

void Matrix::setblk(int r, int c, Matrix A)
{
	if ((m - r + 1 < A.m) || (n - c + 1 < A.n) || r < 1 || c < 1)
	{
		return;
	}
	for (int i = 1; i <= A.m; i++)
	{
		for (int j = 1; j <= A.n; j++)
		{
			(*this)(i + r - 1, j + c - 1) = A(i, j);
		}
	}
}

Matrix Matrix::operator()(int i1, int i2, int j1, int j2)
{
	return getblk(i1, i2, j1, j2);
}

Matrix operator+(const Matrix &A, const Matrix &B)
{
	Matrix C(A.m, A.n);
	if (A.m != B.m || A.n != B.n)
		return C;
	for (int i = 1; i <= A.m; i++)
		for (int j = 1; j <= A.n; j++)
			C(i, j) = A(i, j) + B(i, j);
	return C;
}

Matrix operator-(const Matrix &A, const Matrix &B)
{
	Matrix C(A.m, A.n);
	if (A.m != B.m || A.n != B.n)
		return C;
	for (int i = 1; i <= A.m; i++)
		for (int j = 1; j <= A.n; j++)
			C(i, j) = A(i, j) - B(i, j);
	return C;
}

Matrix operator*(const double &scalar, const Matrix &A)
{
	Matrix B(A.m, A.n);
	for (int i = 1; i <= A.m; i++)
		for (int j = 1; j <= A.n; j++)
			B(i, j) = scalar * A(i, j);
	return B;
}

Matrix operator*(const Matrix &A, const Matrix &B)
{
	Matrix C(A.m, B.n);
	if (A.n != B.m)
		return C;
	for (int i = 1; i <= A.m; i++)
	{
		for (int j = 1; j <= B.n; j++)
		{
			double sum = 0;
			for (int k = 1; k <= A.n; k++)
			{
				sum += (A(i, k) * B(k, j));
			}
			C(i, j) = sum;
		}
	}
	return C;
}

Matrix Matrix::transpose()
{
	Matrix B(this->n, this->m);
	for (int i = 1; i <= this->n; i++)
		for (int j = 1; j <= this->m; j++)
			B(i, j) = (*this)(j, i);
	return B;
}

Matrix identity(int order)
{
	Matrix A(order, order);
	for (auto i = 1; i <= order; i++)
	{
		A(i, i) = 1.0;
	}
	return A;
}

double vNorm2(Matrix x)
{
	if (x.getM() != 1 && x.getN() != 1)
	{
		return 1;
	}
	int size = x.getM() >= x.getN() ? x.getM() : x.getN();
	double result = 0;
	for (int i = 1; i <= size; i++)
	{
		result += x(i) * x(i);
	}
	return sqrt(result);
}

double vNormInf(Matrix x)
{
	if (x.getM() != 1 && x.getN() != 1)
	{
		return 1;
	}
	int size = x.getM() >= x.getN() ? x.getM() : x.getN();
	double result = 0;
	for (int i = 1; i <= size; i++)
	{
		if (abs(x(i)) >= result)
		{
			result = abs(x(i));
		}
	}
	return result;
}

double normInf(Matrix A)
{
	double result = 0;
	for (int i = 1; i <= A.getM(); i++)
	{
		double sum = 0;
		for (int j = 1; j <= A.getN(); j++)
		{
			sum += abs(A(i, j));
		}
		if (sum >= result)
			result = sum;
	}
	return result;
}

double norm1(Matrix A)
{
	double result = 0;
	for (int j = 1; j <= A.getN(); j++)
	{
		double sum = 0;
		for (int i = 1; i <= A.getM(); i++)
		{
			sum += abs(A(i, j));
		}
		if (sum >= result)
			result = sum;
	}
	return result;
}

double normF(Matrix A)
{
	vector<double> vdbl = A.getVdbl();
	std::for_each(vdbl.begin(), vdbl.end(), [](double &x)
				  { x *= x; });
	return sqrt(std::accumulate(vdbl.begin(), vdbl.end(), 0));
}

Matrix diag(Matrix v)
{
	int m = v.getM(), n = v.getN();
	if (m != 1 && n != 1)
		return identity(m);
	int r = (m >= n) ? m : n;
	Matrix result(r, r);
	for (int i = 1; i <= r; i++)
	{
		result(i, i) = v(i);
	}
	return result;
}

Matrix gDiag(Matrix A)
{
	int m = A.getM(), n = A.getN();
	int r = (m <= n) ? m : n;
	Matrix result(r, 1);
	for (int i = 1; i <= r; i++)
	{
		result(i) = A(i, i);
	}
	return result;
}

Matrix hilbert(int n)
{
	Matrix h(n, n);
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			h(i, j) = 1.0 / (i + j - 1);
		}
	}
	return h;
}

Matrix hankel(Matrix C, Matrix R)
{
	int m1 = C.getM(), n1 = C.getN(), m2 = R.getM(), n2 = R.getN();
	if (n1 != 1 || m2 != 1 || C(m1) != R(1))
		return Matrix(m1, n2);
	Matrix result(m1, n2);
	result.setblk(1, 1, C);
	result.setblk(m1, 1, R);
	for (int i = 1; i < m1; i++)
	{
		for (int j = 2; j <= n2; j++)
		{
			result(i, j) = (j - 1 >= m1 - i) ? R(j - m1 + i) : C(i + j - 1);
		}
	}
	return result;
}

Matrix vander(Matrix C)
{
	int m = C.getM(), n = C.getN();
	if (n != 1)
		return Matrix(m, m);
	Matrix result(m, m);
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			result(i, j) = pow(C(i), m - j);
		}
	}
	return result;
}

Matrix linspace(int n, double lower, double upper)
{
	Matrix result(n, 1);
	double delta = (upper - lower) / (n - 1);
	result(1) = lower;
	result(n) = upper;
	for (int i = 2; i <= n - 1; i++)
	{
		result(i) = result(i - 1) + delta;
	}
	return result;
}

void save(const char *fileName, Matrix A)
{

	int col = 20;
	int m = A.getM(), n = A.getN();
	int size = m * n;
	std::ofstream out(fileName);
	out << m << " " << n << endl;
	vector<double> vdbl = A.getVdbl();
	col = (n < col) ? n : col;
	for (int i = 0; i <= size - 1; i++)
	{
		if (i % col == 0)
			out << endl;
		out << vdbl[i] << " ";
	}
}

Matrix load(const char *fileName)
{

	int m, n;
	std::ifstream is(fileName);
	is >> m >> n;
	vector<double> vdbl;
	double element;
	for (int i = 1; i <= m * n; i++)
	{
		is >> element;
		vdbl.push_back(element);
	}
	return Matrix(m, n, vdbl);
}

Matrix normalize(Matrix A)
{
	int m = A.getM(), n = A.getN();
	Matrix result(m, n);
	for (int j = 1; j <= n; j++)
	{
		double size = vNorm2(A(1, m, j, j));
		if (size != 0)
		{
			result.setblk(1, j, (1 / size) * A(1, m, j, j));
		}
	}
	return result;
}

Matrix combine(Matrix A, Matrix B)
{
	int m1 = A.getM();
	int n1 = A.getN();
	int m2 = B.getM();
	int n2 = B.getN();
	Matrix C(m1, n1 + n2);
	if (m1 != m2)
		return C;
	C.setblk(1, 1, A);
	C.setblk(1, n1 + 1, B);
	return C;
}
