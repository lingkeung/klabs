#include "newcmatrix.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
// #include "matrix.h"

using namespace std;
using namespace kling;

typedef complex<double> C;

cMatrix::cMatrix(const int m, const int n) : m{m}, n{n}
{
    if (m <= 0 || n <= 0)
        return;
    vC = vector<C>(m * n, C());
}

cMatrix::cMatrix(const int m, const int n, vector<C> vC)
{
    if (m * n != vC.size())
        return;
    this->m = m;
    this->n = n;
    this->vC = vC;
}

void cMatrix::print()
{

    for (auto i = 1; i <= m; i++)
    {

        for (auto j = 1; j <= n; j++)
        {
            cout << setw(15);
            cout << (*this)(i, j);
        }
        cout << endl;
    }
}

cMatrix cMatrix::hermitian()
{
    cMatrix B(this->n, this->m);
    for (int i = 1; i <= this->n; i++)
        for (int j = 1; j <= this->m; j++)
            B(i, j) = conj((*this)(j, i));
    return B;
}

C cMatrix::operator()(const int i, const int j) const
{
    int k = (i - 1) * n + (j - 1);
    return vC[k];
}

C &cMatrix::operator()(const int i, const int j)
{
    int k = (i - 1) * n + (j - 1);
    return vC[k];
}

C cMatrix::operator()(const int k) const
{
    return vC[k - 1];
}

C &cMatrix::operator()(const int k)
{
    return vC[k - 1];
}

cMatrix cMatrix::getblk(int i1, int i2, int j1, int j2)
{
    int m = abs(i2 - i1 + 1);
    int n = abs(j2 - j1 + 1);
    cMatrix result(m, n);
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

void cMatrix::setblk(int r, int c, cMatrix A)
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

cMatrix cMatrix::operator()(int i1, int i2, int j1, int j2)
{
    return getblk(i1, i2, j1, j2);
}

cMatrix kling::operator+(const cMatrix &A, const cMatrix &B)
{
    cMatrix cC(A.m, A.n);
    if (A.m != B.m || A.n != B.n)
        return cC;
    for (int i = 1; i <= A.m; i++)
        for (int j = 1; j <= A.n; j++)
            cC(i, j) = A(i, j) + B(i, j);
    return cC;
}

cMatrix kling::operator-(const cMatrix &A, const cMatrix &B)
{
    cMatrix cC(A.m, A.n);
    if (A.m != B.m || A.n != B.n)
        return cC;
    for (int i = 1; i <= A.m; i++)
        for (int j = 1; j <= A.n; j++)
            cC(i, j) = A(i, j) - B(i, j);
    return cC;
}

cMatrix kling::operator*(const C &scalar, const cMatrix &A)
{
    cMatrix B(A.m, A.n);
    for (int i = 1; i <= A.m; i++)
        for (int j = 1; j <= A.n; j++)
            B(i, j) = scalar * A(i, j);
    return B;
}

cMatrix kling::operator*(const double &scalar, const cMatrix &A)
{
    cMatrix B(A.m, A.n);
    for (int i = 1; i <= A.m; i++)
        for (int j = 1; j <= A.n; j++)
            B(i, j) = scalar * A(i, j);
    return B;
}

cMatrix kling::operator*(const cMatrix &A, const cMatrix &B)
{
    cMatrix cC(A.m, B.n);
    if (A.n != B.m)
        return cC;
    for (int i = 1; i <= A.m; i++)
    {
        for (int j = 1; j <= B.n; j++)
        {
            C sum = C(0, 0);
            for (int k = 1; k <= A.n; k++)
            {
                sum = sum + (A(i, k) * B(k, j));
            }
            cC(i, j) = sum;
        }
    }
    return cC;
}

/*C p2r(double amplitude, double theta)
{
    C result(amplitude * cos(theta), amplitude * sin(theta));
    return result;
}

C newpp2r(double amplitude, double theta, double power)
{
        double amp = pow(amplitude, power);
        double arg = theta * power;
        return polar(amp, arg);
}
*/
cMatrix kling::combine(cMatrix A, cMatrix B)
{
    int m1 = A.getM();
    int n1 = A.getN();
    int m2 = B.getM();
    int n2 = B.getN();
    cMatrix C(m1, n1 + n2);
    if (m1 != m2)
        return C;
    C.setblk(1, 1, A);
    C.setblk(1, n1 + 1, B);
    return C;
}

cMatrix kling::cIdentity(int order)
{
    cMatrix A(order, order);
    for (auto i = 1; i <= order; i++)
    {
        A(i, i) = C(1.0, 0);
    }
    return A;
}

cMatrix kling::cplex(Matrix A)
{
    int m = A.getM();
    int n = A.getN();
    cMatrix result(m, n);
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            result(i, j) = C(A(i, j), 0);
        }
    }
    return result;
}

double kling::vNorm2(cMatrix x)
{
    double result = 0;
    int m = x.getM();
    int n = x.getN();
    if (m != 1 && n != 1)
    {
        return 1;
    }
    int size = (m > n) ? m : n;
    for (int i = 1; i <= size; i++)
    {
        result += x(i).real() * x(i).real() + x(i).imag() * x(i).imag();
    }
    return sqrt(result);
}

void kling::vNormal(cMatrix &v)
{
    v = (1.0 / vNorm2(v)) * v;
}

double kling::vNormInf(cMatrix x)
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

double kling::normInf(cMatrix A)
{
    double result = 0;
    for (int i = 1; i <= A.getM(); i++)
    {
        for (int j = 1; j <= A.getN(); j++)
        {
            if (abs(A(i, j)) >= result)
            {
                result = abs(A(i, j));
            }
        }
    }
    return result;
}

void kling::vSort(cMatrix &v)
{
    int m = v.getM();
    int n = v.getN();
    if (m == 1 || n == 1)
    {
        vector<C> vC = v.getvC();
        sort(vC.begin(), vC.end(), [](C a, C b)
             { return abs(a) > abs(b); });
        v = cMatrix(m, n, vC);
    }
}

cMatrix kling::distinct(cMatrix v)
{
    vector<C> vec;
    vSort(v);
    int m = v.getM();
    int n = v.getN();
    if ((m == 1 && n > 1) || (n == 1 && m > 1))
    {
        int size = (m > n) ? m : n;
        for (int i = 1; i <= size; i++)
        {
            int k = 0;
            for (int j = i + 1; j <= size; j++)
            {
                if (v(i) == v(j))
                {
                    k++;
                }
            }
            vec.push_back(v(i));
            i += k;
        }
        return cMatrix(vec.size(), 1, vec);
    }
    else
        return v;
}

Matrix kling::real(cMatrix A)
{
    int m = A.getM();
    int n = A.getN();
    Matrix result(m, n);
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            result(i, j) = A(i, j).real();
        }
    }
    return result;
}

bool kling::isReal(cMatrix A)
{
    bool result = true;
    for (int i = 1; i <= A.getM(); i++)
    {
        for (int j = 1; j <= A.getN(); j++)
        {
            if (!(abs(A(i, j).imag()) <= 1e-6))
            {
                result = false;
                return result;
            }
        }
    }
    return result;
}
/*
cMatrix ndft(Matrix x) // naive discrete fourier transform
{
    const double pi = M_PI;
    int N = x.getM();
    double ampWN = 1;
    double argWN = -2 * pi / N; // WN = exp(-i*2*pi/N)
    cMatrix W(N, N);
    for (int n = 0; n <= N - 1; n++)
    {
        for (int k = 0; k <= N - 1; k++)
        {
            W(n + 1, k + 1) = pp2r(ampWN, argWN, n * k); // powered polar to rectangular
        }
    }
    cMatrix cx = cplex(x);
    cMatrix X = W * cx;
    // cMatrix XoN = (1.0 / N) * X;
    // return XoN;
    return X;
}

Matrix nidft(cMatrix X) // naive inverse discrete fourier transform
{
    const double pi = M_PI;
    int N = X.getM();
    double ampWN = 1;
    double argWN = 2 * pi / N; // WN = exp(i*2*pi/N)
    cMatrix W(N, N);
    for (int n = 0; n <= N - 1; n++)
    {
        for (int k = 0; k <= N - 1; k++)
        {
            W(n + 1, k + 1) = pp2r(ampWN, argWN, n * k); // powered polar to rectangular
        }
    }
    cMatrix x = W * X;
    cMatrix xoN = (1.0 / N) * x;
    return real(xoN);
}
*/
cMatrix kling::fdft(Matrix x) // fast discrete fourier transform, N = 2^n
{
    int N = x.getM();
    if (N == 1) // base case
    {
        return kling::cplex(x);
    }
    else
    {
        int m = N / 2; // step 1 split
        Matrix x1(m, 1);
        Matrix x2(m, 1);
        for (int i = 1; i <= m; i++)
        {
            x1(i) = x(2 * i - 1);
            x2(i) = x(2 * i);
        }

        cMatrix X1 = fdft(x1); // step 2 recursive calls
        cMatrix X2 = fdft(x2);

        cMatrix X(N, 1); // step 3 combine
        const double pi = M_PI;
        double ampWN = 1;
        double argWN = -2 * pi / N;
        for (int j = 0; j <= m - 1; j++)
        {
            C WNj = pow(polar(ampWN, argWN), j);
            X(j + 1) = X1(j + 1) + WNj * X2(j + 1);     // these formulas are the heart of fft!
            X(j + 1 + m) = X1(j + 1) - WNj * X2(j + 1); // this is Cooley and Tukey's contribution.
        }
        // X.print();
        return X;
    }
}

cMatrix kling::fidftHelp(cMatrix x) // helper fast inverse discrete fourier transform, N = 2^n
{
    int N = x.getM();
    if (N == 1) // base case
    {
        return (x);
    }
    else
    {
        int m = N / 2; // step 1 split
        cMatrix x1(m, 1);
        cMatrix x2(m, 1);
        for (int i = 1; i <= m; i++)
        {
            x1(i) = x(2 * i - 1);
            x2(i) = x(2 * i);
        }

        cMatrix X1 = fidftHelp(x1); // step 2 recursive calls
        cMatrix X2 = fidftHelp(x2);

        cMatrix X(N, 1); // step 3 combine
        const double pi = M_PI;
        double ampWN = 1;
        double argWN = 2 * pi / N;
        for (int j = 0; j <= m - 1; j++)
        {
            // C WNj = pp2r(ampWN, argWN, j);
            C WNj = pow(polar(ampWN, argWN), j);
            X(j + 1) = X1(j + 1) + WNj * X2(j + 1);     // these formulas are the heart of fft!
            X(j + 1 + m) = X1(j + 1) - WNj * X2(j + 1); // this is Cooley and Tukey's contribution.
        }
        // X.print();
        return X;
    }
}

Matrix kling::fidft(cMatrix X) // fast inverse discrete fourier transform
{
    int N = X.getM();
    return real((1.0 / N) * fidftHelp(X));
}
