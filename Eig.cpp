#include "eig.h"
#include "qr.h"
#include <vector>
#include <iostream> // debug

Matrix qrSchur(Matrix A, int iter)
{
	Matrix A0 = A;
	for (int i = 1; i <= iter; i++)
	{
		Matrix Q, R, Q1, R1;
		qr(A0, Q, R, Q1, R1);
		Matrix A1 = R * Q;
		A0 = A1;
	}
	return A0;
}

void schur2eig(Matrix Schur, cMatrix &eig)
{
	Matrix S = Schur;
	int n = S.getN();
	vector<double> re, im;
	int k = 0;
	double a11, a12, a21, a22, a, b, c, disc;
	for (k = 1; k <= n - 1; k++)
	{
		a11 = S(k, k);
		a12 = S(k, k + 1);
		a21 = S(k + 1, k);
		a22 = S(k + 1, k + 1);
		a = 1;
		b = -1 * (a11 + a22);
		c = (a11 * a22 - a12 * a21);
		disc = b * b - 4 * a * c;
		if (abs(disc) <= 1e-10)
			disc = 0; // critical to avoid downstream problems *****
		if (abs(a21) <= 1e-6)
		{ // setting is critical**********
			// cout << "first if" << endl;
			re.push_back(a11);
			im.push_back(0);
		}
		else if (disc >= 0)
		{
			re.push_back((-1 * b + sqrt(disc)) / 2);
			im.push_back(0);
			re.push_back((-1 * b - sqrt(disc)) / 2);
			im.push_back(0);
			k++;
		}
		else if (disc < 0)
		{
			re.push_back((-1 * b) / 2);
			im.push_back(sqrt(-1 * disc) / 2);
			re.push_back((-1 * b) / 2);
			im.push_back(-1 * sqrt(-1 * disc) / 2);
			k++;
		}
	}
	if (k == n)
	{
		re.push_back(S(n, n));
		im.push_back(0);
	}
	eig = cMatrix(n, 1);
	for (int i = 1; i <= n; i++)
	{
		eig(i) = C(re[i - 1], im[i - 1]);
	}
}

Matrix hessenberg(Matrix A)
{
	int n = A.getN();
	for (int k = 1; k <= n - 2; k++)
	{
		Matrix v(n - k, 1);
		double beta;
		house(A(k + 1, n, k, k), v, beta);
		A.setblk(k + 1, k, (identity(n - k) - beta * v * v.transpose()) * (A(k + 1, n, k, n)));
		A.setblk(1, k + 1, (A(1, n, k + 1, n) * (identity(n - k) - beta * v * v.transpose())));
	}
	return A;
}
void francis(Matrix &H)
{
	int n = H.getN();
	int m = n - 1;
	double s, t, x, y, z;
	s = H(m, m) + H(n, n);
	t = H(m, m) * H(n, n) - H(m, n) * H(n, m);
	x = H(1, 1) * H(1, 1) + H(1, 2) * H(2, 1) - s * H(1, 1) + t;
	y = H(2, 1) * (H(1, 1) + H(2, 2) - s);
	z = H(2, 1) * H(3, 2);
	if (abs(x) <= 1e-6)  // avoid x,y being too close or equal to zero which may result in 'cycling'
	{
		if (x >= 0)
		{
			x += 1e-3;  // value empirically adjusted 
		}
		else
		{
			x -= 1e-3;
		}
	}
	if (abs(y) <= 1e-6)
	{
		if (y >= 0)
		{
			y += 1e-3;
		}
		else
		{
			y -= 1e-3;
		}
	}
	for (int k = 0; k <= n - 3; k++)
	{
		Matrix xyz(3, 1, vector<double>{x, y, z});
		Matrix v;
		double beta;
		house(xyz, v, beta);
		int q = (k >= 1) ? k : 1;		   
		H.setblk(k + 1, q, (identity(3) - beta * v * v.transpose()) * H(k + 1, k + 3, q, n));
		int r = (k + 4 <= n) ? k + 4 : n;
		H.setblk(1, k + 1, H(1, r, k + 1, k + 3) * (identity(3) - beta * v * v.transpose()));
		x = H(k + 2, k + 1);
		y = H(k + 3, k + 1);
		if (k < n - 3)
		{
			z = H(k + 4, k + 1);
		}
	}

	Matrix xy(2, 1, vector<double>{x, y});
	Matrix v;
	double beta;
	house(xy, v, beta);
	H.setblk(n - 1, n - 2, (identity(2) - beta * v * v.transpose()) * H(n - 1, n, n - 2, n));
	H.setblk(1, n - 1, H(1, n, n - 1, n) * (identity(2) - beta * v * v.transpose()));
	// return H;
}

void lmn(Matrix H, int &i1, int &i2, int &j1, int &j2, int &m)
{
	double tol = 1e-10;
	int r = H.getM();
	m = 0;
	int check = 0;
	int row = r + 1; // ensures initially m = r - row + 1 = 0
	for (int i = r; i >= 2; i--)
	{
		if (abs(H(i, i - 1)) <= tol)
		{
			row = i;
			check = 0;
		}
		else
		{
			check++;
			if (check >= 2)
				break;
		}
	}
	m = r - row + 1;

	if (m == r - 1)
	{
		m = r;
	}

	int l = 0;
	int lcheck = 0;
	int lrow = 1; // ensures initially l = lrow - 1 = 0
	for (int i = 2; i < row; i++)
	{ // this for block ensures maximum m and minimum l
		if (abs(H(i, i - 1)) <= tol)
		{
			lrow = i;
			break;
		}
	}
	l = lrow - 1;
	int n = r - l - m; // size of H22
	if (n <= 2)
		m = r; // among other effects this handles case of R11 is 2x2 block
	i1 = l + 1;
	i2 = i1 + n - 1;
	j1 = i1;
	j2 = j1 + n - 1;
}

void deflate(Matrix &H, double tol)
{
	int r = H.getM();
	for (int i = 2; i <= r; i++)
	{ // forces small sub-diagonal elements to zero
		if (abs(H(i, i - 1)) <= tol * (abs(H(i, i)) + abs(H(i - 1, i - 1))))
		{
			H(i, i - 1) = 0;
		}
	}
	for (int j = 1; j <= r; j++)
	{ // maintains original zero elements of Hessenberg zero
		for (int i = j + 2; i <= r; i++)
		{
			H(i, j) = 0;
		}
	}
}

void keig(Matrix A, double tol, cMatrix &eig, int iter)
{
	int r = A.getM();
	if (r == 1)
	{
		eig = cMatrix(1, 1, vector<C>{C(A(1, 1), 0)});
	}
	else if (r == 2)
	{
		schur2eig(A, eig);
	}
	else
	{
		Matrix H = hessenberg(A);
		for (int i = 1; i <= iter; i++)
		{
			deflate(H, tol);
			int i1, i2, j1, j2, m;
			lmn(H, i1, i2, j1, j2, m);
			if (m == r)
			{
				break;
			}
			Matrix H22 = H(i1, i2, j1, j2);
			francis(H22);
			H.setblk(i1, j1, H22);
		}

		schur2eig(H, eig);
	}
}

cMatrix eigVal2Vec(Matrix A, C eigval)
{
	int n = A.getN();
	double offset = -1.0001;
	cMatrix cA = cplex(A) + offset * eigval * cIdentity(n);
	cMatrix P, L, U;
	C det;
	cPlu(cA, P, L, U, det);
	cMatrix e(n, 1, vector<C>(n, C(1, 0)));
	cMatrix v1 = bsub(U, e);
	cMatrix v2;
	double tol = 1e-3;
	int trials = 0;
	int maxt = 100;
	for (trials = 0; trials < maxt; trials++)
	{
		v2 = cAxb(cA, v1);
		v2 = (1 / vNorm2(v2)) * v2;
		v1 = v2;
		cMatrix Ax = cplex(A) * v1;
		cMatrix lamdax = eigval * v1;
		if (vNormInf(Ax - lamdax) < tol)
			break;
	}
	return v1;
}

cMatrix eigVal2VecN(Matrix A, C eigval)
{
	int n = A.getN();
	cMatrix cA = cplex(A);
	cA = cA - (eigval * cIdentity(n));
	return cNull(cA);
}