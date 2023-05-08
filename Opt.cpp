#include <iostream>
#include <cmath>
#include "matrix.h"
#include <functional>
#include "plu.h"

using namespace std;

double fun(double x)
{
	return 3 * x * x * x - 4 * x + 2;
}

double fun2(double x)
{
	return (x - 10) * (x - 10);
}

void ab(std::function<double(double)> f, double x0, double h, double &a, double &b)
{

	double x1 = x0;
	double f1 = f(x1);
	double x2 = x1 + h;
	double f2 = f(x2);
	double x3 = 0;
	double f3 = 0;

	if (f1 > f2)
	{
		h = 2 * h;
	}
	else
	{
		h = -1 * h;
		x3 = x1;
		f3 = f1;
		x1 = x2;
		f1 = f2;
		x2 = x3;
		f2 = f3;
	}
	while (1)
	{
		x3 = x0 + h;
		f3 = f(x3);
		if (f2 < f3)
			break;
		else
		{
			h = 2 * h;
			x1 = x2;
			f1 = f2;
			x2 = x3;
			f2 = f3;
		}
	}
	if (h > 0)
	{
		a = x1;
		b = x3;
	}
	else
	{
		a = x3;
		b = x1;
	}
}

double golden(double f(double), double a, double b, double esp)
{

	double x1 = a + 0.382 * (b - a), f1 = f(x1);
	double x2 = a + 0.618 * (b - a), f2 = f(x2);
	while (1)
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + 0.382 * (b - a);
			f1 = f(x1);
		}
		else
		{
			a = x1;
			x1 = x2;
			f1 = f2;
			x2 = a + 0.618 * (b - a);
			f2 = f(x2);
		}
		if (abs(b - a) <= esp)
		{
			return (a + b) / 2;
		}
	}
}

double ming(double f(double), double x0, double h, double esp)
{

	double a, b;
	ab(f, x0, h, a, b);
	return golden(f, a, b, esp);
}

double quit(std::function<double(double)> f, double a, double b, double esp)
{

	double c = (a + b) / 2;
	double x1 = a, x2 = c, x3 = b;
	double f1 = f(x1), f2 = f(x2), f3 = f(x3);

	while (1)
	{
		double C1 = (x2 * x2 - x3 * x3) * f1 + (x3 * x3 - x1 * x1) * f2 + (x1 * x1 - x2 * x2) * f3;
		double C2 = (x2 - x3) * f1 + (x3 - x1) * f2 + (x1 - x2) * f3;
		double xp = 0.5 * C1 / C2;
		double fp = f(xp);

		if (abs(x2 - xp) <= esp && abs(f2 - fp) <= esp)
		{
			return (fp <= f2) ? xp : x2;
		}
		else
		{
			if (fp <= f2)
			{
				if (xp <= x2)
				{
					x3 = x2;
					x2 = xp;
					f3 = f2;
					f2 = fp;
				}
				else
				{
					x1 = x2;
					x2 = xp;
					f1 = f2;
					f2 = fp;
				}
			}
			else
			{
				if (xp <= x2)
				{
					x1 = xp;
					f1 = fp;
				}
				else
				{
					x3 = xp;
					f3 = fp;
				}
			}
		}
	}
}

double minq(std::function<double(double)> f, double x0, double h, double esp)
{
	double a, b;
	ab(f, x0, h, a, b);
	return quit(f, a, b, esp);
}

double f(Matrix x)
{
	return x(1) * x(1) + 2 * x(2) * x(2) - 2 * x(1) * x(2) - 4 * x(1);
}

Matrix del(Matrix x)
{
	int m = x.getM();
	Matrix out(m, 1);
	out(1) = 2 * x(1) - 2 * x(2) - 4;
	out(2) = -2 * x(1) + 4 * x(2);
	return out;
}

Matrix del2(Matrix x)
{
	int m = x.getM();
	Matrix out(m, m);
	out(1, 1) = 2;
	out(1, 2) = -2;
	out(2, 1) = -2;
	out(2, 2) = 4;
	return out;
}

int gradient(double f(Matrix), Matrix del(Matrix), Matrix X0, double esp, Matrix &Xmin, double &fmin)
{

	Matrix X = X0;
	int t = 0;
	for (t = 1; t <= 100; t++)
	{
		Matrix S = del(X);
		auto obj = [=](double a)
		{ return f(X + a * S); };
		double mina = minq(obj, 1, 1, esp);
		X = X + mina * S;
		if (vNormInf(del(X)) <= esp)
			break;
	}
	Xmin = X;
	fmin = f(X);
	return t;
}

int newton(double f(Matrix), Matrix del(Matrix), Matrix del2(Matrix), Matrix X0, double esp, Matrix &Xmin, double &fmin)
{

	Matrix X = X0;
	int t = 0;
	for (t = 1; t <= 100; t++)
	{
		Matrix A = del2(X);
		Matrix b = -1 * del(X);
		Matrix S = axb(A, b);
		auto obj = [=](double a)
		{ return f(X + a * S); };
		double mina = minq(obj, 1, 1, esp);
		X = X + mina * S;
		if (vNormInf(del(X)) <= esp)
			break;
	}
	Xmin = X;
	fmin = f(X);
	return t;
}

int congrad(double f(Matrix), Matrix del(Matrix), Matrix X0, double esp, Matrix &Xmin, double &fmin)
{

	int k = 0;
	int m = X0.getM();
	Matrix X1 = X0;
	Matrix S = -1 * del(X1);
	int t = 0;
	for (t = 1; t <= m * 100; t++)
	{
		auto obj = [=](double a)
		{ return f(X1 + a * S); };
		double mina = minq([=](double a)
						   { return f(X1 + a * S); },
						   1, 1, esp);
		Matrix X10 = X1;
		X1 = X1 + mina * S;
		if (vNormInf(del(X1)) <= esp)
		{
			break;
		}
		else if (k == m)
		{
			S = -1 * del(X1);
			continue;
		}
		else
		{
			double beta = vNorm2(del(X1)) / vNorm2(del(X10));
			beta = beta * beta;
			S = -1 * del(X1) + beta * S;
		}
		k++;
	}
	Xmin = X1;
	fmin = f(X1);
	return t;
}

int powell(double f(Matrix), Matrix X0, double esp, Matrix &Xmin, double &fmin)
{

	int m = X0.getM();
	Matrix SS = identity(m);
	vector<Matrix> X;
	X.push_back(X0);
	int t = 0;
	for (t = 1; t <= 100; t++)
	{
		for (int i = 1; i <= m; i++)
		{
			Matrix S = SS(1, m, i, i);
			auto obj = [=](double a)
			{ return f(X[i - 1] + a * S); };
			double mina = minq(obj, 1, 1, esp);
			X.push_back(X[i - 1] + mina * S);
		}
		if (vNormInf(X.back() - X[0]) <= esp && abs(f(X.back() - X[0])) <= esp)
		{
			break;
		}
		double dXmax = -1e10;
		int idx = 0;
		for (int j = 1; j <= m - 1; j++)
		{
			double dX = f(X[j - 1]) - f(X[j]);
			if (dX > dXmax)
			{
				dXmax = dX;
				idx = j;
			}
		}
		Matrix Send = X.back() - X[0];
		Matrix Xrefn = 2 * X.back() - X[0];
		double f1 = f(X[0]), f2 = f(X.back()), f3 = f(Xrefn);
		double y1 = dXmax * (f1 - f3) * (f1 - f3) / 2;
		double y2 = (f1 - 2 * f2 + f3) * (f1 - f2 - dXmax) * (f1 - f2 - dXmax);

		if (f3 < f1 && y2 < y1)
		{
			auto obj = [=](double a)
			{ return f(X.back() + a * Send); };
			double mina = minq(obj, 1, 1, 1e-6);
			Matrix Xend = (X.back() + mina * Send);

			for (int k = 1; k <= m; k++)
			{
				if (k >= idx && k < m - 1)
				{
					SS.setblk(1, k, SS(1, m, k + 1, k + 1));
				}
			}
			SS.setblk(1, m, Send);
			X.clear();
			X.push_back(Xend);
		}
		else
		{
			Matrix Xend = (f2 < f3) ? X.back() : Xrefn;
			X.clear();
			X.push_back(Xend);
		}
	}
	Xmin = X.back();
	fmin = f(X.back());
	return t;
}

Matrix mkh(int size, vector<double> data)
{
	Matrix H(size, size);
	if (data.size() == size * (size + 1) / 2)
	{
		int k = 0;
		for (int i = 1; i <= size; i++)
		{
			for (int j = i; j <= size; j++)
			{
				H(j, i) = H(i, j) = data[k];
				k++;
			}
		}
	}
	return H;
}
int main()
{

	Matrix H1 = mkh(3, vector<double>{180, 72, 220, 120, -60, 140});
	Matrix C(3, 1);
	Matrix Ain(1, 3, vector<double>{-0.92, -0.64, -0.41});
	Matrix Aeq(1, 3, vector<double>{1, 1, 1});
	Matrix Bin(1, 1, vector<double>{-0.65});
	Matrix Beq(1, 1, vector<double>{1});

	/*
		Matrix Xmin;
		double fmin;
		Matrix X0(2,1,vector<double>{1,1});
		double esp = 1e-6;
		int t = gradient(f,del,X0,esp,Xmin,fmin);
		Xmin.print();
		cout << "\nfmin = " << fmin << endl;
		cout << "after " << t << " iterations" << endl;
		t = newton(f,del,del2,X0,esp,Xmin,fmin);
		Xmin.print();
		cout << "\nfmin = " << fmin << endl;
		cout << "after " << t << " iterations" << endl;
		t = congrad(f,del,X0,esp,Xmin,fmin);
		Xmin.print();
		cout << "\nfmin = " << fmin << endl;
		cout << "after " << t << " m-cycles" << endl;
		t = powell(f,X0,esp,Xmin,fmin);
		Xmin.print();
		cout << "\nfmin = " << fmin << endl;
		cout << "after " << t << " iterations" << endl;
	*/

	cout << "This is typed with my new Acer mechanical keyboard" << endl;
	return 0;
}