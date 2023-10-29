#include "cn.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

void C::print(int precision)
{
	cout << fixed << setprecision(precision);
	if (b >= 0)
		cout << "(" << a << " + " << b << "i)" << endl;
	else
		cout << "(" << a << " - " << abs(b) << "i)" << endl;
}

C operator+(const C &w, const C &z)
{
	return C(w.a + z.a, w.b + z.b);
}

C operator-(const C &w, const C &z)
{
	return C(w.a - z.a, w.b - z.b);
}

C operator*(const C &w, const C &z)
{
	return C(w.a * z.a - w.b * z.b, w.b * z.a + w.a * z.b);
}

C operator*(const double &s, const C &z)
{
	return C(s * z.a, s * z.b);
}

C operator/(const C &w, const C &z)
{
	return (1 / (z.a * z.a + z.b * z.b)) * C(w.a * z.a + w.b * z.b, w.b * z.a - w.a * z.b);
}

C operator~(const C &z)
{
	return C(z.a, -1 * z.b);
}

double abs(C cn)
{
	double a = cn.a;
	double b = cn.b;
	double r = sqrt(a * a + b * b);
	return r;
}

bool operator==(const C &w, const C &z)
{
	double tol = 1e-8;
	double dbl1 = abs(w.a - z.a);
	double dbl2 = abs(w.b - z.b);
	bool result = (dbl1 < tol) && (dbl2 < tol);
	return result;
}

bool isReal(C cn)
{
	return (abs(cn.b) < 1e-8) ? true : false;
}

C p2r(double amplitude, double theta)
{
	C result(amplitude * cos(theta), amplitude * sin(theta));
	return result;
}

C pp2r(double amplitude, double theta, double power)
{
		double amp = pow(amplitude, power);
		double arg = theta * power;
		return p2r(amp, arg);
}