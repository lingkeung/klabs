#ifndef CN_H
#define CN_H

class C
{
public:
	double a;
	double b;

public:
	C() : a{0}, b{0} {}
	C(double a, double b) : a{a}, b{b} {}
	void print(int precision = 2);
	friend C operator+(const C &w, const C &z);
	friend C operator-(const C &w, const C &z);
	friend C operator*(const C &w, const C &z);
	friend C operator*(const double &s, const C &z);
	friend C operator/(const C &w, const C &z);
	friend C operator~(const C &z);
	friend bool operator==(const C &w, const C &z);
};

double abs(C cn);
bool isReal(C cn);
C p2r(double amplitude, double theta);
C pp2r(double amplitude, double theta, double power);

#endif