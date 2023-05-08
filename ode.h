#ifndef ODE_H
#define ODE_H

#include <vector>

using namespace std;

class ODE
{
	double t0;
	double t1;
	int npt;
	vector<double> X0;
	vector<double (*)(double, double)> vf1;
	vector<double (*)(double, double, double)> vf2;
	vector<double (*)(double, double, double, double)> vf3;
	vector<double (*)(double, double, double, double, double)> vf4;
	vector<double (*)(double, double, double, double, double, double, double, double, double)> vf8;
	vector<double (*)(double, double, double, double, double, double, double, double, double,
					  double, double, double, double)>
		vf12;

public:
	vector<vector<double>> tX;

public:
	ODE(vector<double> X, vector<double (*)(double, double)> vf)
	{
		X0 = X;
		vf1 = vf;
	}
	ODE(vector<double> X, vector<double (*)(double, double, double)> vf)
	{
		X0 = X;
		vf2 = vf;
	}
	ODE(vector<double> X, vector<double (*)(double, double, double, double)> vf)
	{
		X0 = X;
		vf3 = vf;
	}
	ODE(vector<double> X, vector<double (*)(double, double, double, double, double)> vf)
	{
		X0 = X;
		vf4 = vf;
	}
	ODE(vector<double> X, vector<double (*)(double, double, double, double, double,
											double, double, double, double)>
							  vf)
	{
		X0 = X;
		vf8 = vf;
	}
	ODE(vector<double> X, vector<double (*)(double, double, double, double, double,
											double, double, double, double,
											double, double, double, double)>
							  vf)
	{
		X0 = X;
		vf12 = vf;
	}
	void rk41();
	void rk42();
	void rk43();
	void rk44();
	void rk48();
	void rk412();

	void set_t0(double t) { t0 = t; }
	void set_t1(double t) { t1 = t; }
	void set_npt(int n) { npt = n; }

	const char *op;
	const char *title;

	void ptx(int idx, const char *title = "title", const char *op = "AP*");

	void pxy(int idx, int idy, const char *title = "title", const char *op = "AP*");

	void pxyz(int idx, int idy, int idz, const char *title = "title", const char *op = "P0");
};

#endif
