#include <iostream>
#include <cmath>
#include <vector>
#include "graph.h"
#include "TApplication.h"
#include "eqn.h"
#include "ode.h"

using namespace std;

double f(double x)
{
	return (x * x - 2);
}

double df(double x)
{
	return 2 * x;
}

double dx1(double t, double x1, double x2)
{
	return x2;
}

double dx2(double t, double x1, double x2)
{
	return -2 * x1 + 2 * x2 + exp(2 * t);
}

double dx21(double t, double x1, double x2)
{
	return -2 * x1 - 2 * x2;
}

double ldx(double t, double x, double y, double z)
{
	return 10 * (y - x);
}
double ldy(double t, double x, double y, double z)
{
	return 28 * x - y - x * z;
}
double ldz(double t, double x, double y, double z)
{
	return x * y - 8.0 / 3.0 * z;
}

double pdy1(double t, double y1, double y2)
{
	return y2;
}

double pdy2(double t, double y1, double y2)
{
	return -9.81 / 1 * sin(y1);
}
double dx(double t, double x)
{
	return t * t * t - x / t;
}
double odx(double t, double x, double vx, double y, double vy)
{
	return vx;
}
double odvx(double t, double x, double vx, double y, double vy)
{
	return -1 * 3 * x / pow((x * x + y * y), 3.0 / 2);
}
double ody(double t, double x, double vx, double y, double vy)
{
	return vy;
}
double odvy(double t, double x, double vx, double y, double vy)
{
	return -1 * 3 * y / pow((x * x + y * y), 3.0 / 2);
}

void StandaloneApplication(int argc, char **argv)
{

	/*vector<double(*)(double,double)> vf0{dx};
	vector<double>X00{2.0/5};
	ODE ode(X00,vf);
	ode.set_t0(1);
	ode.set_t1(2);
	ode.set_npt(11);
	ode.rk41();
	for (int i = 0; i < 11; i++) {
		cout << ode.tX[1][i] << endl;
	}
	ode.ptx(1,"y;t;y");*/

	vector<double (*)(double, double, double)> vf;
	vf.push_back(dx1);
	vf.push_back(dx21);
	double x10 = 0, x20 = 1;
	vector<double> X0;
	X0.push_back(x10);
	X0.push_back(x20);

	ODE ode(X0, vf);
	ode.set_t0(0);
	ode.set_t1(10);
	ode.set_npt(101);

	ode.rk42();
	ode.ptx(1, "y;t;y");
	ode.ptx(2, "dy/dt;t;dy/dt");
	ode.pxy(1, 2, "Phase Diagram;y(t);dy(t)/dt");

	vector<double (*)(double, double, double, double)> vfl{ldx, ldy, ldz};
	vector<double> X0l{0, 1, 0};
	ODE l1(X0l, vfl);
	l1.set_t0(0);
	l1.set_t1(50);
	l1.set_npt(50 / 0.01 + 1);
	l1.rk43();
	l1.pxyz(1, 2, 3, "Lorenz;x;y;z", "P0");

	vector<double (*)(double, double, double)> vfp{pdy1, pdy2};
	vector<double> X0p{3.1416 / 2, 0};
	ODE p(X0p, vfp);
	p.set_t0(0);
	p.set_t1(10);
	p.set_npt(10 / 0.01 + 1);
	p.rk42();
	p.ptx(1, "Pendulum;t;y", "APC");
	p.ptx(2, "dy/dt;t;dy/dt", "APC");
	p.pxy(1, 2, "Phase Diagram;y;dy/dt", "APC");

	vector<double (*)(double, double, double, double, double)> vfo{odx, odvx, ody, odvy};
	vector<double> X0o{0, 1, 2, 0};
	ODE o(X0o, vfo);
	o.set_t0(0);
	o.set_t1(100);
	o.set_npt(10001);
	o.rk44();
	o.ptx(4, "Orbit4;x;y;", "APC");
	o.ptx(3, "Orbit3;x;y;", "APC");
	o.ptx(2, "Orbit2;x;y;", "APC");
	o.ptx(1, "Orbit1;x;y;", "APC");
	o.pxy(1, 3, "Orbit;x;y;", "APC");
}
/*
int main(int argc, char** argv) {

   TApplication app("ROOT Application", &argc, argv);
   StandaloneApplication(app.Argc(), app.Argv());
   app.Run();
   return 0;

}*/
