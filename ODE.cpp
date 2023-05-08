#include "ode.h"
#include "matrix.h"
#include "graph.h"
#include "eqn.h" // rk4 solvers

using namespace std;

void ODE::rk41()
{
	tX = rk4(t0, X0[0], t1, npt, vf1[0]);
}
void ODE::rk42()
{
	tX = rk4(t0, X0[0], X0[1], t1, npt, vf2);
}
void ODE::rk43()
{
	tX = rk4(t0, X0[0], X0[1], X0[2], t1, npt, vf3);
}
void ODE::rk44()
{
	tX = rk4(t0, X0, t1, npt, vf4);
}
void ODE::rk48()
{
	tX = rk4(t0, X0, t1, npt, vf8);
}
void ODE::rk412()
{
	tX = rk4(t0, X0, t1, npt, vf12);
}

void ODE::ptx(int idx, const char *title, const char *op)
{
	int nva = tX.size() - 1;
	int npt = tX[0].size();
	Matrix t(npt, 1, tX[0]), x(npt, 1);
	for (int i = 1; i <= nva; i++)
	{
		if (idx == i)
		{
			x = Matrix(npt, 1, tX[i]);
			break;
		}
	}
	Graph g(npt, t, x);
	g.set_title(title);
	g.set_opx(op);
	g.xy();
}

void ODE::pxy(int idx, int idy, const char *title, const char *op)
{
	int nva = tX.size() - 1;
	int npt = tX[0].size();
	Matrix x(npt, 1), y(npt, 1);
	int count = 0;
	for (int i = 1; i <= nva; i++)
	{
		if (idx == i)
		{
			x = Matrix(npt, 1, tX[i]);
			count++;
		}
		if (idy == i)
		{
			y = Matrix(npt, 1, tX[i]);
			count++;
		}
		if (count == 2)
			break;
	}
	Graph g(npt, x, y);
	g.set_title(title);
	g.set_opx(op);
	g.xy();
}

void ODE::pxyz(int idx, int idy, int idz, const char *title, const char *op)
{
	int nva = tX.size() - 1;
	int npt = tX[0].size();
	Matrix x(npt, 1), y(npt, 1), z(npt, 1);
	int count = 0;
	for (int i = 1; i <= nva; i++)
	{
		if (idx == i)
		{
			x = Matrix(npt, 1, tX[i]);
			count++;
		}
		if (idy == i)
		{
			y = Matrix(npt, 1, tX[i]);
			count++;
		}
		if (idz == i)
		{
			z = Matrix(npt, 1, tX[i]);
			count++;
		}
		if (count == 3)
			break;
	}
	Graph g(npt, x, y, z);
	g.set_title(title);
	g.set_opx(op);
	g.xyz();
}
