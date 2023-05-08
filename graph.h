#ifndef GRAPH_H
#define GRAPH_H
#include "matrix.h"

using namespace std;

class Graph
{

private:
	int nx;
	Matrix x;
	Matrix y;
	Matrix z;
	const char *opx = "AP";

	int nu;
	Matrix u;
	Matrix v;
	Matrix w;
	const char *opu = "AP";

	const char *title = "Graph";
	int marker = 4;

public:
	Graph(){};
	Graph(Matrix d);
	Graph(int nx, Matrix x, Matrix y);
	Graph(int nx, Matrix x, Matrix y, Matrix z);
	Graph(int nx, int nu, Matrix x, Matrix y, Matrix u, Matrix v);
	Graph(int nx, int nu, Matrix x, Matrix y, Matrix z, Matrix u, Matrix v, Matrix w);

	void set_opx(const char *op) { opx = op; }
	void set_opu(const char *op) { opu = op; }
	void set_title(const char *str) { title = str; }
	void set_marker(int mark) { marker = mark; }

	int get_nx() { return nx; }
	Matrix get_x() { return x; }
	Matrix get_y() { return y; }
	Matrix get_z() { return z; }

	int get_nu() { return nu; }
	Matrix get_u() { return u; }
	Matrix get_v() { return v; }
	Matrix get_w() { return w; }

	void xy();
	void yz();
	void zx();
	void xyz();
	void xyuv();
	void xyzuvw();
	void explicitY(const char *name, const char *expr, double lower, double upper);
	void implicitY(const char *name, const char *expr, double xlower, double xupper, double ylower, double yupper);
	void polar();
	void d();
	void hd(int bins, const char *title = "Histogram-y"); // draw histogram of y
};

void histogram(Matrix y, int bins, double min, double max, const char *title);
void histogram(int size, int arr[], int bins, double min, double max, const char *title);

#endif