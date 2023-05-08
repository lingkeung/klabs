#include "graph.h"
#include "matrix.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphPolar.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "stats.h"

using namespace std;

Graph::Graph(Matrix d)
{
	nx = d.getM() * d.getN();
	x = linspace(nx, 1, nx);
	y = d;
}

Graph::Graph(int nx, Matrix x, Matrix y)
{
	this->nx = nx;
	this->x = x;
	this->y = y;
}

Graph::Graph(int nx, Matrix x, Matrix y, Matrix z)
{

	this->nx = nx;
	this->x = x;
	this->y = y;
	this->z = z;
}

Graph::Graph(int nx, int nu, Matrix x, Matrix y, Matrix u, Matrix v)
{

	this->nx = nx;
	this->x = x;
	this->y = y;
	this->nu = nu;
	this->x = x;
	this->y = y;
	this->u = u;
	this->v = v;
}

Graph::Graph(int nx, int nu, Matrix x, Matrix y, Matrix z, Matrix u, Matrix v, Matrix w)
{

	this->nx = nx;
	this->x = x;
	this->y = y;
	this->nu = nu;
	this->x = x;
	this->y = y;
	this->z = z;
	this->u = u;
	this->v = v;
	this->w = w;
}

void Graph::xy()
{
	double ax[nx], ay[nx];
	for (int i = 0; i < nx; i++)
	{
		ax[i] = x(i + 1);
		ay[i] = y(i + 1);
	}
	auto c = new TCanvas();
	auto g = new TGraph(nx, ax, ay);
	g->SetTitle(title);
	g->SetMarkerColor(kBlue);
	g->SetLineColor(kBlue);
	g->SetMarkerStyle(marker);
	g->GetXaxis()->SetTitleOffset(1.5);
	g->GetYaxis()->SetTitleOffset(1.5);
	g->Draw(opx);
}

void Graph::yz()
{
	double ay[nx], az[nx];
	for (int i = 0; i < nx; i++)
	{
		ay[i] = y(i + 1);
		az[i] = z(i + 1);
	}
	auto c = new TCanvas();
	auto g = new TGraph(nx, ay, az);
	g->SetTitle(title);
	g->SetMarkerColor(kBlue);
	g->SetLineColor(kBlue);
	g->SetMarkerStyle(marker);
	g->GetXaxis()->SetTitleOffset(1.5);
	g->GetYaxis()->SetTitleOffset(1.5);
	g->Draw(opx);
}

void Graph::zx()
{
	double ax[nx], az[nx];
	for (int i = 0; i < nx; i++)
	{
		ax[i] = x(i + 1);
		az[i] = z(i + 1);
	}
	auto c = new TCanvas();
	auto g = new TGraph(nx, az, ax);
	g->SetTitle(title);
	g->SetMarkerColor(kBlue);
	g->SetLineColor(kBlue);
	g->SetMarkerStyle(marker);
	g->GetXaxis()->SetTitleOffset(1.5);
	g->GetYaxis()->SetTitleOffset(1.5);
	g->Draw(opx);
}

void Graph::xyuv()
{
	double ax[nx], ay[nx], au[nu], av[nu];
	for (int i = 0; i < nx; i++)
	{
		ax[i] = x(i + 1);
		ay[i] = y(i + 1);
	}
	for (int i = 0; i < nu; i++)
	{
		au[i] = u(i + 1);
		av[i] = v(i + 1);
	}
	auto c = new TCanvas();
	auto g = new TGraph(nx, ax, ay);
	g->SetTitle(title);
	g->SetMarkerColor(kBlue);
	g->SetLineColor(kBlue);
	g->SetMarkerStyle(marker);
	g->GetXaxis()->SetTitleOffset(1.5);
	g->GetYaxis()->SetTitleOffset(1.5);
	g->Draw(opx);
	auto g2 = new TGraph(nu, au, av);
	g2->Draw(opu);
}

void Graph::xyz()
{
	double ax[nx], ay[nx], az[nx];
	for (int i = 0; i < nx; i++)
	{
		ax[i] = x(i + 1);
		ay[i] = y(i + 1);
		az[i] = z(i + 1);
	}
	auto c = new TCanvas();
	auto dt = new TGraph2D();
	for (int i = 0; i < nx; i++)
	{
		dt->SetPoint(i, ax[i], ay[i], az[i]);
	}
	dt->SetTitle(title);
	dt->SetMarkerColor(kBlue);
	dt->SetLineColor(kBlue);
	dt->SetMarkerStyle(marker);
	dt->GetXaxis()->SetTitleOffset(1.5);
	dt->GetYaxis()->SetTitleOffset(1.5);
	dt->GetZaxis()->SetTitleOffset(1.5);
	dt->Draw(opx);
}

void Graph::xyzuvw()
{
	double ax[nx], ay[nx], az[nx], au[nu], av[nu], aw[nu];
	for (int i = 0; i < nx; i++)
	{
		ax[i] = x(i + 1);
		ay[i] = y(i + 1);
		az[i] = z(i + 1);
		au[i] = u(i + 1);
		av[i] = v(i + 1);
		aw[i] = w(i + 1);
	}

	auto c = new TCanvas();

	auto dt = new TGraph2D();
	for (int i = 0; i < nx; i++)
	{
		dt->SetPoint(i, ax[i], ay[i], az[i]);
	}
	dt->SetTitle(title);
	dt->SetMarkerColor(kBlue);
	dt->SetLineColor(kBlue);
	dt->SetMarkerStyle(marker);
	dt->GetXaxis()->SetTitleOffset(1.5);
	dt->GetYaxis()->SetTitleOffset(1.5);
	dt->GetZaxis()->SetTitleOffset(1.5);
	dt->Draw(opx);

	auto dt2 = new TGraph2D();
	for (int i = 0; i < nu; i++)
	{
		dt2->SetPoint(i, au[i], av[i], aw[i]);
	}
	dt2->Draw(opu);
}

void Graph::explicitY(const char *name, const char *expr, double lower, double upper)
{
	auto c = new TCanvas();
	TF1 f1(name, expr, lower, upper);
	f1.DrawClone();
}

void Graph::implicitY(const char *name, const char *expr, double xlower, double xupper, double ylower, double yupper)
{
	auto c = new TCanvas();
	TF2 f2(name, expr, xlower, xupper, ylower, yupper);
	f2.DrawClone();
}

void Graph::polar()
{

	double ax[nx], ay[nx];
	for (int i = 0; i < nx; i++)
	{
		ax[i] = x(i + 1);
		ay[i] = y(i + 1);
	}
	auto c = new TCanvas();
	auto g = new TGraphPolar(nx, ax, ay);
	g->SetTitle(title);
	g->SetLineWidth(3);
	g->SetLineColor(2);
	g->Draw("L");
}

void Graph::d()
{
	this->xy();
}

void Graph::hd(int bins, const char *title)
{
	double min = minimum(y);
	double max = maximum(y);
	auto h = new TH1D("h", title, bins, min, max);
	int size = y.getM() * y.getN();
	for (int i = 1; i <= size; i++)
	{
		h->Fill(y(i));
	}
	h->Draw();
}

void histogram(Matrix y, int bins, double min, double max, const char *title)
{
	auto h = new TH1D("h", title, bins, min, max);
	int size = y.getM() * y.getN();
	for (int i = 1; i <= size; i++)
	{
		h->Fill(y(i));
	}
	h->Draw();
}

void histogram(int size, int arr[], int bins, double min, double max, const char *title)
{
	auto h = new TH1D("h", title, bins, min, max);
	for (int i = 0; i < size; i++)
	{
		h->Fill(arr[i]);
	}
	h->Draw();
}
