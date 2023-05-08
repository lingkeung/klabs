#ifndef LR_H
#define LR_H

#include <vector>
#include "matrix.h"
#include "stats.h"
#include "prob.h"
#include <iostream>
#include <iomanip>
#include "data.h"
#include "qr.h"
#include "qrc.h"
#include <ctime>

using namespace std;

class LR
{

public:
	string fileName;
	vector<int> yx;
	Data d;
	int m, p;
	Matrix y;
	Matrix X;
	double sigma;
	Matrix beta;
	Matrix tStat;
	Matrix tPval;
	double R_squared;
	double adj_R_squared;
	double FStat;
	double FPval;

public:
	LR(string fileName, vector<int> yx);
	void print();
	double yf(vector<double> vxf);
	void yflh(double &lo, double &hi, vector<double> vxf, double percent);
	void yblh(double &lo, double &hi, vector<double> vxf, double percent);
	void print(vector<double> vxf);
};

LR::LR(string fileName, vector<int> yx)
{

	this->fileName = fileName;

	this->yx = yx;

	d = Data(fileName);

	m = d.D.getM();

	p = yx.size();

	y = d.D(1, m, yx[0], yx[0]);
	X = Matrix(m, p);
	X.setblk(1, 1, Matrix(m, 1, vector<double>(m, 1)));
	for (int i = 1; i < p; i++)
	{
		X.setblk(1, i + 1, d.D(1, m, yx[i], yx[i]));
	}

	QR q(X);
	beta = q.axb(y);

	Matrix var_beta;
	Matrix C = X.transpose() * X;
	QR c(C);
	C = c.inv();
	Matrix e = y - X * beta;
	double norme = vNorm2(e);
	double RSS = norme * norme;
	double sigmaSq = RSS / (m - p);
	vector<double> cjj;
	for (int j = 1; j <= p; j++)
	{
		cjj.push_back(C(j, j));
	}
	Matrix Cjj(p, 1, cjj);
	var_beta = sigmaSq * Cjj;

	sigma = sqrt(sigmaSq);
	tStat = Matrix(p, 1);
	for (int i = 1; i <= p; i++)
	{
		tStat(i) = beta(i) / sigma / sqrt(Cjj(i));
	}

	tPval = Matrix(p, 1);
	for (int i = 1; i <= p; i++)
	{
		tPval(i) = 2 * (1 - tdistribution_cdf(abs(tStat(i)), m - p));
	}

	double ybar = mean(y);
	double TSS = (y.transpose() * y)(1) - (m * ybar * ybar);

	R_squared = 1 - (RSS / TSS);

	adj_R_squared = 1 - (1 - R_squared) * (m - 1) / (m - p);

	FStat = (m - p) * R_squared / (p - 1) / (1 - R_squared);

	FPval = 1 - fdistribution_cdf(FStat, p - 1, m - p);
}

void LR::print()
{

	time_t now = time(NULL);

	int columnWidth1 = 15, w1 = columnWidth1;
	int columnWidth2 = 15, w2 = columnWidth2;
	int columnWidth3 = 15, w3 = columnWidth3;
	int width = 2 * w1 + w2 + w3;
	int precision = 6, pre = precision;

	cout << endl;
	cout << setw(2 * w1 + w2) << right << "REGRESSION MODEL REPORT" << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << "Dependent variable: " << d.a[yx[0] - 1] << endl;
	cout << "Method: Least Squares" << endl;
	cout << "Date & Time: " << ctime(&now);
	cout << "File Name: " << fileName << endl;
	cout << "Included observations: " << m << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << fixed << setprecision(pre);
	cout << setw(w1) << left << "Variable" << setw(w1) << right << "Coefficient" << setw(w2) << right << "t-Statistic" << setw(w3) << right << "p-value" << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << setw(w1) << left << "Constant" << setw(w1) << right << beta(1) << setw(w2) << right << tStat(1) << setw(w3) << right << tPval(1) << endl;
	for (int i = 1; i < p; i++)
	{
		cout << setw(w1) << left << d.a[yx[i] - 1] << setw(w1) << right << beta(i + 1) << setw(w2) << right << tStat(i + 1) << setw(w3) << right << tPval(i + 1) << endl;
	}
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << setw(w1) << left << "R_squared" << setw(w1) << right << R_squared << endl;
	cout << setw(w1) << left << "Adj. R_squared" << setw(w1) << right << adj_R_squared << endl;
	cout << setw(w1) << left << "F-Statistic" << setw(w1) << right << FStat << endl;
	cout << setw(w1) << left << "p-value(F-Stat)" << setw(w1) << right << FPval << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
}

double LR::yf(vector<double> vxf)
{
	Matrix xf(1, p, vxf);
	return (xf * beta)(1);
}

void LR::yflh(double &lo, double &hi, vector<double> vxf, double percent)
{
	Matrix xf(1, p, vxf);
	double yh = (xf * beta)(1);
	QR q(X.transpose() * X);
	double se = sigma * sqrt(1 + (xf * q.axb(xf.transpose()))(1));
	double half_alpha = (1 - percent / 100.0) / 2;
	lo = yh + tdistribution_quantile(half_alpha, m - p) * se;
	hi = yh - tdistribution_quantile(half_alpha, m - p) * se;
}

void LR::yblh(double &lo, double &hi, vector<double> vxf, double percent)
{
	Matrix xf(1, p, vxf);
	double yh = (xf * beta)(1);
	QR q(X.transpose() * X);
	double se = sigma * sqrt((xf * q.axb(xf.transpose()))(1));
	double half_alpha = (1 - percent / 100.0) / 2;
	lo = yh + tdistribution_quantile(half_alpha, m - p) * se;
	hi = yh - tdistribution_quantile(half_alpha, m - p) * se;
}

void LR::print(vector<double> vxf)
{

	time_t now = time(NULL);
	vector<double> percent{90, 95, 98, 99};
	vector<double> vyf, vyfl, vyfh, vybl, vybh;
	double yfl, yfh, ybl, ybh;
	for (int i = 0; i <= percent.size(); i++)
	{

		yflh(yfl, yfh, vxf, percent[i]);
		yblh(ybl, ybh, vxf, percent[i]);
		vyf.push_back(yf(vxf));
		vyfl.push_back(yfl);
		vyfh.push_back(yfh);
		vybl.push_back(ybl);
		vybh.push_back(ybh);
	}

	int columnWidth1 = 20, w1 = columnWidth1;
	int columnWidth2 = 20, w2 = columnWidth2;
	int columnWidth3 = 20, w3 = columnWidth3;
	int columnWidth4 = 20, w4 = columnWidth4;
	int width = w1 + w2 + w3 + w4;
	int precision = 6, pre = precision;

	cout << endl;
	cout << setw(w1 + w2 + w3) << right << "REGRESSION MODEL & FORECAST REPORT" << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << "Dependent variable: " << d.a[yx[0] - 1] << endl;
	cout << "Method: Least Squares" << endl;
	cout << "Date & Time: " << ctime(&now);
	cout << "File Name: " << fileName << endl;
	cout << "Included observations: " << m << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << fixed << setprecision(pre);
	cout << setw(w1) << left << "Variable" << setw(w2) << right << "Coefficient" << setw(w3) << "p-value" << setw(w4) << "Fcast Pt(Xf)" << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << setw(w1) << left << "Constant" << setw(w2) << right << beta(1) << setw(w3) << tPval(1) << setw(w4) << 1.0 << endl;
	for (int i = 1; i < p; i++)
	{
		cout << setw(w1) << left << d.a[yx[i] - 1] << setw(w2) << right << beta(i + 1) << setw(w3) << tPval(i + 1) << setw(w4) << vxf[i] << endl;
	}
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << setw(w1) << left << "Point(yf) Estimate" << setw(w2) << right << "% Confidence Lvl" << setw(w3) << " yf Range Est" << setw(w4) << "yf_bar Range Est" << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
	for (int i = 0; i < percent.size(); i++)
	{
		cout << setw(w1) << left << vyf[i] << setw(w2) << right << percent[i] << setw(w3) << vyfl[i] << setw(w4) << vybl[i] << endl;
		cout << setw(w1) << left << " " << setw(w2) << right << " " << setw(w3) << vyfh[i] << setw(w4) << vybh[i] << endl;
	}
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;

	cout << setw(w1) << left << "R_squared" << right << setw(w2) << R_squared << endl;
	cout << setw(w1) << left << "Adj. R_squared" << right << setw(w2) << adj_R_squared << endl;
	cout << setw(w1) << left << "F-Statistic" << right << setw(w2) << FStat << endl;
	cout << setw(w1) << left << "p-value(F-Stat)" << right << setw(w2) << FPval << endl;
	for (int i = 1; i <= width; i++)
	{
		cout << "=";
	}
	cout << endl;
}

#endif
