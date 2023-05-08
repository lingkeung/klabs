#ifndef PCA_H
#define PCA_H

#include <vector>
#include "stats.h"
#include <iostream>
#include "data.h"
#include "svd.h"
#include "svdc.h"
#include <cmath>
#include <vector>
#include <string>

using namespace std;

class PCA
{
public:
	string fileName;
	int m, n;
	vector<string> attribute;
	Matrix D;
	Matrix X;
	Matrix sigmaX;
	Matrix totalLoadingFactor;
	Matrix xLoadingFactor;
	Matrix component;
	Matrix componentVariance;
	Matrix varianceExplained;
	Matrix varianceCumulative;

public:
	PCA(string fN);
	void pCompRank();
	void pCompXcorr(int kComponents);
};

PCA::PCA(string fN)
{

	fileName = fN;
	Data d(fileName);
	m = d.m;
	n = d.n;
	attribute = d.a;
	D = d.D;
	X = stdz(D);
	sigmaX = (1.0 / (m - 1)) * (X.transpose() * X); // covariance matrix

	SVD s(sigmaX);
	totalLoadingFactor = s.getU(); // get eigenvectors matrix

	component = X * totalLoadingFactor; // apply loading

	componentVariance = Matrix(1, n);
	for (int i = 1; i <= n; i++)
	{
		componentVariance(i) = cov(component)(i, i); // (same as eignenvalues)
	}

	xLoadingFactor = Matrix(n, n);
	for (int i = 1; i <= n; i++)
	{
		xLoadingFactor.setblk(1, i, sqrt(componentVariance(i)) * totalLoadingFactor(1, n, i, i));
	}
	xLoadingFactor = xLoadingFactor.transpose();

	double sum = 0;
	for (int i = 1; i <= n; i++)
	{
		sum += componentVariance(i);
	}
	varianceExplained = Matrix(1, n);
	varianceExplained = (1 / sum) * componentVariance;
	varianceExplained = 100.0 * varianceExplained;

	varianceCumulative = Matrix(1, n);
	varianceCumulative(1) = varianceExplained(1);
	for (int i = 2; i <= n; i++)
	{
		varianceCumulative(i) = varianceCumulative(i - 1) + varianceExplained(i);
	}
}

void PCA::pCompRank()
{

	int w1 = 9, w2 = 3, w3 = 15, w4 = 20, w5 = 20;
	double pre = 6;
	int w = w1 + w2 + w3 + w4 + w5;
	cout << fixed << setprecision(pre) << endl;
	cout << endl;
	cout << setw(w1 + w2 + w3 + w4) << right << "Principal Components Variance" << endl;
	;
	for (int i = 1; i <= w; i++)
	{
		cout << "=";
	}
	cout << endl;
	cout << setw(w1) << left << "Component" << setw(w2) << right << "  " << setw(w3) << "Variance" << setw(w4) << "% Explained" << setw(w5) << "% Cumulative" << endl;
	for (int i = 1; i <= w; i++)
	{
		cout << "=";
	}
	cout << endl;
	for (int i = 1; i <= n; i++)
	{
		cout << setw(w1) << left << "Component" << setw(w2) << right << i << setw(w3) << componentVariance(i) << setw(w4) << varianceExplained(i)
			 << setw(w5) << varianceCumulative(i) << endl;
	}
	for (int i = 1; i <= w; i++)
	{
		cout << "=";
	}
	cout << endl;
}

// method to print xLoadingFactor
void PCA::pCompXcorr(int kComponents)
{

	int k = kComponents;

	// append cumulative contributions
	Matrix output(k + 1, n);
	output.setblk(1, 1, xLoadingFactor(1, k, 1, n));
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= k; i++)
		{
			output(k + 1, j) = output(k + 1, j) + output(i, j) * output(i, j);
		}
	}

	// print appended matrix in a nice table
	cout << fixed << setprecision(6) << endl;

	int columnWidth0 = 15, w0 = columnWidth0;
	int columnWidth = 15, w = columnWidth;
	int intColWidth = 2, wi = intColWidth;
	int tblWidth = w0 + wi + n * w; // width of table with 2+n columns
	const char(*title) = "Correlation(Component, X) Contribution";

	cout << setw(w0 + wi + n * w / 2) << right << title << endl;

	for (int i = 1; i <= tblWidth; i++)
	{
		cout << "=";
	}
	cout << endl;

	cout << setw(w0) << left << " " << setw(wi) << " ";
	for (int i = 0; i < n; i++)
	{
		cout << setw(w) << right << attribute[i];
	}
	cout << endl;

	for (int i = 1; i <= tblWidth; i++)
	{
		cout << "=";
	}
	cout << endl;

	for (int i = 1; i <= k; i++)
	{
		cout << setw(w0) << left << "Component " << setw(wi) << i;
		for (int j = 1; j <= n; j++)
		{
			cout << setw(w) << right << output(i, j);
		}
		cout << endl;
	}

	for (int i = 1; i <= tblWidth; i++)
	{
		cout << "=";
	}
	cout << endl;

	cout << setw(w0) << left << "% Explained" << setw(wi) << " ";
	for (int i = 1; i <= n; i++)
	{
		cout << setw(w) << right << 100 * output(k + 1, i);
	}
	cout << endl;

	for (int i = 1; i <= tblWidth; i++)
	{
		cout << "=";
	}
	cout << endl;
}

#endif