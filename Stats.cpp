#include "stats.h"
#include <vector>
#include "matrix.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include "prob.h"

using namespace std;

Matrix stdz(Matrix A)
{

	int m = A.getM(), n = A.getN();
	vector<double> means;
	vector<double> stds;
	for (int j = 1; j <= n; j++)
	{
		double sum = 0;
		for (int i = 1; i <= m; i++)
		{
			sum += A(i, j);
		}
		means.push_back(sum / m);
	}
	for (int j = 1; j <= n; j++)
	{
		double sum = 0;
		for (int i = 1; i <= m; i++)
		{
			sum += pow((A(i, j) - means[j - 1]), 2);
		}
		stds.push_back(sqrt(sum / (m - 1)));
	}

	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			A(i, j) = (A(i, j) - means[j - 1]) / stds[j - 1];
		}
	}
	return A;
}

Matrix cov(Matrix A)
{

	int m = A.getM(), n = A.getN();
	vector<double> means;
	vector<double> stds;
	for (int j = 1; j <= n; j++)
	{
		double sum = 0;
		for (int i = 1; i <= m; i++)
		{
			sum += A(i, j);
		}
		means.push_back(sum / m);
	}
	for (int j = 1; j <= n; j++)
	{
		double sum = 0;
		for (int i = 1; i <= m; i++)
		{
			sum += pow((A(i, j) - means[j - 1]), 2);
		}
		stds.push_back(sqrt(sum / (m - 1)));
	}

	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			A(i, j) = A(i, j) - means[j - 1];
		}
	}
	return (1.0 / (m - 1)) * (A.transpose() * A);
}

Matrix corr(Matrix A)
{

	int m = A.getM(), n = A.getN();
	vector<double> means;
	vector<double> stds;
	for (int j = 1; j <= n; j++)
	{
		double sum = 0;
		for (int i = 1; i <= m; i++)
		{
			sum += A(i, j);
		}
		means.push_back(sum / m);
	}
	for (int j = 1; j <= n; j++)
	{
		double sum = 0;
		for (int i = 1; i <= m; i++)
		{
			sum += pow((A(i, j) - means[j - 1]), 2);
		}
		stds.push_back(sqrt(sum / (m - 1)));
	}
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= m; i++)
		{
			A(i, j) = (A(i, j) - means[j - 1]) / stds[j - 1];
		}
	}
	Matrix Corr = (1.0 / (m - 1)) * (A.transpose() * A);

	double T, rho;
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j < i; j++)
		{
			rho = Corr(i, j);
			T = rho / sqrt((1 - rho * rho) / (m - 2));
			Corr(j, i) = 2 * (1 - tdistribution_cdf(T, m - 2));
		}
	}
	return Corr;
}

int count(Matrix data)
{
	return data.getM() * data.getN();
}

double maximum(Matrix data)
{
	vector<double> vdbl = data.getVdbl();
	return *(std::max_element(vdbl.begin(), vdbl.end()));
}

double minimum(Matrix data)
{
	vector<double> vdbl = data.getVdbl();
	return *(std::min_element(vdbl.begin(), vdbl.end()));
}

double mean(Matrix data)
{
	vector<double> vdbl = data.getVdbl();
	double sum = std::accumulate(vdbl.begin(), vdbl.end(), 0.0);
	return sum / (data.getM() * data.getN());
}

double variance(Matrix data)
{
	int m = data.getM(), n = data.getN();
	vector<double> v = data.getVdbl();
	double average = mean(data);
	for_each(v.begin(), v.end(), [a = average](double &x)
			 { x = (x) * (x); });
	double sum = std::accumulate(v.begin(), v.end(), 0.0) - m * n * average * average;
	return sum / (m * n - 1);
}

double stdev(Matrix data)
{
	return sqrt(variance(data));
}

double range(Matrix data)
{
	return maximum(data) - minimum(data);
}

double median(Matrix data)
{
	int n = count(data);
	vector<double> v = data.getVdbl();
	std::sort(v.begin(), v.end());
	return (n % 2 == 0) ? (v[n / 2 - 1] + v[n / 2]) / 2 : v[(n + 1) / 2 - 1];
}

double percentile(double p, Matrix data)
{
	int n = count(data);
	vector<double> v = data.getVdbl();
	std::sort(v.begin(), v.end());
	double P = (n - 1) * p / 100;
	int F = floor(P), C = ceil(P);
	if (abs(P - F) < 1e-6)
		return v[F];
	else if (abs(C - P) < 1e-6)
		return v[C];
	else
	{
		double f = P - F;
		return v[F] + f * (v[F + 1] - v[F]);
	}
}

double skew(Matrix data)
{
	int n = count(data);
	double avg = mean(data);
	vector<double> v3 = data.getVdbl();
	for_each(v3.begin(), v3.end(), [=](double &x)
			 { x = (x * x * x); });
	vector<double> v2 = data.getVdbl();
	for_each(v2.begin(), v2.end(), [=](double &x)
			 { x = x * x; });
	vector<double> v1 = data.getVdbl();
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum = sum + v3[i] - 3 * avg * v2[i] + 2 * v1[i] * pow(avg, 2);
	}
	return sum / pow(stdev(data), 3) * n / (n - 1) / (n - 2);
}

double skew2(Matrix data)
{
	int n = count(data);
	double avg = mean(data);
	vector<double> v = data.getVdbl();

	for_each(v.begin(), v.end(), [=](double &x)
			 { x = pow((x - avg), 3); });
	return std::accumulate(v.begin(), v.end(), 0.0) / pow(stdev(data), 3) * n / (n - 1) / (n - 2);
}

double kurt(Matrix data)
{
	int n = count(data);
	double avg = mean(data);
	vector<double> v4 = data.getVdbl();
	for_each(v4.begin(), v4.end(), [=](double &x)
			 { x = (x * x * x * x); });
	vector<double> v3 = data.getVdbl();
	for_each(v3.begin(), v3.end(), [=](double &x)
			 { x = (x * x * x); });
	vector<double> v2 = data.getVdbl();
	for_each(v2.begin(), v2.end(), [=](double &x)
			 { x = x * x; });
	vector<double> v1 = data.getVdbl();
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum = sum + v4[i] - 4 * avg * v3[i] + 6 * avg * avg * v2[i] - 3 * avg * avg * avg * v1[i];
	}
	double s = stdev(data);
	return sum / pow(s, 4) * n * (n + 1) / (n - 1) / (n - 2) / (n - 3) - 3 * (n - 1) * (n - 1) / (n - 2) / (n - 3);
}

double kurt2(Matrix data)
{
	int n = count(data);
	double avg = mean(data);
	vector<double> v = data.getVdbl();

	for_each(v.begin(), v.end(), [=](double &x)
			 { x = pow((x - avg), 4); });
	return std::accumulate(v.begin(), v.end(), 0.0) / pow(stdev(data), 4) * n * (n + 1) / (n - 1) / (n - 2) / (n - 3) - 3 * (n - 1) * (n - 1) / (n - 2) / (n - 3);
}

void statsPrint(Matrix data)
{

	cout << fixed << setprecision(2) << endl;
	cout << setw(20) << "Count" << setw(20) << count(data) << endl;
	cout << setw(20) << "Min" << setw(20) << minimum(data) << endl;
	cout << setw(20) << "Q1(25%)" << setw(20) << percentile(25, data) << endl;
	cout << setw(20) << "Mean" << setw(20) << mean(data) << endl;
	cout << setw(20) << "Median" << setw(20) << percentile(50, data) << endl;
	cout << setw(20) << "Q3(75%)" << setw(20) << percentile(75, data) << endl;
	cout << setw(20) << "Max" << setw(20) << maximum(data) << endl;
	cout << setw(20) << "Range" << setw(20) << range(data) << endl;
	cout << setw(20) << "Var" << setw(20) << variance(data) << endl;
	cout << setw(20) << "Std" << setw(20) << stdev(data) << endl;
	cout << setw(20) << "Skew" << setw(20) << skew(data) << endl;
	cout << setw(20) << "Kurt" << setw(20) << kurt(data) << endl;
}

double cov(Matrix x, Matrix y)
{

	double result = 0;
	if (count(x) != count(y))
		return result;
	else
	{
		double xbar = mean(x), ybar = mean(y);
		int size = count(x);
		for (int i = 1; i <= size; i++)
		{
			result = result + (x(i) - xbar) * (y(i) - ybar);
		}
		return result / (size - 1);
	}
}

double corr(Matrix x, Matrix y)
{
	return cov(x, y) / stdev(x) / stdev(y);
}