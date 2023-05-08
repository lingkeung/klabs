#include "rnt.h"
#include <chrono>
#include <random>
#include "matrix.h"
#include <ctime>

using namespace std;

void Stopwatch::start()
{
	start_time = std::chrono::system_clock::now();
}

void Stopwatch::stop()
{
	stop_time = std::chrono::system_clock::now();
}

double Stopwatch::etime()
{
	std::chrono::duration<double> diff = stop_time - start_time;
	return diff.count();
}

void delay(double ms)
{
	Stopwatch sw;
	sw.start();
	int i;
	for (i = 0; i < 1e100; i++)
	{
		sw.stop();
		double et = sw.etime();
		if (1000 * et >= ms)
		{
			break;
		}
	}
}

void Random::set_seed(int s)
{
	seed = s;
	generator.seed(seed);
}

void Random::set_seed()
{
	seed = time(NULL);
	generator.seed(seed);
}

double Random::real(char dist, double p1, double p2)
{
	double rv = 0;
	switch (dist)
	{
	case 'u':
	{
		uniform_real_distribution<double> u(p1, p2);
		rv = u(generator);
		break;
	}
	case 'n':
	{
		normal_distribution<double> n(p1, p2);
		rv = n(generator);
		break;
	}
	default:
		return rv;
	}
	return rv;
}

void Random::array(int size, char dist, double p1, double p2, double arr[])
{

	for (int i = 0; i < size; i++)
	{
		arr[i] = real(dist, p1, p2);
	}
}

Matrix Random::matrix(int m, int n, char dist, double p1, double p2)
{

	Matrix result(m, n);
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			result(i, j) = real(dist, p1, p2);
		}
	}
	return result;
}

int Random::integer(int a, int b)
{
	int rv = 0;
	uniform_int_distribution<int> i(a, b);
	rv = i(generator);
	return rv;
}

void Random::iarray(int size, int a, int b, int arr[])
{
	for (int i = 0; i < size; i++)
	{
		arr[i] = integer(a, b);
	}
}

int Random::poisson(double lamda)
{
	std::poisson_distribution<int> distribution(lamda);
	return distribution(generator);
}

void Random::parray(int size, double lamda, int arr[])
{
	for (int i = 0; i < size; i++)
	{
		arr[i] = poisson(lamda);
	}
}

bool Random::bernoulli(double prob)
{
	std::bernoulli_distribution distribution(prob);
	return distribution(generator);
}

void Random::bearray(int size, double prob, bool arr[])
{
	for (int i = 0; i < size; i++)
	{
		arr[i] = bernoulli(prob);
	}
}

int Random::binomial(int n, double p)
{
	std::binomial_distribution<int> distribution(n, p);
	return distribution(generator);
}

void Random::biarray(int size, int n, double p, int arr[])
{
	for (int i = 0; i < size; i++)
	{
		arr[i] = binomial(n, p);
	}
}

Matrix Random::pmatrix(int m, int n, double lamda)
{
	Matrix out(m, n);
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			out(i, j) = poisson(lamda);
		}
	}
	return out;
}

Matrix Random::bematrix(int m, int n, double prob)
{
	Matrix out(m, n);
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			out(i, j) = bernoulli(prob);
		}
	}
	return out;
}

Matrix Random::bimatrix(int s, int t, int n, double p)
{
	Matrix out(s, t);
	for (int i = 1; i <= s; i++)
	{
		for (int j = 1; j <= t; j++)
		{
			out(i, j) = binomial(n, p);
		}
	}
	return out;
}