#include "prob.h"

using namespace std;

double normal_pdf(double x, double sigma, double mean)
{
	return ROOT::Math::normal_pdf(x, sigma, mean);
}

double beta_pdf(double x, double a, double b)
{
	return ROOT::Math::beta_pdf(x, a, b);
}

double exponential_pdf(double x, double lamda, double x0)
{
	return ROOT::Math::exponential_pdf(x, lamda, x0);
}

double fdistribution_pdf(double x, double n, double m, double x0)
{
	return ROOT::Math::fdistribution_pdf(x, n, m, x0);
}

double tdistribution_pdf(double x, double r, double x0)
{
	return ROOT::Math::tdistribution_pdf(x, r, x0);
}

double chisquared_pdf(double x, double r, double x0)
{
	return ROOT::Math::chisquared_pdf(x, r, x0);
}

double binomial_pdf(unsigned int k, double p, unsigned n)
{
	return ROOT::Math::binomial_pdf(k, p, n);
}

double poisson_pdf(unsigned int k, double lamda)
{
	return ROOT::Math::poisson_pdf(k, lamda);
}

double uniform_pdf(double x, double a, double b, double x0)
{
	return ROOT::Math::uniform_pdf(x, a, b);
}

double normal_cdf(double x, double sigma, double mean)
{
	return ROOT::Math::normal_cdf(x, sigma, mean);
}

double beta_cdf(double x, double a, double b)
{
	return ROOT::Math::beta_cdf(x, a, b);
}

double exponential_cdf(double x, double lamda, double x0)
{
	return ROOT::Math::exponential_cdf(x, lamda, x0);
}

double fdistribution_cdf(double x, double n, double m, double x0)
{
	return ROOT::Math::fdistribution_cdf(x, n, m, x0);
}

double tdistribution_cdf(double x, double r, double x0)
{
	return ROOT::Math::tdistribution_cdf(x, r, x0);
}

double chisquared_cdf(double x, double r, double x0)
{
	return ROOT::Math::chisquared_cdf(x, r, x0);
}

double binomial_cdf(unsigned int k, double p, unsigned n)
{
	return ROOT::Math::binomial_cdf(k, p, n);
}

double poisson_cdf(unsigned int k, double lamda)
{
	return ROOT::Math::poisson_cdf(k, lamda);
}

double uniform_cdf(double x, double a, double b, double x0)
{
	return ROOT::Math::uniform_cdf(x, a, b);
}

double normal_quantile(double z, double sigma)
{
	return ROOT::Math::gaussian_quantile(z, sigma);
}

double beta_quantile(double z, double a, double b)
{
	return ROOT::Math::beta_quantile(z, a, b);
}

double exponential_quantile(double z, double lamda)
{
	return ROOT::Math::exponential_quantile(z, lamda);
}

double fdistribution_quantile(double z, double n, double m)
{
	return ROOT::Math::fdistribution_quantile(z, n, m);
}

double tdistribution_quantile(double z, double r)
{
	return ROOT::Math::tdistribution_quantile(z, r);
}

double chisquared_quantile(double z, double r)
{
	return ROOT::MathMore::chisquared_quantile(z, r);
}

double uniform_quantile(double z, double a, double b)
{
	return ROOT::Math::uniform_quantile(z, a, b);
}
