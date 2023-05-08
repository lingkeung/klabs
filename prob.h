#ifndef PROB_H
#define PROB_H

#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/QuantFuncMathMore.h"

using namespace std;

double normal_pdf(double x, double sigma=1, double mean=0);

double beta_pdf(double x, double a, double b);

double exponential_pdf(double x, double lamda, double x0=0);

double fdistribution_pdf(double x, double n, double m, double x0=0);

double tdistribution_pdf(double x, double r, double x0=0);

double chisquared_pdf(double x, double r, double x0=0);

double binomial_pdf(unsigned int k, double p, unsigned n);

double poisson_pdf(unsigned int k, double lamda);

double uniform_pdf(double x, double a, double b, double x0=0);

double normal_cdf(double x, double sigma=1, double mean=0);

double beta_cdf(double x, double a, double b);

double exponential_cdf(double x, double lamda, double x0=0);

double fdistribution_cdf(double x, double n, double m, double x0=0);

double tdistribution_cdf(double x, double r, double x0=0);

double chisquared_cdf(double x, double r, double x0=0);

double binomial_cdf(unsigned int k, double p, unsigned n);

double poisson_cdf(unsigned int k, double lamda);

double uniform_cdf(double x, double a, double b, double x0=0);

double normal_quantile(double z, double sigma=1);

double beta_quantile(double z, double a, double b);

double exponential_quantile(double z, double lamda);

double fdistribution_quantile(double z, double n, double m);

double tdistribution_quantile(double z, double r);

double chisquared_quantile(double z, double r);

double uniform_quantile(double z, double a, double b);

#endif