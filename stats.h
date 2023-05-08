#ifndef STATS_H
#define STATS_H

#include "matrix.h"

using namespace std;


int count(Matrix data);

double maximum(Matrix data);

double minimum(Matrix data);

double mean(Matrix data);

double variance(Matrix data);

double stdev(Matrix data);

double range(Matrix data);

double median(Matrix data);

double percentile(double p, Matrix data);

double skew(Matrix data);

double skew2(Matrix data);

double kurt(Matrix data);

double kurt2(Matrix data);

void statsPrint(Matrix data);

double cov(Matrix x, Matrix y); 

double corr(Matrix x, Matrix y);

Matrix stdz(Matrix A);

Matrix cov(Matrix A); 

Matrix corr(Matrix A);


#endif