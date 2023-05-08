#ifndef RNT_H
#define RNT_H
#include <chrono>
#include <random>
#include "matrix.h"

using namespace std;

class Stopwatch
{
private:
	chrono::system_clock::time_point start_time;
	chrono::system_clock::time_point stop_time;

public:
	void start();
	void stop();
	double etime();
	string gtu() { return time_unit; }

private:
	string time_unit = " sec ";
};

class Random
{
private:
	default_random_engine generator;
	int seed;

public:
	int get_seed() { return seed; }
	void set_seed(int s);
	void set_seed();

	double real(char dist, double p1, double p2);
	Matrix matrix(int m, int n, char dist, double p1, double p2);
	void array(int size, char dist, double p1, double p2, double arr[]);

	int integer(int a, int b);
	void iarray(int size, int a, int b, int arr[]);

	int poisson(double lamda);
	void parray(int size, double lamda, int arr[]);

	bool bernoulli(double prob);
	void bearray(int size, double prob, bool arr[]);

	int binomial(int n, double p);
	void biarray(int size, int n, double p, int arr[]);

	Matrix pmatrix(int m, int n, double lamda);
	Matrix bematrix(int m, int n, double prob);
	Matrix bimatrix(int s, int t, int n, double p);
};

void delay(double ms);

#endif