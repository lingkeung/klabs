#ifndef EQN_H
#define EQN_H
#include <vector>

using namespace std;

double newton(double f(double), double df(double), double x0, double tol = 1e-6);

vector<vector<double>> euler(double t0, double x0, double tend, int n, double xprime(double, double));

vector<vector<double>> euler(double t0, double x10, double x20, double tend, int n, vector<double (*)(double, double, double)> vf);

vector<vector<double>> rk4(double t0, double x0, double tend, int n, double xprime(double, double));

vector<vector<double>> rk4(double t0, double x10, double x20, double tend, int n, vector<double (*)(double, double, double)> vf);

vector<vector<double>> rk4(double t0, double x10, double x20, double x30, double tend, int n, vector<double (*)(double, double, double, double)> vf);

vector<vector<double>> rk4(double t0, double x10, double x20, double x30, double x40, double tend, int n, vector<double (*)(double, double, double, double, double)> vf);

vector<vector<double>> rk4(double t0, vector<double> X0, double tend, int n, vector<double (*)(double, double, double, double, double)> vf);

vector<vector<double>> rk4(double t0, vector<double> X0, double tend, int n, vector<double (*)(double, double, double, double, double, double, double, double, double)> vf);

vector<vector<double>> rk4(double t0, vector<double> X0, double tend, int n, vector<double (*)(double, double, double, double, double, double, double, double, double, double, double, double, double)> vf);

#endif