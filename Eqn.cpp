#include "eqn.h"
#include <vector>
#include <cmath>

using namespace std;

double newton(double f(double), double df(double), double x0, double tol)
{
    double x = x0;
    double x1 = 0;

    for (int i = 0; i < 100; i++)
    {
        x1 = x - f(x) / df(x);
        if (abs(x1 - x) < tol)
        {
            return x1;
        }
        x = x1;
    }
    return x1;
}

vector<vector<double>> euler(double t0, double x0, double tend, int n, double xprime(double, double))
{

    vector<double> t(n, 0);
    vector<double> x(n, 0);
    t[0] = t0;
    x[0] = x0;
    double dt = (tend - t0) / (n - 1);
    for (int i = 1; i < n; i++)
    {
        t[i] = t[i - 1] + dt;
        x[i] = x[i - 1] + xprime(t[i - 1], x[i - 1]) * dt;
    }
    vector<vector<double>> result = {t, x};
    return result;
}

vector<vector<double>> euler(double t0, double x10, double x20, double tend, int n, vector<double (*)(double, double, double)> vf)
{

    vector<double> t(n, 0);
    vector<double> x1(n, 0);
    vector<double> x2(n, 0);
    t[0] = t0;
    x1[0] = x10;
    x2[0] = x20;
    double dt = (tend - t0) / (n - 1);
    for (int i = 1; i < n; i++)
    {
        t[i] = t[i - 1] + dt;
        x1[i] = x1[i - 1] + vf[0](t[i - 1], x1[i - 1], x2[i - 1]) * dt;
        x2[i] = x2[i - 1] + vf[1](t[i - 1], x1[i - 1], x2[i - 1]) * dt;
    }
    vector<vector<double>> result = {t, x1, x2};
    return result;
}

vector<vector<double>> rk4(double t0, double x0, double tend, int n, double xprime(double, double))
{

    vector<double> t(n, 0);
    vector<double> x(n, 0);
    t[0] = t0;
    x[0] = x0;
    double dt = (tend - t0) / (n - 1);
    for (int i = 1; i < n; i++)
    {
        t[i] = t[i - 1] + dt;
        double K1 = xprime(t[i - 1], x[i - 1]);
        double K2 = xprime((t[i - 1] + dt / 2), (x[i - 1] + dt * K1 / 2));
        double K3 = xprime((t[i - 1] + dt / 2), (x[i - 1] + dt * K2 / 2));
        double K4 = xprime((t[i - 1] + dt), (x[i - 1] + dt * K3));
        x[i] = x[i - 1] + dt / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
    }
    vector<vector<double>> result = {t, x};
    return result;
}

vector<vector<double>> rk4(double t0, double x10, double x20, double tend, int n, vector<double (*)(double, double, double)> vf)
{

    vector<double> t(n, 0);
    vector<double> x1(n, 0);
    vector<double> x2(n, 0);
    t[0] = t0;
    x1[0] = x10;
    x2[0] = x20;
    double dt = (tend - t0) / (n - 1);
    vector<double> K1(2, 0);
    vector<double> K2(2, 0);
    vector<double> K3(2, 0);
    vector<double> K4(2, 0);
    for (int i = 1; i < n; i++)
    {
        t[i] = t[i - 1] + dt;
        for (int j = 0; j < 2; j++)
        {
            double k1 = dt * vf[j](t[i - 1], x1[i - 1], x2[i - 1]);
            K1[j] = k1;
            double k2 = dt * vf[j]((t[i - 1] + dt / 2), (x1[i - 1] + K1[0] / 2), (x2[i - 1] + K1[1] / 2));
            K2[j] = k2;
            double k3 = dt * vf[j]((t[i - 1] + dt / 2), (x1[i - 1] + K2[0] / 2), (x2[i - 1] + K2[1] / 2));
            K3[j] = k3;
            double k4 = dt * vf[j]((t[i - 1] + dt), (x1[i - 1] + K3[0]), (x2[i - 1] + K3[1]));
            K4[j] = k4;
        }
        x1[i] = x1[i - 1] + (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]) / 6;
        x2[i] = x2[i - 1] + (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]) / 6;
    }
    vector<vector<double>> result = {t, x1, x2};
    return result;
}

vector<vector<double>> rk4(double t0, double x10, double x20, double x30, double tend, int n,
                           vector<double (*)(double, double, double, double)> vf)
{

    vector<double> t(n, 0);
    vector<double> x1(n, 0);
    vector<double> x2(n, 0);
    vector<double> x3(n, 0);
    t[0] = t0;
    x1[0] = x10;
    x2[0] = x20;
    x3[0] = x30;
    double dt = (tend - t0) / (n - 1);
    vector<double> K1(3, 0);
    vector<double> K2(3, 0);
    vector<double> K3(3, 0);
    vector<double> K4(3, 0);
    for (int i = 1; i < n; i++)
    {
        t[i] = t[i - 1] + dt;
        for (int j = 0; j < 3; j++)
        {
            double k1 = vf[j](t[i - 1], x1[i - 1], x2[i - 1], x3[i - 1]);
            K1[j] = k1;
            double k2 = vf[j]((t[i - 1] + dt / 2), (x1[i - 1] + dt * K1[0] / 2), (x2[i - 1] + dt * K1[1] / 2),
                              (x3[i - 1] + dt * K1[2] / 2));
            K2[j] = k2;
            double k3 = vf[j]((t[i - 1] + dt / 2), (x1[i - 1] + dt * K2[0] / 2), (x2[i - 1] + dt * K2[1] / 2),
                              (x3[i - 1] + dt * K2[2] / 2));
            K3[j] = k3;
            double k4 = vf[j]((t[i - 1] + dt), (x1[i - 1] + dt * K3[0]), (x2[i - 1] + dt * K3[1]),
                              (x3[i - 1] + dt * K3[2]));
            K4[j] = k4;
        }
        x1[i] = x1[i - 1] + dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
        x2[i] = x2[i - 1] + dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
        x3[i] = x3[i - 1] + dt / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
    }
    vector<vector<double>> result = {t, x1, x2, x3};
    return result;
}

vector<vector<double>> rk4(double t0, double x10, double x20, double x30, double x40, double tend, int n, vector<double (*)(double, double, double, double, double)> vf)
{

    vector<double> t(n, 0);
    vector<double> x1(n, 0);
    vector<double> x2(n, 0);
    vector<double> x3(n, 0);
    vector<double> x4(n, 0);
    t[0] = t0;
    x1[0] = x10;
    x2[0] = x20;
    x3[0] = x30;
    x4[0] = x40;
    double dt = (tend - t0) / (n - 1);
    vector<double> K1(4, 0);
    vector<double> K2(4, 0);
    vector<double> K3(4, 0);
    vector<double> K4(4, 0);
    for (int i = 1; i < n; i++)
    {
        t[i] = t[i - 1] + dt;
        for (int j = 0; j < 4; j++)
        {
            double k1 = vf[j](t[i - 1], x1[i - 1], x2[i - 1], x3[i - 1], x4[i - 1]);
            K1[j] = k1;
            double k2 = vf[j]((t[i - 1] + dt / 2), (x1[i - 1] + dt * K1[0] / 2), (x2[i - 1] + dt * K1[1] / 2),
                              (x3[i - 1] + dt * K1[2] / 2), (x4[i - 1] + dt * K1[3] / 2));
            K2[j] = k2;
            double k3 = vf[j]((t[i - 1] + dt / 2), (x1[i - 1] + dt * K2[0] / 2), (x2[i - 1] + dt * K2[1] / 2),
                              (x3[i - 1] + dt * K2[2] / 2), (x4[i - 1] + dt * K2[3] / 2));
            K3[j] = k3;
            double k4 = vf[j]((t[i - 1] + dt), (x1[i - 1] + dt * K3[0]), (x2[i - 1] + dt * K3[1]),
                              (x3[i - 1] + dt * K3[2]), (x4[i - 1] + dt * K3[3]));
            K4[j] = k4;
        }
        x1[i] = x1[i - 1] + dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
        x2[i] = x2[i - 1] + dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
        x3[i] = x3[i - 1] + dt / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
        x4[i] = x4[i - 1] + dt / 6 * (K1[3] + 2 * K2[3] + 2 * K3[3] + K4[3]);
    }
    vector<vector<double>> result = {t, x1, x2, x3, x4};
    return result;
}

vector<vector<double>> rk4(double t0, vector<double> X0, double tend, int n,
                           vector<double (*)(double, double, double, double, double)> vf)
{

    vector<vector<double>> x(5, vector<double>(n, 0));
    x[0][0] = t0;
    for (int i = 1; i <= 4; i++)
    {
        x[i][0] = X0[i - 1];
    }
    double dt = (tend - t0) / (n - 1);
    vector<double> K1(4, 0);
    vector<double> K2(4, 0);
    vector<double> K3(4, 0);
    vector<double> K4(4, 0);

    for (int i = 1; i < n; i++)
    {
        x[0][i] = x[0][i - 1] + dt;
        for (int j = 0; j < 4; j++)
        {
            double k1 = vf[j](x[0][i - 1], x[1][i - 1], x[2][i - 1], x[3][i - 1], x[4][i - 1]);
            K1[j] = k1;
            double k2 = vf[j]((x[0][i - 1] + dt / 2), (x[1][i - 1] + dt * K1[0] / 2), (x[2][i - 1] + dt * K1[1] / 2),
                              (x[3][i - 1] + dt * K1[2] / 2), (x[4][i - 1] + dt * K1[3] / 2));
            K2[j] = k2;
            double k3 = vf[j]((x[0][i - 1] + dt / 2), (x[1][i - 1] + dt * K2[0] / 2), (x[2][i - 1] + dt * K2[1] / 2),
                              (x[3][i - 1] + dt * K2[2] / 2), (x[4][i - 1] + dt * K2[3] / 2));
            K3[j] = k3;
            double k4 = vf[j]((x[0][i - 1] + dt), (x[1][i - 1] + dt * K3[0]), (x[2][i - 1] + dt * K3[1]),
                              (x[3][i - 1] + dt * K3[2]), (x[4][i - 1] + dt * K3[3]));
            K4[j] = k4;
        }
        x[1][i] = x[1][i - 1] + dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
        x[2][i] = x[2][i - 1] + dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
        x[3][i] = x[3][i - 1] + dt / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
        x[4][i] = x[4][i - 1] + dt / 6 * (K1[3] + 2 * K2[3] + 2 * K3[3] + K4[3]);
    }
    return x;
}

vector<vector<double>> rk4(double t0, vector<double> X0, double tend, int n,
                           vector<double (*)(double, double, double, double, double, double, double, double, double)> vf)
{

    vector<vector<double>> x(9, vector<double>(n, 0));
    x[0][0] = t0;
    for (int i = 1; i <= 8; i++)
    {
        x[i][0] = X0[i - 1];
    }
    double dt = (tend - t0) / (n - 1);
    vector<double> K1(8, 0);
    vector<double> K2(8, 0);
    vector<double> K3(8, 0);
    vector<double> K4(8, 0);

    for (int i = 1; i < n; i++)
    {
        x[0][i] = x[0][i - 1] + dt;
        for (int j = 0; j < 8; j++)
        {
            double k1 = vf[j](x[0][i - 1], x[1][i - 1], x[2][i - 1], x[3][i - 1], x[4][i - 1],
                              x[5][i - 1], x[6][i - 1], x[7][i - 1], x[8][i - 1]);
            K1[j] = k1;
            double k2 = vf[j]((x[0][i - 1] + dt / 2), (x[1][i - 1] + dt * K1[0] / 2), (x[2][i - 1] + dt * K1[1] / 2),
                              (x[3][i - 1] + dt * K1[2] / 2), (x[4][i - 1] + dt * K1[3] / 2),
                              (x[5][i - 1] + dt * K1[4] / 2), (x[6][i - 1] + dt * K1[5] / 2),
                              (x[7][i - 1] + dt * K1[6] / 2), (x[8][i - 1] + dt * K1[7] / 2));
            K2[j] = k2;
            double k3 = vf[j]((x[0][i - 1] + dt / 2), (x[1][i - 1] + dt * K2[0] / 2), (x[2][i - 1] + dt * K2[1] / 2),
                              (x[3][i - 1] + dt * K2[2] / 2), (x[4][i - 1] + dt * K2[3] / 2),
                              (x[5][i - 1] + dt * K2[4] / 2), (x[6][i - 1] + dt * K2[5] / 2),
                              (x[7][i - 1] + dt * K2[6] / 2), (x[8][i - 1] + dt * K2[7] / 2));
            K3[j] = k3;
            double k4 = vf[j]((x[0][i - 1] + dt), (x[1][i - 1] + dt * K3[0]), (x[2][i - 1] + dt * K3[1]),
                              (x[3][i - 1] + dt * K3[2]), (x[4][i - 1] + dt * K3[3]),
                              (x[5][i - 1] + dt * K3[4]), (x[6][i - 1] + dt * K3[5]),
                              (x[7][i - 1] + dt * K3[6]), (x[8][i - 1] + dt * K3[7]));
            K4[j] = k4;
        }
        x[1][i] = x[1][i - 1] + dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
        x[2][i] = x[2][i - 1] + dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
        x[3][i] = x[3][i - 1] + dt / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
        x[4][i] = x[4][i - 1] + dt / 6 * (K1[3] + 2 * K2[3] + 2 * K3[3] + K4[3]);
        x[5][i] = x[5][i - 1] + dt / 6 * (K1[4] + 2 * K2[4] + 2 * K3[4] + K4[4]);
        x[6][i] = x[6][i - 1] + dt / 6 * (K1[5] + 2 * K2[5] + 2 * K3[5] + K4[5]);
        x[7][i] = x[7][i - 1] + dt / 6 * (K1[6] + 2 * K2[6] + 2 * K3[6] + K4[5]);
        x[8][i] = x[8][i - 1] + dt / 6 * (K1[7] + 2 * K2[7] + 2 * K3[7] + K4[7]);
    }
    return x;
}

vector<vector<double>> rk4(double t0, vector<double> X0, double tend, int n,
                           vector<double (*)(double, double, double, double, double, double, double, double, double,
                                             double, double, double, double)>
                               vf)
{

    vector<vector<double>> x(13, vector<double>(n, 0));
    x[0][0] = t0;
    for (int i = 1; i <= 12; i++)
    {
        x[i][0] = X0[i - 1];
    }
    double dt = (tend - t0) / (n - 1);
    vector<double> K1(12, 0);
    vector<double> K2(12, 0);
    vector<double> K3(12, 0);
    vector<double> K4(12, 0);
    for (int i = 1; i < n; i++)
    {
        x[0][i] = x[0][i - 1] + dt;
        for (int j = 0; j < 8; j++)
        {
            double k1 = vf[j](x[0][i - 1], x[1][i - 1], x[2][i - 1], x[3][i - 1], x[4][i - 1],
                              x[5][i - 1], x[6][i - 1], x[7][i - 1], x[8][i - 1],
                              x[9][i - 1], x[10][i - 1], x[11][i - 1], x[12][i - 1]);
            K1[j] = k1;
            double k2 = vf[j]((x[0][i - 1] + dt / 2), (x[1][i - 1] + dt * K1[0] / 2), (x[2][i - 1] + dt * K1[1] / 2),
                              (x[3][i - 1] + dt * K1[2] / 2), (x[4][i - 1] + dt * K1[3] / 2),
                              (x[5][i - 1] + dt * K1[4] / 2), (x[6][i - 1] + dt * K1[5] / 2),
                              (x[7][i - 1] + dt * K1[6] / 2), (x[8][i - 1] + dt * K1[7] / 2),
                              (x[9][i - 1] + dt * K1[8] / 2), (x[10][i - 1] + dt * K1[9] / 2),
                              (x[11][i - 1] + dt * K1[10] / 2), (x[12][i - 1] + dt * K1[11] / 2));
            K2[j] = k2;
            double k3 = vf[j]((x[0][i - 1] + dt / 2), (x[1][i - 1] + dt * K2[0] / 2), (x[2][i - 1] + dt * K2[1] / 2),
                              (x[3][i - 1] + dt * K2[2] / 2), (x[4][i - 1] + dt * K2[3] / 2),
                              (x[5][i - 1] + dt * K2[4] / 2), (x[6][i - 1] + dt * K2[5] / 2),
                              (x[7][i - 1] + dt * K2[6] / 2), (x[8][i - 1] + dt * K2[7] / 2),
                              (x[9][i - 1] + dt * K2[8] / 2), (x[10][i - 1] + dt * K2[9] / 2),
                              (x[11][i - 1] + dt * K2[10] / 2), (x[12][i - 1] + dt * K2[11] / 2));
            K3[j] = k3;
            double k4 = vf[j]((x[0][i - 1] + dt), (x[1][i - 1] + dt * K3[0]), (x[2][i - 1] + dt * K3[1]),
                              (x[3][i - 1] + dt * K3[2]), (x[4][i - 1] + dt * K3[3]),
                              (x[5][i - 1] + dt * K3[4]), (x[6][i - 1] + dt * K3[5]),
                              (x[7][i - 1] + dt * K3[6]), (x[8][i - 1] + dt * K3[7]),
                              (x[9][i - 1] + dt * K3[8]), (x[10][i - 1] + dt * K3[9]),
                              (x[11][i - 1] + dt * K3[10]), (x[12][i - 1] + dt * K3[11]));
            K4[j] = k4;
        }
        x[1][i] = x[1][i - 1] + dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
        x[2][i] = x[2][i - 1] + dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
        x[3][i] = x[3][i - 1] + dt / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
        x[4][i] = x[4][i - 1] + dt / 6 * (K1[3] + 2 * K2[3] + 2 * K3[3] + K4[3]);
        x[5][i] = x[5][i - 1] + dt / 6 * (K1[4] + 2 * K2[4] + 2 * K3[4] + K4[4]);
        x[6][i] = x[6][i - 1] + dt / 6 * (K1[5] + 2 * K2[5] + 2 * K3[5] + K4[5]);
        x[7][i] = x[7][i - 1] + dt / 6 * (K1[6] + 2 * K2[6] + 2 * K3[6] + K4[5]);
        x[8][i] = x[8][i - 1] + dt / 6 * (K1[7] + 2 * K2[7] + 2 * K3[7] + K4[7]);
        x[9][i] = x[9][i - 1] + dt / 6 * (K1[8] + 2 * K2[8] + 2 * K3[8] + K4[8]);
        x[10][i] = x[10][i - 1] + dt / 6 * (K1[9] + 2 * K2[9] + 2 * K3[9] + K4[9]);
        x[11][i] = x[11][i - 1] + dt / 6 * (K1[10] + 2 * K2[10] + 2 * K3[10] + K4[10]);
        x[12][i] = x[12][i - 1] + dt / 6 * (K1[11] + 2 * K2[11] + 2 * K3[11] + K4[11]);
    }
    return x;
}
