#include <iostream>
#include <complex>
#include <vector>
#include "newcmatrix.h"
#include "newcplu.h"
#include "neweig.h"
#include "neweigen.h"

using namespace std;
using namespace kling;

typedef complex<double> C;

int main()
{
    Matrix A(3, 3, vector<double>{1, 2, 3, 5, 3, 7, 3, 8, 5});
    Eigen a(A);
    if (a.getEigIsReal())
    {
        auto rLamda = a.getrLamda();
        auto rS = a.getrS();
        auto rSinv = a.getrSinv();

        (rS * rLamda * rSinv).print();
    }
    else
    {
        auto Lamda = a.getLamda();
        auto S = a.getS();
        auto Sinv = a.getSinv();

        (S * Lamda * Sinv).print();
    }
    return 0;
}