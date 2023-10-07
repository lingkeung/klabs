#include "eig.h"
#include "eigen.h"
#include "svd.h"
#include "svdc.h"
#include "qr.h"
#include "qrc.h"

// #include "cmatrix.h"
//  #include <iostream>

void francisOS(Matrix &H)
{
    int n = H.getN();
    int m = n - 1;
    double s, t, x, y, z;
    s = H(m, m) + H(n, n);
    t = H(m, m) * H(n, n) - H(m, n) * H(n, m);
    x = H(1, 1) * H(1, 1) + H(1, 2) * H(2, 1) - s * H(1, 1) + t + 0.00001;
    y = H(2, 1) * (H(1, 1) + H(2, 2) - s) + 0.00001;
    z = H(2, 1) * H(3, 2);
    cout << x << " " << y << " " << z << endl; // debug
    for (int k = 0; k <= n - 3; k++)
    {
        Matrix xyz(3, 1, vector<double>{x, y, z});
        Matrix v;
        double beta;
        house(xyz, v, beta);
        cout << "xyz =" << endl;           // debug
        xyz.print();                       // debug
        cout << "v = " << endl;            // debug
        v.print();                         // debug
        cout << "beta = " << beta << endl; // debug
        int q = (k > 1) ? k : 1;           // k>=1 changed to k>1 // debug
        H.setblk(k + 1, q, (identity(3) - beta * v * v.transpose()) * H(k + 1, k + 3, q, n));
        int r = (k + 4 <= n) ? k + 4 : n;
        H.setblk(1, k + 1, H(1, r, k + 1, k + 3) * (identity(3) - beta * v * v.transpose()));
        x = H(k + 2, k + 1);
        y = H(k + 3, k + 1);
        if (k < n - 3)
        {
            z = H(k + 4, k + 1);
        }
    }

    Matrix xy(2, 1, vector<double>{x, y});
    Matrix v;
    double beta;
    house(xy, v, beta);
    H.setblk(n - 1, n - 2, (identity(2) - beta * v * v.transpose()) * H(n - 1, n, n - 2, n));
    H.setblk(1, n - 1, H(1, n, n - 1, n) * (identity(2) - beta * v * v.transpose()));
    // return H;
}

int main()
{
    // Matrix A(3, 3, {1,1.5,0.5,1.5,7,5.5,0.5,5.5,3});
    Matrix A(3, 3, {3.0, -1, 0.00000, -1, 3, -1, 0.00000, -1, 3.0});
    Matrix B = A;
    Eigen e(A);
    cMatrix Lamda = e.getLamda();
    Lamda.print();
    SVD s(B);
    Matrix Sigma = s.getSigma();
    Matrix U = s.getU();
    Matrix VT = s.getVT();
    (U * Sigma * VT).print();
    Sigma.print();

    return 0;
}
