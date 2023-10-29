// #include "svd.h"
#include "svdc.h"
#include "cn.h"
#include "cmatrix.h"



/* bool isSymmetric(Matrix A)
{
    bool result = true;
    for (int i = 1; i <= A.getM(); i++)
    {
        for (int j = 1; j <= A.getN(); j++)
        {
            if (abs(A(i, j) - A(j, i)) >= 1e-10)
            {
                result = false;
                return result;
            }
        }
    }
    cout << result << endl;
    return result;
}

void svds(Matrix A)
{
    cout << "SVD Analysis" << endl;
    cout << "Matrix A = " << endl;
    A.print();
    SVD s(A);
    Matrix U = s.getU();
    Matrix Sigma = s.getSigma();
    Matrix VT = s.getVT();
    cout << "SVD decomposition A = U*Sigma*VT : " << endl;
    cout << "U = " << endl;
    U.print();
    cout << "Sigma = " << endl;
    Sigma.print();
    cout << "VT = " << endl;
    VT.print();
    cout << "The SVD decompostion of A reveals the following :" << endl;
    int rank(double tol = 1e-8);

    Matrix null(double tol = 1e-8);
    Matrix range(double tol = 1e-8);
    Matrix sumRank1(int k1, int k2, double tol = 1e-8);
    Matrix pinv(double tol = 1e-8);
    void polar(Matrix & Q, Matrix & S);
    Matrix axb(Matrix b, bool &c, bool &u, Matrix &P, Matrix &Q, double tol = 1e-8);
    double cond(double tol = 1e-8);
    Matrix inv(double tol = 1e-8);
}
 */
int main()
{
    // Matrix A(3, 3, {1, 1.5, 0.5, 1.5, 7, 5.5, 0.5, 5.5, 3});
    //   Matrix A(3, 3, {3, -1, 0, -1, 3, -1, 0, -1, 3});
    //   Matrix A(3, 3, {0, -2, 2, -2, -3, 4, 2, 4, -3});
    //   Matrix A(3,3,{2,2,0,2,1,2,0,2,0});
    //   Matrix A(3, 3, {0, -1, 0, -1, 0, -1, 0, -1, 0});
    //   Matrix A(3, 4, {1, -5, 2, -3, 5, 3, 6, -1, 2, 4, 2, 1});
    //   Matrix A(4, 5, {1, 0, 1, -1, -3, 1, 2, -1, 0, -1, 4, 6, -2, -4, 3, 2, -2, 4, -7, 4});
    //   Matrix A(3, 3, {-1, 1, 0, -4, 3, 0, 1, 0, 2});
    //  Matrix A(2, 2, {1, -1, 1, 1});
    //  Matrix A(3,3,{4.0/5, -3.0/5,0,3.0/5,4.0/5,0,1,2,2});
    //  Matrix A(3, 3, {0, -1, 1, -1, 0, 1, 1, 1, 0});
    //  Matrix B(3,3,{-2,1,0,-4,2,0,1,0,1});
    //  svds(A);
    //  cout << fixed << setprecision(15);
    /* const double pi = M_PI;
    int N = 4;
    double ampW = 1;
    double argW = -2 * pi / N;
    cMatrix W(N, N);
    for (int n = 0; n <= N - 1; n++)
    {
        for (int k = 0; k <= N - 1; k++)
        {
            W(n + 1, k + 1) = pp2r(ampW, argW, n * k);
        }
    }
    cMatrix x(N, 1, {C(1, 0), C(1, 0), C(0, 0), C(0, 0)});
    cMatrix X = W * x;
    X.print(); */
    Matrix x(4,1,{2,0,-2,0});
    ndft(x).print();

    return 0;
}
