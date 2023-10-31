// #include "svd.h"
#include "svdc.h"
#include "cn.h"
#include "cmatrix.h"

/* void svds(Matrix A)
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
    Matrix x(8, 1, {0, 1, 2, 3, 4, 5, 6, 7});
    // Matrix xx(4, 1, {8, 4, 8, 0});
    // Matrix matlab(8, 1, {1, 3, 5, 7, 9, 11, 13, 15});
    cMatrix result = fdft(x);
    result.print();
    Matrix result1 = nidft(result);
    result1.print();
    Matrix result2 = fidft(result);
    result2.print();
    // cMatrix result2 = ndft(x);
    // result2.print();
    // cMatrix result3 = fdft(xx);
    // result3.print();
    // cMatrix result4 = fdft(matlab);
    // result4.print();
    return 0;
}
