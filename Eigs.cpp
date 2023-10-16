#include "eig.h"
#include "eigen.h"
//#include "leqs.h"

void eigs(Matrix A)
{
    cout << "Analysis of Ax = \u03BBx" << endl
         << endl;
    cout << "Matrix A = " << endl;
    A.print();
    Eigen e(A);
    if (e.getEigIsReal())
    {
        cout << "The eigen system is real." << endl
             << endl;
        if (e.get_isDiag())
        {
            cout << "Eigenvalues Lamda = " << endl;
            e.getrLamda().print(4, 12);
            cout << "Matrix is diagonizable." << endl
                 << endl;
            cout << "Eigen space S = " << endl;
            e.getrS().print(4, 12);
            cout << "Inverse of S, Sinv = " << endl;
            e.getrSinv().print(4, 12);
        }
        else
        {
            cout << "Eigenvalues Lamda = " << endl;
            Matrix Lamda = e.getrLamda();
            Lamda.print(4, 12);
            cout << "Matrix is not diagonizable." << endl;
            cout << "Eigen space S = " << endl;
            Matrix S = e.getrS();
            S.print(4, 12);
        }
    }
    else
    {
        cout << "The eigen system is complex" << endl;
        if (e.get_isDiag())
        {
            cout << "Eigenvalues Lamda = " << endl;
            e.getLamda().print();
            cout << "Matrix is diagonizable." << endl
                 << endl;
            cout << "Eigen space S = " << endl;
            e.getS().print();
            cout << "Inverse of S, Sinv = " << endl;
            e.getSinv().print();
        }
        else
        {
            cout << "Eigenvalues Lamda = " << endl;
            cMatrix Lamda = e.getLamda();
            Lamda.print();
            cout << "Matrix is not diagonizable." << endl;
            cout << "Eigen space S = " << endl;
            cMatrix S = e.getS();
            S.print();
        }
    }
}

int main()
{
    // Matrix A(3, 3, {1, 1.5, 0.5, 1.5, 7, 5.5, 0.5, 5.5, 3});
    //  Matrix A(3, 3, {3, -1, 0, -1, 3, -1, 0, -1, 3});
    //  Matrix A(3, 3, {0, -2, 2, -2, -3, 4, 2, 4, -3});
    //  Matrix A(3,3,{2,2,0,2,1,2,0,2,0});
    //  Matrix A(3, 3, {0, -1, 0, -1, 0, -1, 0, -1, 0});
    //  Matrix A(3, 4, {1, -5, 2, -3, 5, 3, 6, -1, 2, 4, 2, 1});
    //  Matrix A(4, 5, {1, 0, 1, -1, -3, 1, 2, -1, 0, -1, 4, 6, -2, -4, 3, 2, -2, 4, -7, 4});
    //  Matrix A(3, 3, {-1, 1, 0, -4, 3, 0, 1, 0, 2});
    //Matrix A(2, 2, {1, -1, 1, 1});
    Matrix A(3,3,{4.0/5, -3.0/5,0,3.0/5,4.0/5,0,1,2,2});
    eigs(A);
    // Matrix B(3,3,{-2,1,0,-4,2,0,1,0,1});
    // Matrix ns = nullspace(B);
    // ns.print();
    //(B*ns).print();
    return 0;
}
