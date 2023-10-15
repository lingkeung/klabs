#include "eig.h"
#include "eigen.h"

void eigs(Matrix A)
{
    cout << "Analysis of Ax = \u03BBx" << endl
         << endl;
    Eigen e(A);
    if (e.getEigIsReal())
    {
        cout << "The eigen system is real." << endl
             << endl;
        if (e.get_isDiag())
        {
            cout << "Matrix is diagonizable." << endl
                 << endl;
            cout << "Eigen space S = " << endl;
            e.getrS().print();
            cout << "Eigen values Lamda = " << endl;
            e.getrLamda().print();
            cout << "Inverse of S, Sinv = " << endl;
            e.getrSinv().print();
        }
        else
        {
            cout << "Matrix is not diagonizable." << endl;
            cout << "Eigen space S = " << endl;
            e.getrS().print();
            cout << "Eigen values Lamda = " << endl;
            e.getrLamda().print();
        }
    }
    else
    {
        cout << "The eigen system is complex" << endl;
    }
}

int main()
{
    //Matrix A(3, 3, {1, 1.5, 0.5, 1.5, 7, 5.5, 0.5, 5.5, 3});
    // Matrix A(3, 3, {3, -1, 0, -1, 3, -1, 0, -1, 3});
    // Matrix A(3, 3, {0, -2, 2, -2, -3, 4, 2, 4, -3});
    // Matrix A(3,3,{2,2,0,2,1,2,0,2,0});
    // Matrix A(3, 3, {0, -1, 0, -1, 0, -1, 0, -1, 0});
    // Matrix A(3, 4, {1, -5, 2, -3, 5, 3, 6, -1, 2, 4, 2, 1});
    // Matrix A(4, 5, {1, 0, 1, -1, -3, 1, 2, -1, 0, -1, 4, 6, -2, -4, 3, 2, -2, 4, -7, 4});
    Matrix A(3,3,{-1,1,0,-4,3,0,1,0,2});
    eigs(A);
    return 0;
}
