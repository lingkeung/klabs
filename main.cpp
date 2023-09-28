#include "leqs.h"
// #include <iostream>
vector<int> ieqns(Matrix A)
{
    Matrix AT = A.transpose();
    return icols(AT);
}

int main()
{
    Matrix A(4, 5, vector<double>{1, 0, 1, -1, -3, 1, 2, -1, 0, -1, 4, 6, -2, -4, 3, 2, -2, 4, -7, 4});
    Matrix b(4, 1, vector<double>{-2, 1, 7, 1});
    //  Matrix b(4,1);
    // Matrix A(4,4,vector<double>{1,0,-3,-6,2,-5,1,1,0,-1,2,2,1,-7,4,6});
    // Matrix b(4,1,vector<double>{9,8,-5,0});
    // Matrix A(3, 2, vector<double>{1, 1.2, 1, 1.4, 1, 1.8});
    // Matrix b(3, 1, vector<double>{1.6, 1.7, 2});
    // Matrix A(4, 4, vector<double>{1, 2, 3, 4, 2, 2, 1, 1, 2, 4, 6, 8, 4, 4, 2, 2});
    // Matrix b(4, 1, vector<double>{1, 2, 3, 4});
    // Matrix Ab(4,5,vector<double>{1,2,3,4,1,2,2,1,1,2,2,4,6,8,3,4,4,2,2,4});
    // cout << icols2(A).size() << endl;
    // cout << icols(A).size() << endl;
    // Matrix b(4,1,vector<double>{1,2,3,4});
    // legs(A, b);
    A.print();
    vector<int> ivar = icols(A); // basic variables
    vector<int> ieqn = ieqns(A); // independent equations of homogenous system
    int r = ivar.size();         // rank(A);
    Matrix R(r, r);              // to form coefficient matrix of reduced linear system
    R.print();                   // debug
    for (int i = 1; i <= r; i++)
    {
        for (int j = 1; j <= r; j++)
        {
            R(i, j) = A(ieqn[i - 1], ivar[j - 1]);
        }
    }
    R.print(); // debug

    vector<int> fvar; // find free variables; to be improved using STL search algorithm
    int n = A.getN();
    for (int i = 1; i <= n; i++)
    {
        bool isIn = false;
        for (int j = 0; j < r; j++)
        {
            if (i == ivar[j])
            {
                isIn = true;
            }
            continue;
        }
        if (!isIn)
        {
            fvar.push_back(i);
            cout << "i = " << i << endl;
        }
    }

    Matrix B(r, n - r);

    for (int j = 1; j <= (n - r); j++)
    {
        for (int i = 1; i <= r; i++)

        {
            B(i, j) = -1 * A(ieqn[i - 1], fvar[j - 1]);
        }
    }
    B.print();
    QR q(R);
    for (int j = 1; j <= n - r; j++)
    {
        B.setblk(1, j, q.axb(B(1, r, j, j)));
    }
    B.print();
    Matrix nspace(n, n - r);
    for (int j = 1; j <= (n - r); j++)
    {
        for (int i = 1; i <= n; i++)
        {
            for (int k = 0; k < (n - r); k++)
            {
                if (i == fvar[k])
                {
                    nspace(i, j) = 1;
                    fvar[k] = 0;
                    i = n + 1;
                }
            }
        }
    }

    nspace.print();
    for (int j = 1; j <= (n - r); j++)
    {
        for (int k = 0; k < r; k++)
        {
            nspace(ivar[k], j) = B(k + 1, j);
        }
    }
    nspace.print();
    Matrix randv(2, 1, vector<double>{1, 2});
    (A * (nspace * randv)).print();

    return 0;
}
