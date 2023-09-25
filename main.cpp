#include "matrix.h"
#include "plu.h"
#include "qrc.h"
#include "qr.h"
#include "svd.h"
#include "svdc.h"

bool hsoln(Matrix A, Matrix b)
{
    bool result = false;
    int m = A.getM();
    int n = A.getN();
    if (m != b.getM() || b.getN() != 1)
    {
        return result;
    }
    Matrix Ab(m, n + 1);
    Ab.setblk(1, 1, A);
    Ab.setblk(1, n + 1, b);
    if (icols(A).size() == icols(Ab).size())
    {
        result = true;
    }
    return result;
}

bool usoln(Matrix A)
{
    bool result = false;
    if (icols(A).size() == A.getN())
    {
        result = true;
    }
    return result;
}

int main()
{
    // Matrix A(4, 5, vector<double>{1, 0, 1, -1, -3, 1, 2, -1, 0, -1, 4, 6, -2, -4, 3, 2, -2, 4, -7, 4});
    // Matrix b(4, 1, vector<double>{-2, 1, 7, 1});
    // Matrix b(4,1);
    // Matrix A(4,4,vector<double>{1,0,-3,-6,2,-5,1,1,0,-1,2,2,1,-7,4,6});
    // Matrix b(4,1,vector<double>{9,8,-5,0});
    Matrix A(3, 2, vector<double>{1, 1.2, 1, 1.4, 1, 1.8});
    Matrix b(3, 1, vector<double>{1.6, 1.7, 2});
    int n = A.getN();
    int m = A.getM();

    if (usoln(A))
    {
        QR q(A);
        Matrix result = q.axb(b);

        if (hsoln(A, b) && m == n)
        {
            cout << "Unique exact solution = " << endl;
            result.print();
        }
        else if (!hsoln(A, b) && m > n)
        {
            cout << "Unique least square solution = " << endl;
            result.print();
        }
    }
    else
    {
        SVD s(A);
        Matrix xp = s.pinv() * b;
        Matrix null = s.null();

        if (hsoln(A, b))
        {
            cout << "Exact particular solution = " << endl;
            xp.print();
            cout << "and" << endl << endl;
            cout << "Exact nullspace = " << endl;
            null.print();
        }
        else if (!hsoln(A, b) && m > n)
        {
            cout << "Minimum norm least square solution = " << endl;
            xp.print();
        }
    }

    return 0;
}
