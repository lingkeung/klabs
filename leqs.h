#include "matrix.h"
#include "plu.h"
#include "qrc.h"
#include "qr.h"
#include "svd.h"
#include "svdc.h"

bool hsoln(Matrix A, Matrix b);
bool usoln(Matrix A);
void legs(Matrix A, Matrix b);

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

void legs(Matrix A, Matrix b)
{
    int n = A.getN();
    int m = A.getM();
    bool e, u;
    u = usoln(A);    // ~ exact solution is unique if existent
    e = hsoln(A, b); // ~ exact solution exists

    if (u)
    {
        QR q(A);
        Matrix result = q.axb(b);

        if (e && m == n)
        {
            cout << "Unique exact solution = " << endl;
            result.print();
        }
        else if (!e && m > n)
        {
            cout << "Unique least squares solution = " << endl;
            result.print();
        }
    }
    else
    {
        SVD s(A);
        Matrix xp = s.pinv() * b;
        Matrix null = s.null();

        if (e)
        {
            cout << "Exact particular solution = " << endl;
            xp.print();
            cout << "and" << endl
                 << endl;
            cout << "Exact nullspace = " << endl;
            null.print();
        }
        else if (!e)
        {
            cout << "Minimum norm least squares solution = " << endl;
            xp.print();
        }
    }
}