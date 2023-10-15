#include "matrix.h"
#include "plu.h"
#include "qrc.h"
#include "qr.h"
#include "svd.h"
#include "svdc.h"
#include "leqs.h"

/*bool hsoln(Matrix A, Matrix b);
bool usoln(Matrix A);
void legs(Matrix A, Matrix b);
Matrix nullspace(Matrix A);
Matrix xp(Matrix A, Matrix b);*/

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

void leqs(Matrix A, Matrix b)
{
    int n = A.getN();
    int m = A.getM();
    bool e, u;
    u = usoln(A);    // ~ exact or least squares solution is unique if existent
    e = hsoln(A, b); // ~ exact solution exists

    if (u)
    {
        QR q(A);
        Matrix result = q.axb(b);

        if (e && m == n)
        {
            cout << "Unique exact solution = " << endl;
            result.print(4, 12);
        }
        else if (!e && m > n)
        {
            cout << "Unique least squares solution = " << endl;
            result.print(4, 12);
        }
    }
    else
    {
        SVD s(A);
        Matrix xpart = s.pinv() * b;
        // Matrix xpart = xp(A,b);
        Matrix xnull = nullspace(A);
        // Matrix xnull = snull(A);      // alternate solution

        if (e)
        {
            cout << "Exact particular solution = " << endl;
            xpart.print(4, 12);
            cout << "and" << endl
                 << endl;
            cout << "Exact nullspace = " << endl;
            xnull.print(4, 12);
        }
        else if (!e)
        {
            cout << "Minimum norm least squares solution = " << endl;
            xpart.print(4, 12);
        }
    }
}

Matrix nullspace(Matrix A) // find nullspace of underdetermined system
{
    int n = A.getN();
    int m = A.getM();
    vector<int> ivar = icols(A); // basic variables
    vector<int> ieqn;            // independent equations of homogenous system
    if (m != n)
    {
        ieqn = icols(A.transpose());
    }
    else
    {
        ieqn = ivar;
    }

    int r = ivar.size(); // rank(A);
    Matrix R(r, r);      // skeleton coefficient matrix of reduced system
    // R.print();      // debug
    for (int i = 1; i <= r; i++)
    {
        for (int j = 1; j <= r; j++)
        {
            R(i, j) = A(ieqn[i - 1], ivar[j - 1]); // coefficients of basic variables of independent eqn's
        }
    }
    // R.print();                   // debug
    vector<int> fvar;            // free variables
    for (int i = 1; i <= n; i++) // identify free variables
    {
        bool isIn = false;
        for (int j = 0; j < r; j++)
        {
            if (i == ivar[j])
            {
                isIn = true;
                break;
            }
        }
        if (!isIn)
        {
            fvar.push_back(i);
            // cout << "i = " << i << endl; // debug
        }
    }

    Matrix B(r, n - r); // constants vectors for reduced system
                        // formed by assigning 1 to one free variable, rest 0's.

    for (int j = 1; j <= (n - r); j++)
    {
        for (int i = 1; i <= r; i++)

        {
            B(i, j) = -1 * A(ieqn[i - 1], fvar[j - 1]); // negative of coefficient of "1" free variable
        }
    }
    // B.print(); // debug
    QR q(R); // solve reduced system using QR decomposition
    for (int j = 1; j <= n - r; j++)
    {
        B.setblk(1, j, q.axb(B(1, r, j, j))); // overwrite B with solution of linear system
    }
    // B.print();                         // debug
    Matrix nspace(n, n - r);           // construct nullspace
    for (int j = 1; j <= (n - r); j++) // choose one free variable as 1, rest 0's
    {
        nspace(fvar[j - 1], j) = 1;
    }
    // nspace.print();                    // debug
    for (int j = 1; j <= (n - r); j++) // place elements of solution of reduced system in correct position
    {
        for (int k = 0; k < r; k++)
        {
            nspace(ivar[k], j) = B(k + 1, j); // B contains solutions of reduced system
        }
    }
    // nspace.print(); // debug
    return nspace;
}

Matrix xp(Matrix A, Matrix b) // find particular solution of underdetermined system

{
    vector<int> ivar = icols(A);             // basic variables
    vector<int> ieqn = icols(A.transpose()); // independent equations
    int r = ivar.size();                     // rank(A);
    Matrix R(r, r);                          // skeleton coefficient matrix of reduced system
    // R.print();                            // debug
    for (int i = 1; i <= r; i++)
    {
        for (int j = 1; j <= r; j++)
        {
            R(i, j) = A(ieqn[i - 1], ivar[j - 1]); // coefficients of basic variables of independent eqn's
        }
    }
    // R.print();        // debug
    vector<int> fvar; // free variables
    int n = A.getN();
    for (int i = 1; i <= n; i++) // identify free variables
    {
        bool isIn = false;
        for (int j = 0; j < r; j++)
        {
            if (i == ivar[j])
            {
                isIn = true;
                break;
            }
        }
        if (!isIn)
        {
            fvar.push_back(i);
            // cout << "i = " << i << endl; // debug
        }
    }

    Matrix br(r, 1);
    for (int i = 1; i <= r; i++)
    {
        br(i) = b(ieqn[i - 1]);
    }
    // br.print(); // debug
    QR q(R); // solve reduced system using QR decomposition
    br = q.axb(br);
    // br.print(); // debug
    Matrix xpart(n, 1);

    for (int k = 0; k < r; k++)
    {
        xpart(ivar[k]) = br(k + 1); // br contains solutions of reduced system
    }

    // xpart.print(); // debug

    return xpart;
}

Matrix snull(Matrix A)
{
    SVD s(A);
    return s.null();
}
