#include "matrix.h"
#include "plu.h"

vector<int> icolssss(Matrix A)
{
    vector<int> result;
    int m = A.getM(), n = A.getN();

    int p = 0;
    for (int k = 1; k <= m - 1; k++)
    {
        double max = 0;
        for (int i = k; i <= m; i++)
        {
            if (abs(A(i, k)) > max)
            {
                max = abs(A(i, k));
                p = i;
            }
        }
        if (p != k)
        {
            rexch(p, k, A);
        }

        if (abs(A(k, k)) >= 1e-6)
        {
            A.setblk(k + 1, k, (1 / A(k, k)) * A(k + 1, m, k, k)); //store multipliers (L)
            A.setblk(k + 1, k + 1, A(k + 1, m, k + 1, n) - A(k + 1, m, k, k) * A(k, k, k + 1, n));
            for (int i = k + 1; i <= m; i++)
            {
                A(i, k) = 0; //restore zeros (erase multipliers)
            }
            result.push_back(k);
        }
    }
    if (abs(A(m, m)) >= 1e-6)
    {
        result.push_back(m); //check for last column pivot
    }
    //A.print(); // for debugging
    return result;
}

int main()
{
    Matrix A(8, 2, vector<double>{1, 0, 1, -1, 1, 2, -1, 0, 4, 6, -2, -4, 2, -2, 4, -7});
    vector<int> idx = icols(A);
    cout << idx[0] << " " << idx[1] << " " << idx[2] << endl;
    Matrix B(3, 5, vector<double>{0,2,-4,-1,-4,5,3,1,7,0,5,-10,2,3,0});
    vector<int> idx2 = icols(B);
    cout << idx2[0] << " " << idx2[1] << " " << idx2[2] << endl;
    return 0;
}
