#include "eig.h"
#include "eigen.h"
#include "svd.h"
#include "svdc.h"
#include "qr.h"
#include "qrc.h"
// #include "matrix.h"

int main()
{
    // Matrix A(3, 3, {1, 1.5, 0.5, 1.5, 7, 5.5, 0.5, 5.5, 3});
    //  Matrix A(3, 3, {3.0, -1, 0.00000, -1, 3, -1, 0.00000, -1, 3.0});
    //   Matrix A(3, 3, {0, -2, 2, -2, -3, 4, 2, 4, -3});
    //  Matrix A(3, 3, {2, 2, 0, 2, 1, 2, 0, 2, 0});
    Matrix A(3, 3, vector<double>{1, 2, 3, 5, 3, 7, 3, 8, 5});
    // Matrix A(3, 3, vector<double>{4.0 / 5, -3.0 / 5, 0, 3.0 / 5, 4.0 / 5, 0, 1, 2, 2});
    Eigen a(A);
    if (a.getEigIsReal())
    {
        cout << "Eigen system is real" << endl;
        auto rLamda = a.getrLamda();
        auto rS = a.getrS();
        auto rSinv = a.getrSinv();

        (rS * rLamda * rSinv).print();
    }
    else
    {
        cout << "Eigen system is complex" << endl;
        auto Lamda = a.getLamda();
        auto S = a.getS();
        auto Sinv = a.getSinv();

        (S * Lamda * Sinv).print();
        (S * Lamda * Sinv).print();
    }
    a.getS().print();
    a.getrS().print();

    // Matrix A(3,3,{0,-1,1,-1,0,1,1,1,0});
    /*
   char isPositive;
   cout << "Matrix A is positive semidefinite or definite , y or n ?" << endl;
   cin >> isPositive;
   if (isPositive == 'y')
   {
       SVD s(A);
       Matrix Lamda = s.getSigma();
       cout << "Lamda = " << endl;
       Lamda.print();
       cout << "Eigenspace = " << endl;
       s.getU().print();
       cout << "Eigenspace inverse = " << endl;
       s.getVT().print();
   }
   else
   {
       Eigen e(A);
       cMatrix Lamda = e.getLamda();
       Matrix rLamda;
       if (e.realLamda(rLamda))
       {
           cout << "Lamda is real = " << endl;
           rLamda.print();
           Matrix rS = real(e.getS());
           cout << "Eigenspace is real = " << endl;
           rS = normalize(rS);
           rS.print();

           if (e.get_isDiag())
           {
               QR q(rS);
               Matrix rSinv = q.inv();
               cout << "Matrix can be diagonalized; eigenspace inverse = " << endl;
               rSinv.print();
               //(rS*rLamda*rSinv).print(); // debug
           }
       }
       else
       {
           cout << "Eigenspace = " << endl;
           e.getS().print();
           cout << "Lamda = " << endl;
           Lamda.print();
           if (e.get_isDiag())
           {
               cMatrix Sinv = e.getSinv();
               cout << "Matrix can be diagonalized; eigenspace inverse = " << endl;
               Sinv.print();
           }
       }
   }
*/
    return 0;
}
