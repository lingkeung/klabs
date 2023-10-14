#include "eig.h"
#include "eigen.h"
#include "svd.h"
#include "svdc.h"
#include "qr.h"
#include "qrc.h"
#include "plu.h"
#include "leqs.h"

int main()
{
    //Matrix A(3, 3, {1, 1.5, 0.5, 1.5, 7, 5.5, 0.5, 5.5, 3});
    Matrix A(3, 3, {3.0, -1, 0.00000, -1, 3, -1, 0.00000, -1, 3.0});
    // Matrix A(3, 3, {0, -2, 2, -2, -3, 4, 2, 4, -3});
    //  Matrix A(3,3,{2,2,0,2,1,2,0,2,0});
     //Matrix A(3, 3, {0, -1, 0, -1, 0, -1, 0, -1, 0});
    Eigen e(A);
    Matrix Lamda = e.getrLamda();
    Matrix S = e.getrS();
    Matrix Sinv = e.getrSinv();
    cout << "Check S*Lamda*Sinv" << endl;
    (S * Lamda * Sinv).print();
    cout << "= A" << endl;
    A.print();
   
    return 0;
}
