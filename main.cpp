#include "matrix.h"
#include "plu.h"
#include "qrc.h"
#include "qr.h"
#include "svd.h"
#include "svdc.h"


int main()
{
    Matrix A(4,5, vector<double>{1,0,1,-1,-3,1,2,-1,0,-1,4,6,-2,-4,3,2,-2,4,-7,4});
    Matrix b(4,1,vector<double>{-2,1,7,1});
    A.print();
    b.print();
    SVD s(A);
    Matrix xp = s.pinv()*b;
    xp.print();
    Matrix n = s.null();
    n.print();
    Matrix a(2,1,vector<double>{1.6,26});
    Matrix x = xp + n*a;
    (A*x).print();
    return 0;
}
