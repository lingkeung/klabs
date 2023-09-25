#include "leqs.h"
//#include <iostream>

int main()
{
    //Matrix A(4, 5, vector<double>{1, 0, 1, -1, -3, 1, 2, -1, 0, -1, 4, 6, -2, -4, 3, 2, -2, 4, -7, 4});
    //Matrix b(4, 1, vector<double>{-2, 1, 7, 1});
    // Matrix b(4,1);
    //Matrix A(4,4,vector<double>{1,0,-3,-6,2,-5,1,1,0,-1,2,2,1,-7,4,6});
    //Matrix b(4,1,vector<double>{9,8,-5,0});
    //Matrix A(3, 2, vector<double>{1, 1.2, 1, 1.4, 1, 1.8});
    //Matrix b(3, 1, vector<double>{1.6, 1.7, 2});
    Matrix A(4,4,vector<double>{1,2,3,4,2,2,1,1,2,4,6,8,4,4,2,2});
    Matrix b(4,1,vector<double>{1,2,3,4});
    //Matrix Ab(4,5,vector<double>{1,2,3,4,1,2,2,1,1,2,2,4,6,8,3,4,4,2,2,4});
    //cout << icols(A).size() << endl;
    //cout << icols(Ab).size() << endl;
    //Matrix b(4,1,vector<double>{1,2,3,4});
    legs(A,b);
    return 0;
}
