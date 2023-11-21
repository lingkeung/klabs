#include <iostream>
#include <complex>
#include <vector>
#include "newcmatrix.h"

using namespace std;

int main()
{
    newcMatrix mat(2, 2);
    mat.print();
    cout << endl;
    newcMatrix mat1(2, 2, vector<complex<double>>{{1, 2.001}, {3, 4.53}, {3, 5}, {5, 7}});
    mat1.print();
    cout << endl;
    mat1(1, 1) = complex<double>{2, 3};
    cout << "mat1 = " << endl;
    mat1.print();
    cout << endl;
    newcMatrix v(3, 1, vector<complex<double>>{{1, 2}, {3, 4}, {5, 6}});
    cout << v(2) << endl;
    v(1) = complex<double>{9, 9};
    v.print();
    cout << endl;
    mat1.getblk(1, 1, 1, 2).print();
    mat.setblk(1, 1, mat1);
    mat.print();
    cout << endl;
    cout << "row 1 of mat = " << endl;
    mat(1, 1, 1, 2).print();
    (mat + mat1).print();
    (mat - mat1).print();
    (complex<double>{2, 1} * mat).print();
    (mat * mat1).print();
    cout << "3 * mat = " << endl;
    (3*mat).print();
    return 0;
}