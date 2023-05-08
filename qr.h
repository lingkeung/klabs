#ifndef QR_H
#define QR_H

#include "matrix.h"

void house(Matrix x, Matrix &v, double &beta);
void qr(Matrix A, Matrix &Q, Matrix &R, Matrix &Q1, Matrix &R1);
Matrix qrAxb(Matrix A, Matrix b);
Matrix qrBasis(Matrix A);
void qrGivens(Matrix A, Matrix &Q, Matrix &R, Matrix &Q1, Matrix &R1);

class Givens /* Note 1 */
{
public:
    double c;
    double s;
    Matrix rot;

public:
    Givens(){};
    Givens(double c, double s);
};

Givens givens(double a, double b);

void grm(Givens G, int i, int k, Matrix &A);

void glm(Givens G, int i, int k, Matrix &A);

/* Note 1 : The follow codes demonstrate the use Givens and related functions to zero one of the elements of an order 2 column or row vector.

    QUOTE

    Matrix A(2, 1, vector<double>{3, 4}); //columen vector
    Givens g = givens(3, 4); // construct Givens object g.s = sin(theta)
    glm(g, 1, 2, A); // glm() rotates vector by theta clockwise (left multiply by g.rot)
    cout << "A = "; 
    A.print();  // check below the second component is zeroed

    Matrix B(1, 2, vector<double>{3, 4}); // row vector
    g = givens(3, 4); // construct Givens object g.s = sin(theta)
    grm(g, 1, 2, B);// grm() rotates vector theta clockwise (right multiply by g.rotT)
    cout << "B = "; 
    B.print();  // check below the second component is zeroed

    Matrix C(2, 1, vector<double>{3, 4}); // column vector
    g = givens(4, 3); // construct Givens object g.s = sin(phi( = pi/2 - theta))
    g.s = -1 * g.s; // change sign of phi
    glm(g, 1, 2, C); // glm rotates vector by phi clockwise
    cout << "C = "; // overall vector is rotated anticlockwise by phi (left multiply by g.rotT)
    C.print(); // check below the first component is zeroed

    Matrix D(1, 2, vector<double>{3, 4}); // row vector
    g = givens(4, 3); //construct Givens object g.s = sin(phi( = pi/2 - theta))
    g.s = -1 * g.s; //change sign of phi
    grm(g, 1, 2, D); // grm rotates vector by phi clockwise
    cout << "D = "; // overall vector is rotated anticlockwise by phi (right multiply by g.rot)
    D.print(); // check below the first component is zeroed.

    UNQUOTE
    
    The above codes produce the following result.
   
    A = 
        5.00
        0.00
    B = 
        5.00        0.00
    C = 
        0.00
        5.00
    D = 
        0.00        5.00

    Although not recommended the same effects can be achieved by using direct matrix muliplication. For example the case of Matric C above, as indicated, can be achieved with the codes below.

    QUOTE  

    Matrix C(2, 1, vector<double>{3, 4});
    Givens g = givens(4, 3);  // returns angle between second component and vector
    Matrix result = (g.rot).transpose() * C; // anticlockwise rotation of vector
    result.print();  // first component is zeroed

    UNQUOTE

    The result is as follows:

    C = 
        0.00
        5.00    
     */

#endif