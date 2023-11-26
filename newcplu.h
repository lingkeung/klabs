#ifndef NEWCPLU_H
#define NEWCPLU_H

#include "newcmatrix.h"

using namespace kling;

cMatrix fsub(cMatrix L, cMatrix b); 

cMatrix bsub(cMatrix U, cMatrix b);

void rexch(int i, int j, cMatrix &A);

/*double abs(C cn);*/

bool cPlu(cMatrix A, cMatrix &P, cMatrix &L, cMatrix &U, C &det);

cMatrix cAxb(cMatrix A, cMatrix b);

C det(cMatrix A);

cMatrix cInverse(cMatrix A);

cMatrix cNull(cMatrix A);

#endif