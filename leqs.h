#ifndef LEQS_H
#define LEQS_H
#include "matrix.h"
#include "plu.h"
#include "qrc.h"
#include "qr.h"
#include "svd.h"
#include "svdc.h"

bool hsoln(Matrix A, Matrix b);
bool usoln(Matrix A);
void leqs(Matrix A, Matrix b);
Matrix nullspace(Matrix A);
Matrix xp(Matrix A, Matrix b);
Matrix snull(Matrix A);

#endif
