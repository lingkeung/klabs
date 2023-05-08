#ifndef LP_H
#define LP_H

#include "simplex.h"
#include "matrix.h"
#include <vector>
#include <list>

using namespace std;

class LP
{
public:
	Matrix A;
	Matrix b;
	Matrix c;
	Matrix s;
	Matrix x;
	int m;
	int n;
	int status;
	double minCost;
	Matrix sP;
	int t;

public:
	LP() {}
	LP(Matrix A, Matrix b, Matrix c, Matrix s);
	int branch(int i, LP &lpf, LP &lpc);
	vector<LP> branch(vector<int> idx);
	bool isInt(vector<int> idx);
	void print(int precision = 2, int width = 10);
};

void pLPi(LP lp, LP lpi, vector<int> idx, int iter);

vector<LP> lpi(Matrix A, Matrix b, Matrix c, Matrix s, vector<int> idx, char chr = 's');

#endif