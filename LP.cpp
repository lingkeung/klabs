#include "lp.h"
#include "simplex.h"
#include "matrix.h"
#include <vector>
#include "iostream"
#include <algorithm>
#include <deque>

using namespace std;

LP::LP(Matrix A, Matrix b, Matrix c, Matrix s)
{

	this->A = A;
	this->b = b;
	this->c = c;
	this->s = s;
	m = A.getM();
	n = A.getN();
	lps(A, b, c, s, x, minCost, status, sP, t);
}

int LP::branch(int i, LP &lpf, LP &lpc)
{
	if (status == 0)
	{
		double flr = floor(x(i));
		double cel = ceil(x(i));
		Matrix c1 = c;
		Matrix b1(m + 1, 1);
		Matrix s1(m + 1, 1);
		Matrix A1(m + 1, n);
		b1.setblk(1, 1, b);
		b1(m + 1) = flr;
		s1.setblk(1, 1, s);
		s1(m + 1) = 1;
		A1.setblk(1, 1, A);
		A1(m + 1, i) = 1;
		lpf = LP(A1, b1, c1, s1);
		Matrix b2 = b1;
		Matrix s2 = s1;
		b2(m + 1) = cel;
		s2(m + 1) = -1;
		Matrix A2 = A1;
		Matrix c2 = c1;
		lpc = LP(A2, b2, c2, s2);
		return 0;
	}
	else
		return 1;
}

vector<LP> LP::branch(vector<int> idx)
{

	int size = idx.size();
	vector<LP> result;
	LP lpf, lpc;
	int i;
	for (int k = 0; k < size; k++)
	{
		i = idx[k];
		if (!(abs(floor(x(i)) - x(i)) < 1e-6 || abs(ceil(x(i)) - x(i)) < 1e-6))
		{
			branch(i, lpf, lpc);
			result.push_back(lpf);
			result.push_back(lpc);
		}
	}

	return result;
}

bool LP::isInt(vector<int> idx)
{
	int size = idx.size();
	bool bol = false;
	for (int k = 0; k < size; k++)
	{
		int i = idx[k];
		if (abs(floor(x(i)) - x(i)) < 1e-6 || abs(ceil(x(i)) - x(i)) < 1e-6)
		{
			bol = true;
		}
		else
		{
			bol = false;
			break;
		}
	}
	return bol;
}

void LP::print(int precision, int width)
{
	if (status == 0)
	{
		cout << endl;
		cout << "KLABS LINEAR PROGRAMMING REPORT" << endl;
		cout << "===============================" << endl;
		cout << "\nMINIMIZE : " << endl
			 << endl;
		;
		cout << "    ";
		int k = 0;
		for (k = 1; k <= n; k++)
		{
			if (abs(c(k)) > 1e-6)
			{
				cout << c(k) << " X" << k;
				break;
			}
		}
		for (int i = k + 1; i <= n; i++)
		{
			if (abs(c(i)) > 1e-6)
			{
				const char *sign = (abs(c(i)) == c(i) ? "  +  " : "  -  ");
				cout << sign << abs(c(i)) << " X" << i;
			}
		}
		cout << "  = Z,     OR MAXIMIZE -Z" << endl;

		cout << endl
			 << endl;
		cout << "SUBJECT TO :" << endl
			 << endl;

		for (int i = 1; i <= m; i++)
		{
			cout << setw(4) << i << ")  ";
			int p = 0;
			bool zero = true;
			for (p = 1; p <= n; p++)
			{
				if (abs(A(i, p)) > 1e-6)
				{
					zero = false;
					if (abs(A(i, p)) != 1)
						cout << A(i, p) << " X" << p;
					else
						cout << " X" << p;
					break;
				}
			}
			if (zero)
				cout << 0;
			else
			{
				for (int j = p + 1; j <= n; j++)
				{
					if (abs(A(i, j)) > 1e-6)
					{
						const char *sign = (abs(A(i, j)) == A(i, j) ? "  +  " : "  -  ");
						if (abs(A(i, j)) != 1)
							cout << sign << abs(A(i, j)) << " X" << j;
						else
							cout << sign << " X" << j;
					}
				}
			}

			const char *relation;
			if (s(i) == 0)
				relation = "   =  ";
			else if (s(i) == -1)
				relation = "  >=  ";
			else
				relation = "  <=  ";
			cout << relation << b(i) << endl;
		}
		cout << "\n        Xi >= 0 (i = 1, ... ," << n << ")" << endl
			 << endl
			 << endl;

		if (status != 0)
			cout << "FAILED TO FIND LP OPTIMUM" << endl
				 << endl;
		else
		{
			cout << "LP OPTIMUM FOUND AFTER " << t << " PIVOTS / BRANCHINGS :" << endl
				 << endl;
			;
			cout << "    OBJECTIVE FUNCTION VALUE = " << minCost << endl
				 << endl;
		}
		cout << "    VARIABLE" << setw(12) << "VALUE" << endl;
		for (int i = 1; i <= n; i++)
		{
			if (i > 9)
				cout << "   X" << i << "          ";
			else
				cout << "        X" << i;
			cout << "    " << setw(11) << (abs(x(i)) <= 1e-10 ? 0 : x(i)) << endl;
		}
		cout << endl;
		cout << "    CONSTRAINT       SHADOW PRICE       SHADOW PRICE (maximize)" << endl;
		for (int i = 1; i <= m; i++)
		{
			cout << "    " << setw(5) << i;
			cout << "    " << setw(18) << (abs(sP(i)) <= 1e-10 ? 0 : sP(i)) << setw(18) << (abs(sP(i)) <= 1e-10 ? 0 : -1 * sP(i)) << endl;
		}
		cout << endl;
	}
}

void pLPi(LP lp, LP lpi, vector<int> idx, int iter)
{

	Matrix A = lp.A;
	Matrix b = lp.b;
	Matrix c = lp.c;
	Matrix s = lp.s;
	int m = lp.m;
	int n = lp.n;

	cout << endl;
	cout << "           KLABS LINEAR PROGRAMMING REPORT" << endl;
	cout << "           ===============================" << endl
		 << endl;

	cout << "1. lpi(.) RETURNED Vector<LP> OF FEASIBLE AND OPTIMAL SOLUTIONS" << endl
		 << endl
		 << endl;

	cout << "2. REPORT ON THE LAST FOUND OPTIMAL SOLUTION :" << endl
		 << endl;
	cout << "\nMINIMIZE : " << endl
		 << endl;
	;
	cout << "    ";
	int k = 0;
	for (k = 1; k <= n; k++)
	{
		if (abs(c(k)) > 1e-6)
		{
			cout << c(k) << " X" << k;
			break;
		}
	}
	for (int i = k + 1; i <= n; i++)
	{
		if (abs(c(i)) > 1e-6)
		{
			const char *sign = (abs(c(i)) == c(i) ? "  +  " : "  -  ");
			cout << sign << abs(c(i)) << " X" << i;
		}
	}
	cout << "  =  Z,     OR MAXIMIZE -Z" << endl;

	cout << endl
		 << endl;
	cout << "SUBJECT TO :" << endl
		 << endl;

	for (int i = 1; i <= m; i++)
	{
		cout << setw(4) << i << ")  ";
		int p = 0;
		for (p = 1; p <= n; p++)
		{
			if (abs(A(i, p)) > 1e-6)
			{
				if (abs(A(i, p)) != 1)
					cout << A(i, p) << " X" << p;
				else
					cout << " X" << p;
				break;
			}
		}
		for (int j = p + 1; j <= n; j++)
		{
			if (abs(A(i, j)) > 1e-6)
			{
				const char *sign = (abs(A(i, j)) == A(i, j) ? "  +  " : "  -  ");
				if (abs(A(i, j)) != 1)
					cout << sign << abs(A(i, j)) << " X" << j;
				else
					cout << sign << " X" << j;
			}
		}
		const char *relation;
		if (s(i) == 0)
			relation = "   =  ";
		else if (s(i) == -1)
			relation = "  >=  ";
		else
			relation = "  <=  ";
		cout << relation << b(i) << endl;
	}
	cout << "\n        Xi >= 0 (i = 1, ... ," << n << ")" << endl;
	cout << "        Integer Constraint : ";
	for_each(idx.begin(), idx.end(), [](int e)
			 { cout << e << ","; });
	cout << endl
		 << endl
		 << endl;

	if (lp.status || lpi.status)
		cout << "FAILED TO FIND LP OPTIMUM" << endl
			 << endl;
	else
	{
		cout << "WITHOUT INTEGER CONSTRAINT : " << endl
			 << endl;
		cout << "    LP OPTIMUM FOUND AFTER " << lp.t << " PIVOTS" << endl;
		cout << "    OBJECTIVE FUNCTION VALUE = " << lp.minCost << endl
			 << endl;
		cout << "\nWITH INTEGER CONSTRAINT : " << endl
			 << endl;
		cout << "    LP OPTIMUM FOUND AFTER " << iter << " BRANCHES" << endl;
		cout << "    OBJECTIVE FUNCTION VALUE = " << lpi.minCost << endl
			 << endl;
	}
	Matrix x = lpi.x, sP = lpi.sP;
	cout << "    VARIABLE" << setw(12) << "VALUE" << endl;
	for (int i = 1; i <= n; i++)
	{
		if (i > 9)
			cout << "       X" << i << "          ";
		else
			cout << "        X" << i;
		cout << setw(12) << (abs(x(i)) <= 1e-10 ? 0 : x(i)) << endl;
	}
	cout << endl;
	cout << "    CONSTRAINT" << setw(18) << "SHADOW PRICE" << setw(28) << "SHADOW PRICE (maximize)" << endl;
	for (int i = 1; i <= m; i++)
	{
		cout << "    " << setw(5) << i;
		cout << "    " << setw(15) << (abs(sP(i)) <= 1e-10 ? 0 : sP(i)) << setw(18) << (abs(sP(i)) <= 1e-10 ? 0 : -1 * sP(i)) << endl;
	}
	cout << endl;
}

vector<LP> lpi(Matrix A, Matrix b, Matrix c, Matrix s, vector<int> idx, char chr)
{

	if (idx.empty())
	{
		vector<LP> result;
		LP lp(A, b, c, s);
		result.push_back(lp);
		if (chr != 's')
		{
			cout << "\nNO INTEGER CONSTRAINT, LP RELEXATION RETURNED." << endl;
			lp.print();
		}
		return result;
	}
	LP lp(A, b, c, s);
	vector<LP> result;
	if (lp.isInt(idx))
	{
		result.push_back(lp);
		if (chr != 's')
		{
			cout << "\nLP RELEXATION OPTIMAL IS AN INTEGER SOLUTION !" << endl;
			lp.print();
		}
		return result;
	}
	deque<LP> tbi;
	double ubound = lp.minCost;
	double bound = 1e100; // any large number
	tbi.push_back(lp);
	for (int i = 1; i <= 100; i++)
	{
		LP lp0 = tbi.front();
		tbi.pop_front();
		if (lp0.status == 0)
		{
			vector<LP> bran = lp0.branch(idx);
			int nb = bran.size();
			for (int i = 0; i < nb; i++)
			{
				LP lp1 = bran[i];
				if ((lp1.isInt(idx)) && ((lp1.minCost) <= bound) && ((lp1.minCost) >= ubound) && lp1.status == 0)
				{
					bound = lp1.minCost;
					result.push_back(lp1);
				}
				if ((!lp1.isInt(idx)) && ((lp1.minCost) <= bound) && ((lp1.minCost) >= ubound) && lp1.status == 0)
				{
					tbi.push_back(lp1);
				}
			}
		}
		if (!result.empty() && tbi.empty())
		{
			if (chr != 's')
			{
				pLPi(lp, result.back(), idx, i);
			}
			return result;
		}
		else if (result.empty() && tbi.empty())
		{
			if (chr != 's')
			{
				cout << "No feasible solutions found. Empty vector<LP> returned." << endl;
			}
			return result;
		}
	}

	if (chr != 's')
	{
		cout << "Search for optimum stopped after 100 branchings. " << endl;
		cout << "vector<LP> returned may contain some feasible solutions or it may be empty" << endl;
	}
	return result;
}
