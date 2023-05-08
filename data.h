#ifndef DATA_H
#define DATA_H

#include "matrix.h"
#include <string>
#include <vector>
#include <utility> //for std::make_pair()
#include <fstream>
#include <map>

void savePlus(string fileName, Matrix D, vector<string> a)
{

	int col = 20;
	int m = D.getM(), n = D.getN();
	int size = m * n;
	std::ofstream out(fileName);
	out << m << " " << n << endl;
	for (int i = 0; i <= n - 1; i++)
	{
		out << a[i];
	}
	out << endl;
	vector<double> vdbl = D.getVdbl();
	col = (n < col) ? n : col;
	for (int i = 0; i <= size - 1; i++)
	{
		if (i % col == 0)
			out << endl;
		out << vdbl[i] << " ";
	}
}

void loadPlus(string fileName, Matrix &D, vector<string> &a, int &m, int &n)
{

	try
	{
		std::ifstream is(fileName);
		if (is.is_open() == 0)
		{
			throw 0;
		}
		is >> m >> n;
		string s;
		for (int i = 1; i <= n; i++)
		{
			is >> s;
			a.push_back(s);
		}
		vector<double> vdbl;
		double element;
		for (int i = 1; i <= m * n; i++)
		{
			is >> element;
			vdbl.push_back(element);
		}
		D = Matrix(m, n, vdbl);
	}

	catch (int e)
	{
		cout << "void loadPlus() cannot open " << fileName << ". Is this file name valid ?" << endl;
	}
}

class Data
{

public:
	int m;
	int n;
	Matrix D;
	vector<string> a;
	map<string, Matrix> mp;
	string fileName;

public:
	Data();
	Data(string s);
};

Data::Data()
{
	;
}

Data::Data(string s)
{
	fileName = s;
	loadPlus(fileName, D, a, m, n);
	Matrix mtemp(m, 1);
	std::pair<string, Matrix> ptemp;
	for (int i = 1; i <= n; i++)
	{
		mtemp.setblk(1, 1, D(1, m, i, i));
		ptemp = std::make_pair(a[i - 1], mtemp);
		mp.insert(ptemp);
	}
}

#endif