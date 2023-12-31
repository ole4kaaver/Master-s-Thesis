#include <vector>
#include <iostream>
#include <set>
using namespace std;
#include "Objects.h"
#define test

class FEM
{
public:
	const int maxIter = 10000;
	const double eps = 1e-10;
	const double eps1 = 1e-8;
	const double epsS = 1e-8;
	vector<int> ia, ja;
	vector<double> al, au, di;
	vector<vector<double>> matrixG;
	vector<vector<double>> elementConnections;
	vector<double> globalB, localB;
	vector<double> y;

#ifdef test
	double Test(double x);
	double F(double x);
#endif

	void Portrait(vector <FiniteElement> elements);

	void BuildingGLocal(int number, Init Object);

	void BuildingBLocal(int number, Init Object);

	void LUDecomposition(int n);

	void BoundaryConditions1(vector <FiniteElement> elem);

	void BoundaryConditions2(double Q0);

	void StraightRun(int n);

	void Reverse(int n);

	FEM(Init Object)
	{
		Portrait(Object.elements);
		matrixG.resize(2);
		for (int ii = 0; ii < 2; ii++)
			matrixG[ii].resize(2, 0.0);
		localB.resize(2);
		al.resize(ia.back(), 0.0);
		au.resize(ia.back(), 0.0);
		di.resize(Object.elements.size() + 1, 0.0);
		globalB.resize(Object.elements.size() + 1, 0.0);
		double h;
		for (int itemNumber = 0; itemNumber < Object.elements.size(); itemNumber++) // �� ������� ��������
		{
			h = Object.elements[itemNumber].xEnd - Object.elements[itemNumber].xBegin;
			BuildingGLocal(itemNumber, Object);
			BuildingBLocal(itemNumber, Object);

			for (int i = 0; i < 2; i++)
			{
				globalB[itemNumber + i] += localB[i];
				int i_i = elementConnections[itemNumber][i];
				di[itemNumber + i] += matrixG[i][i];
				for (int j = 0; j < i; j++)
				{
					int j_j = elementConnections[itemNumber][j];
					if (i_i < j_j)
						swap(i_i, j_j);
					for (itemNumber = ia[i_i]; ja[itemNumber] - j_j != 0; itemNumber++);
					al[itemNumber] += matrixG[i][j]; // ��� ���������� ������� A
				}
			}
		}
		au = al;

		BoundaryConditions2(Object.flow0);
		BoundaryConditions1(Object.elements);
		LUDecomposition(Object.elements.size() + 1);
		StraightRun(Object.elements.size() + 1);
		Reverse(Object.elements.size() + 1);
	}

	FEM()
	{

	}
};