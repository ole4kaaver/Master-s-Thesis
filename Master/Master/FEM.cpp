#include <vector>
#include <fstream>
#include <iostream>
#include <set>
using namespace std;
#include "Objects.cpp"

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

	double K, F; // структурная проницаемость, структурная пористость
	

	void Portrait(Init Object)
	{
		ia.resize(Object.elements.size() + 2);
		ja.resize(Object.elements.size());
		ia[0] = 0;
		for (int i = 0; i < Object.elements.size(); i++)
		{
			ia[i + 1] = i;
			ja[i] = i;
		}
		ia[Object.elements.size() + 1] = Object.elements.size();
		elementConnections.resize(Object.elements.size());
		for (int i = 0; i < Object.elements.size(); i++)
		{
			elementConnections[i].resize(2);
			elementConnections[i][0] = i;
			elementConnections[i][1] = i + 1;
		}
	}

	void BuildingGLocal(int number, double h)
	{
		for (int ii = 0; ii < 2; ii++)
			matrixG[ii].resize(2);
		/*double multiplier = K * (k[0][number] / nu_water[number] + k[1][number] / nu_oil) / (1.0 * h);
		G[0][0] = multiplier;
		G[0][1] = -multiplier;
		G[1][0] = -multiplier;
		G[1][1] = multiplier;*/
	}

	FEM(Init Object)
	{
		Portrait(Object);
		matrixG.resize(2);
		al.resize(ia[ia.size() - 1], 0);
		au.resize(ia[ia.size() - 1], 0);
		di.resize(ia.size() - 1);
		double h;
		for (int itemNumber = 0; itemNumber < Object.elements.size(); itemNumber++) // по каждому элементу
		{
			h = Object.elements[itemNumber].xEnd - Object.elements[itemNumber].xBegin;
			BuildingGLocal(itemNumber, h);

			for (int i = 0; i < 2; i++)
			{
				int i_i = elementConnections[itemNumber][i];
				di[itemNumber + i] += matrixG[i][i];
				for (int j = 0; j < i; j++)
				{
					int j_j = elementConnections[itemNumber][j];
					if (i_i < j_j)
						swap(i_i, j_j);
					for (itemNumber = ia[i_i]; ja[itemNumber] - j_j != 0; itemNumber++);
						al[itemNumber] += matrixG[i][j]; // для глобальной матрицы A
				}
			}
		}
	}

	FEM()
	{

	}
};