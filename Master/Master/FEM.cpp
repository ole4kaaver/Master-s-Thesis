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
	vector<double> globalB;
	vector<double> y;
	double Q0; // начальный поток
	double K, F; // структурная проницаемость, структурная пористость

	double Test(double x)
	{
		return x;
	}	

	void Portrait(vector <FiniteElement> elements)
	{
		ia.resize(elements.size() + 2);
		ja.resize(elements.size());
		ia[0] = 0;
		for (int i = 0; i < elements.size(); i++)
		{
			ia[i + 1] = i;
			ja[i] = i;
		}
		ia[elements.size() + 1] = elements.size();
		elementConnections.resize(elements.size());
		for (int i = 0; i < elements.size(); i++)
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

	void LUDecomposition(int n)
	{
		//cout << endl << endl << "       Start LU_decomposition" << endl;

		int j0 = 0;
		int j = 0;
		int jbeg = 0;
		int jend = 0;
		int j0j = 0;
		int jjbeg = 0;
		int jjend = 0;
		int indau = 0;
		int indal = 0;
		double cu = 0;
		double cl = 0;
		//int ii = 0;

		for (int i = 1; i < n; ++i)
		{
			double sumdi = 0;
			j0 = i - (ia[i + 1] - ia[i]);
			for (int ii = ia[i]; ii < ia[i + 1]; ++ii)
			{
				j = ii - ia[i] + j0;
				jbeg = ia[j];
				jend = ia[j + 1];
				if (jbeg < jend)
				{
					j0j = j - (jend - jbeg);
					jjbeg = max(j0, j0j);
					jjend = min(j, i - 1);
					cl = 0;
					for (int k = 0; k < jjend - jjbeg; k++)
					{
						indau = ia[j] + jjbeg - j0j + k;
						indal = ia[i] + jjbeg - j0 + k;
						cl = cl + au[indau] * al[indal];
					}
					al[ii] = al[ii] - cl;
					cu = 0;
					for (int k = 0; k < jjend - jjbeg; k++)
					{
						indal = ia[j] + jjbeg - j0j + k;
						indau = ia[i] + jjbeg - j0 + k;
						cu = cu + au[indau] * al[indal];
					}
					au[ii] = au[ii] - cu;
				}
				au[ii] = au[ii] / di[j];
				sumdi += al[ii] * au[ii];
			}
			di[i] -= sumdi;
		}
	}

	void BoundaryConditions1(vector <FiniteElement> elem) // стр 237
	{
		double B = 1e+10, sum = 0;
		double x1, x2;
		x1 = elem[0].xBegin;
		x2 = elem[elem.size() - 1].xEnd;
		double res1 = Test(x1) * 101325.0;
		double res2 = Test(x2) * 101325.0;
		/*double res1 = test(x1);
		double res2 = test(x2);*/
		di[elem.size()] = B;
		globalB[elem.size()] = B * res2;
		/*di[0] = B;
		global_b[0] = B * res1;*/
	}

	void BoundaryConditions2()
	{
		globalB[0] += Q0;
		//global_b[0] -= K * (k[0][0] / nu_water[0] + k[1][0] / nu_oil);
		//global_b[0] -= K * (k[0][0] / nu_water[0] + k[1][0] / nu_oil)*(-sin(nodes[0]));
	}

	void StraightRun(int n)
	{
		//cout << endl << "       Start straight_run" << endl;
		y.resize(n);

		double q = 0;
		for (int i = 0; i < n; i++)
		{
			int m = ia[i + 1];
			int w = ia[i];
			int j0 = i - (m - w);
			y[i] = globalB[i];
			int j = 0;
			for (int k = w; k < m; k++, j++)
			{
				q = al[k] * y[j0 + j];;
				y[i] -= q;
			}
			y[i] /= di[i];
		}

		/*cout << endl << "y_vec: ";
		for (int i = 0; i < n; i++)
			cout << y_vec[i] << ' ';
		cout << endl;

		cout << endl << "       End straight_run" << endl;*/
	}

	void Reverse(int n)
	{
		//cout << endl << endl << "       Start reverse" << endl;

		int p = 0;
		double w = 0;
		for (int i = n - 2; i >= 0; i--)
		{
			for (int k = ia[i + 1]; k < ia[i + 2]; k++)
			{
				p = k + i + 1 - ia[i + 2];
				w = au[k] * y[i + 1];
				y[p] -= au[k] * y[i + 1];;
			}
		}

		//fpressure << endl << "res_vec: " << endl;
		/*for (int i = 0; i < n; i++)
		{
			fpressure << fixed;
			fpressure.precision(5);
			fpressure << scientific;
			fpressure << y[i] << endl;
		}*/

		/*double delta = 0;
		int www = 0;
		for (int i = 0; i < n; i++)
		{
			if (i % 1 == 0)
			{
				delta += abs(test(nodes[i]) - q[i]);
				www++;
			}
		}
		delta /= www;
		cout << endl << endl << delta;*/
		//cout << endl << "       End reverse" << endl;
	}

	FEM(vector <FiniteElement> elements)
	{
		Portrait(elements);
		matrixG.resize(2);
		al.resize(ia[ia.size() - 1], 0);
		au.resize(ia[ia.size() - 1], 0);
		di.resize(ia.size() - 1);
		globalB.resize(ia.size() - 1);

		double h;
		for (int itemNumber = 0; itemNumber < elements.size(); itemNumber++) // по каждому элементу
		{
			h = elements[itemNumber].xEnd - elements[itemNumber].xBegin;
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
		au = al;

		BoundaryConditions2();
		BoundaryConditions1(elements);
		LUDecomposition(elements.size() + 1);
		StraightRun(elements.size() + 1);
		Reverse(elements.size() + 1);
	}

	FEM()
	{

	}
};