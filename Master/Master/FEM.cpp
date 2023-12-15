#include "FEM.h"

double FEM::Test(double x)
{
	return x;
}

double FEM::F(double x)
{
	return 0;
}

void FEM::Portrait(vector <FiniteElement> elements)
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

void FEM::BuildingGLocal(int number, Init Object)
{
	double multiplier = Object.permeability * 
		(Object.phases[0].multiplierToPhasePermeability[number] / Object.phases[0].viscosity[number] 
			+ Object.phases[1].multiplierToPhasePermeability[number] / Object.phases[1].viscosity[number])
		/ (Object.elements[number].xEnd - Object.elements[number].xBegin);
	matrixG[0][0] = multiplier;
	matrixG[0][1] = -multiplier;
	matrixG[1][0] = -multiplier;
	matrixG[1][1] = multiplier;
}

void FEM::BuildingBLocal(int number, Init Object)
{
	double multiplier = (Object.elements[number].xEnd - Object.elements[number].xBegin) / 6.0;
	localB[0] = multiplier * (2.0 * F(Object.elements[number].xBegin) + F(Object.elements[number].xEnd));
	localB[1] = multiplier * (F(Object.elements[number].xBegin) + 2.0 * F(Object.elements[number].xEnd));
}

void FEM::LUDecomposition(int n)
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

void FEM::BoundaryConditions1(vector <FiniteElement> elem) // стр 237
{
	double B = 1e+10, sum = 0;
	double x1, x2;
	x1 = elem[0].xBegin;
	x2 = elem.back().xEnd;
	double res1 = Test(x1) * 101325.0;
	double res2 = Test(x2) * 101325.0;
	/*double res1 = test(x1);
	double res2 = test(x2);*/
	di[elem.size()] = B;
	globalB[elem.size()] = B * res2;
	/*di[0] = B;
	global_b[0] = B * res1;*/
}

void FEM::BoundaryConditions2(double Q0)
{
	globalB[0] += Q0;
	//global_b[0] -= K * (k[0][0] / nu_water[0] + k[1][0] / nu_oil);
	//global_b[0] -= K * (k[0][0] / nu_water[0] + k[1][0] / nu_oil)*(-sin(nodes[0]));
}

void FEM::StraightRun(int n)
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

void FEM::Reverse(int n)
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

	cout << endl << "res_vec: " << endl;
	for (int i = 0; i < n; i++)
	{
		cout << fixed;
		cout.precision(5);
		cout << scientific;
		cout << y[i] << endl;
	}

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