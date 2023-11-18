#include <vector>
#include <iostream>
#include "FEM.h"
using namespace std;

class CalcOverflows
{
public:
	void PhaseFlows(Init Object, vector <double> pressure);

	void TimeStep();

	void PhaseVolumesAndSaturations();

	void �omponentProperties();

	CalcOverflows(Init Object, vector <double> pressure)
	{
		PhaseFlows(Object, pressure);
		TimeStep();
		PhaseVolumesAndSaturations();
		�omponentProperties();

		fxi << endl;
		for (int e = 1; e < n; e++)
		{
			fxi << fixed;
			fxi.precision(3);
			fxi << scientific << x[1][e] << " ";
		}

		total_V_out_oil += V_out[1][n - 1];
		total_V_out_water += V_out[0][n - 1];

		ftotal_oil << total_V_out_oil << endl;
		ftotal_water << total_V_out_water << endl;
		fcalc << endl << endl << "total_V_out_oil: " << total_V_out_oil;
		fcalc << endl << endl << "total_V_out_water: " << total_V_out_water;
		fcalc << endl << "��������� ������������ ����� ����: " << summ_v[0];
		fcalc << endl << "��������� ������������ ����� �����: " << summ_v[1];

		fcalc << endl << endl << "������������ ������ ��� ������ ���� �� ������� ������� ��:" << endl;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < V_out[i].size(); j++)
			{
				fcalc << fixed;
				fcalc.precision(5);
				fcalc << scientific << V_out[i][j] << " ";
			}
			fcalc << endl;
		}
	}
};

