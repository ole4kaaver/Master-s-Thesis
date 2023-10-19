#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <set>
using namespace std;

class InitAndFEM
{
public:
	const int max_iter = 10000;
	const double eps = 1e-10;
	const double eps1 = 1e-8;
	const double eps_s = 1e-8;

	int n; // количество узлов
	int ke = 0; // номер кэ
	double a_time, b_time; // начало и конец по сетки по времени
	double a, b, h, k_discharge; // начало и конец области, шаг разбиения, коэф разрядки
	double K, F; // структурная проницаемость, структурная пористость
	double nu_oil; // вязкость нефти
	double Q_0; // начальный поток
	double s_s; // площадь сечения керна
	double dt, dt1, t; // шаг по времени
	double total_V_out_oil = 0, total_V_out_water = 0;
	double s_crit_oil, s_res_oil; // минимальная и остаточные насыщенности
	vector<double> summ_v, M, nu_water, density, k_max, s_transitional, s_res_oil_vec;
	vector<vector<double>> Q_m, V_m, s, V_out, k, alfa, discrete_nu_func, x, amount_of_substance;

	vector<double> nodes;
	vector<double> pp, Q;
	vector<double> global_b;
	vector<vector<double>> G;
	vector<int> ia, ja;
	vector<double> al, au, di, f;
	vector<vector<int>> elems;
	vector<double> y_vec, q;

	ofstream fpressure;
	ofstream fcalc;
	ofstream ftotal_oil;
	ofstream ftotal_water;
	ofstream fxi;
	ofstream ft_grid;
	ofstream fnu_out;
	ofstream fs_res;
	ofstream ftick;
	ofstream fparity;
	ofstream fpress;

	double test(double x);

	double Function_nu_in_point(double x_polymer);

	void Reading_viscosity_of_water(string nu_water);

	void Reading_grid(string grid);

	void Reading_parameters(string parameters);

	void Building_grid();

	void Resize();

	void Init();

	void Portrait();

	void Building_G_local();

	void Boundary_conditions_1(int t_cur);

	void Boundary_conditions_2();

	void LU_decomposition();

	void straight_run();

	void reverse();

	void FEM(int t);
};

