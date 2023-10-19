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

	int n; // ���������� �����
	int ke = 0; // ����� ��
	double a_time, b_time; // ������ � ����� �� ����� �� �������
	double a, b, h, k_discharge; // ������ � ����� �������, ��� ���������, ���� ��������
	double K, F; // ����������� �������������, ����������� ����������
	double nu_oil; // �������� �����
	double Q_0; // ��������� �����
	double s_s; // ������� ������� �����
	double dt, dt1, t; // ��� �� �������
	double total_V_out_oil = 0, total_V_out_water = 0;
	double s_crit_oil, s_res_oil; // ����������� � ���������� ������������
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

