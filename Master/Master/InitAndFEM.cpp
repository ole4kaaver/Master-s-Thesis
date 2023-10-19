#include "InitAndFEM.h"

double InitAndFEM::test(double x)
{
	return x;
}

double InitAndFEM::Function_nu_in_point(double x_polymer)
{
	if (x_polymer < discrete_nu_func[discrete_nu_func.size() / 2][0])
	{
		if (x_polymer < discrete_nu_func[discrete_nu_func.size() / 4][0])
		{
			for (int i = discrete_nu_func.size() / 4 - 1; i >= 0; i--)
			{
				if (x_polymer >= discrete_nu_func[i][0])
				{
					return discrete_nu_func[i][1] + (discrete_nu_func[i + 1][1] - discrete_nu_func[i][1])
						* (x_polymer - discrete_nu_func[i][0]) / (discrete_nu_func[i + 1][0] - discrete_nu_func[i][0]);
				}
			}
		}
		else
		{
			for (int i = discrete_nu_func.size() / 4 + 1; i <= discrete_nu_func.size() / 2; i++)
			{
				if (x_polymer <= discrete_nu_func[i][0])
				{
					return discrete_nu_func[i - 1][1] + (discrete_nu_func[i][1] - discrete_nu_func[i - 1][1])
						* (x_polymer - discrete_nu_func[i - 1][0]) / (discrete_nu_func[i][0] - discrete_nu_func[i - 1][0]);
				}
			}
		}
	}
	else
	{
		if (x_polymer < discrete_nu_func[discrete_nu_func.size() - discrete_nu_func.size() / 4][0])
		{
			for (int i = discrete_nu_func.size() - discrete_nu_func.size() / 4 - 1; i >= discrete_nu_func.size() / 2; i--)
			{
				if (x_polymer >= discrete_nu_func[i][0])
				{
					return discrete_nu_func[i][1] + (discrete_nu_func[i + 1][1] - discrete_nu_func[i][1])
						* (x_polymer - discrete_nu_func[i][0]) / (discrete_nu_func[i + 1][0] - discrete_nu_func[i][0]);
				}
			}
		}
		else
		{
			for (int i = discrete_nu_func.size() - discrete_nu_func.size() / 4 + 1; i < discrete_nu_func.size(); i++)
			{
				if (x_polymer <= discrete_nu_func[i][0])
				{
					return discrete_nu_func[i - 1][1] + (discrete_nu_func[i][1] - discrete_nu_func[i - 1][1])
						* (x_polymer - discrete_nu_func[i - 1][0]) / (discrete_nu_func[i][0] - discrete_nu_func[i - 1][0]);
				}
			}
		}
	}
}

void InitAndFEM::Reading_viscosity_of_water(string nu_water)
{
	ifstream f_nu_water;

	f_nu_water.open(nu_water + ".txt");

	int discrete_nu_func_size;
	f_nu_water >> discrete_nu_func_size;
	discrete_nu_func.resize(discrete_nu_func_size);
	for (int i = 0; i < discrete_nu_func_size; i++)
	{
		discrete_nu_func[i].resize(2);
		f_nu_water >> discrete_nu_func[i][0] >> discrete_nu_func[i][1];
	}
}

void InitAndFEM::Reading_grid(string grid)
{
	ifstream f_grid;

	f_grid.open(grid + ".txt");

	f_grid >> a >> b >> h >> k_discharge;
	f_grid >> a_time >> b_time >> dt;
}

void InitAndFEM::Reading_parameters(string parameters)
{
	ifstream f_param;

	f_param.open(parameters + ".txt");

	f_param >> K;
	f_param >> F;
	f_param >> nu_oil;
	f_param >> M[0] >> M[1] >> M[2] >> M[3];
	f_param >> x[1][0] >> x[2][0];
	f_param >> s_s;
	f_param >> Q_0;
	f_param >> density[0] >> density[1];
	f_param >> s_crit_oil;
	f_param >> s_res_oil;
}

void InitAndFEM::Building_grid()
{
	cout << "Grid:";
	for (double w = a; w < b - eps1; w += h)
	{
		nodes.push_back(w);
		h *= k_discharge;
		cout << endl << w << " ";
	}

	nodes.push_back(b);
	n = nodes.size();
	cout << endl << b << endl;
}

void InitAndFEM::Resize()
{
	s.resize(2), k.resize(2), alfa.resize(2), amount_of_substance.resize(4), x.resize(3);
	for (size_t m = 0; m < 2; m++)
	{
		s[m].resize(n);
		k[m].resize(n);
		alfa[m].resize(n);
	}
	for (size_t j = 0; j < x.size(); j++)
		x[j].resize(n);

	for (size_t j = 0; j < amount_of_substance.size(); j++)
		amount_of_substance[j].resize(n);

	k_max.resize(2);
	G.resize(2);
	M.resize(4);
	summ_v.resize(2);
	nu_water.resize(n);
	density.resize(2);
	q.resize(n);
	di.resize(n);
	ia.resize(n + 1);
	pp.resize(n);
	global_b.resize(n, 0.0);
	Q.resize(n);
	elems.resize(n - 1);
	s_res_oil_vec.resize(n);
	s_transitional.resize(n);

	int j = 0;
	for (int l = 0; l < n - 1; l++)
	{
		elems[l].resize(2);
		elems[l][0] = j;
		elems[l][1] = j + 1;
		j++;
	}

	Q_m.resize(2), V_m.resize(2), V_out.resize(2);
	for (int l = 0; l < 2; l++)
	{
		Q_m[l].resize(n);
		V_m[l].resize(n);
		V_out[l].resize(n);
	}
}

void InitAndFEM::Init()
{
	for (size_t e = 0; e < n; e++)
	{
		s[0][e] = 0.3;
		s[1][e] = 0.7;
		k[0][e] = s[0][e];
		k[1][e] = s[1][e];
		s_res_oil_vec[e] = s_res_oil;
	}
	x[0][0] = 1.0 - x[1][0] - x[2][0];
	nu_water[0] = Function_nu_in_point(x[1][0]);
	for (size_t e = 1; e < n; e++)
	{
		x[2][e] = 0.0;
		x[1][e] = 0.0;
		x[0][e] = 1.0;
		nu_water[e] = Function_nu_in_point(x[1][e]);
	}
	k_max[0] = 1.0;
	k_max[1] = 1.0;
}

void InitAndFEM::Portrait()
{
	vector<set<int>> portrait(n);
	int ja_size = 0;
	for (int k = 0; k < nodes.size(); k++)
	{
		for (int i = 0; i < elems.size(); i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (elems[i][0] == k || elems[i][1] == k)
					if (elems[i][j] < k)// пробегаемся по элементу, если есть узел меньше, то добавляем его
						portrait[k].insert(elems[i][j]);
			}
		}
	}

	//Подсчёт кол-ва элесментов для ja 
	//cout << "\t-----Portrait-----" << endl;
	for (int i = 0; i < portrait.size(); i++)
	{
		//cout << i << ": ";
		for (auto j = portrait[i].begin(); j != portrait[i].end(); j++)
		{
			//cout << *j << " ";
			ja_size++;
		}
		//cout << endl;
	}

	ja.resize(ja_size);

	//. Заполнение векторов ia и ja
	for (int i = 2; i < ia.size(); i++)
		ia[i] = ia[i - 1] + portrait[i - 1].size();

	for (int index = 0, i = 1; i < portrait.size(); i++)
		for (int temp : portrait[i])// темп берёт в себя каждое значение портрета
			ja[index++] = temp;
}

void InitAndFEM::Building_G_local()
{
	for (int ii = 0; ii < 2; ii++)
		G[ii].resize(2);
	G[0][0] = K * (k[0][ke] / nu_water[ke] + k[1][ke] / nu_oil) / (1.0 * h);
	G[0][1] = -K * (k[0][ke] / nu_water[ke] + k[1][ke] / nu_oil) / (1.0 * h);
	G[1][0] = -K * (k[0][ke] / nu_water[ke] + k[1][ke] / nu_oil) / (1.0 * h);
	G[1][1] = K * (k[0][ke] / nu_water[ke] + k[1][ke] / nu_oil) / (1.0 * h);
}

void InitAndFEM::Boundary_conditions_1(int t_cur) // стр 237
{
	double B = 1e+10, sum = 0;
	double x1, x2;
	x1 = nodes[0];
	x2 = nodes[nodes.size() - 1];
	double res1 = test(x1) * 101325.0;
	double res2 = test(x2) * 101325.0;
	/*double res1 = test(x1);
	double res2 = test(x2);*/
	di[nodes.size() - 1] = B;
	global_b[nodes.size() - 1] = B * res2;
	/*di[0] = B;
	global_b[0] = B * res1;*/
}

void InitAndFEM::Boundary_conditions_2()
{
	global_b[0] += Q_0;
	//global_b[0] -= K * (k[0][0] / nu_water[0] + k[1][0] / nu_oil);
	//global_b[0] -= K * (k[0][0] / nu_water[0] + k[1][0] / nu_oil)*(-sin(nodes[0]));
}

void InitAndFEM::LU_decomposition()
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


	/*cout << endl << "di: ";
	for (int i = 0; i < n; i++)
		cout << di[i] << ' ';

	cout << endl << "al: ";
	for (int i = 0; i < ia[n]; i++)
		cout << al[i] << ' ';

	cout << endl << "au: ";
	for (int i = 0; i < ia[n]; i++)
		cout << au[i] << ' ';
	cout << endl;

	cout << endl << "       End LU_decomposition" << endl << endl;*/
}

void InitAndFEM::straight_run()
{
	//cout << endl << "       Start straight_run" << endl;
	y_vec.resize(n);

	double q = 0;
	for (int i = 0; i < n; i++)
	{
		int m = ia[i + 1];
		int w = ia[i];
		int j0 = i - (m - w);
		y_vec[i] = global_b[i];
		int j = 0;
		for (int k = w; k < m; k++, j++)
		{
			q = al[k] * y_vec[j0 + j];;
			y_vec[i] -= q;
		}
		y_vec[i] /= di[i];
	}

	/*cout << endl << "y_vec: ";
	for (int i = 0; i < n; i++)
		cout << y_vec[i] << ' ';
	cout << endl;

	cout << endl << "       End straight_run" << endl;*/
}

void InitAndFEM::reverse()
{
	//cout << endl << endl << "       Start reverse" << endl;

	int p = 0;
	double w = 0;
	for (int i = n - 2; i >= 0; i--)
	{
		for (int k = ia[i + 1]; k < ia[i + 2]; k++)
		{
			p = k + i + 1 - ia[i + 2];
			w = au[k] * y_vec[i + 1];
			y_vec[p] -= au[k] * y_vec[i + 1];;
		}
	}

	fpressure << endl << "res_vec: " << endl;
	for (int i = 0; i < n; i++)
	{
		fpressure << fixed;
		fpressure.precision(5);
		fpressure << scientific;
		fpressure << y_vec[i] << endl;
	}

	q = y_vec;
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

void InitAndFEM::FEM(int t)
{
	Portrait();
	int k = 0;
	al.resize(ia[n], 0); au.resize(ia[n], 0);

	for (int k = 0; k < nodes.size() - 1; k++) // по каждому элементу
	{
		Building_G_local();

		for (int i = 0; i < 2; i++)
		{
			int i_i = elems[k][i];
			di[k + i] += G[i][i];
			for (int j = 0; j < i; j++)
			{
				int j_j = elems[k][j];
				if (i_i < j_j)
					swap(i_i, j_j);
				for (k = ia[i_i]; ja[k] - j_j != 0; k++);
				al[k] += G[i][j]; // для глобальной матрицы A
			}
		}
		ke++;
	}

	/*cout << endl << endl << "al: ";
	for (int i = 0; i < al.size(); i++)
		cout << al[i] << " ";

	cout << endl<< "di: ";
	for (int i = 0; i < di.size(); i++)
		cout << di[i] << " ";

	cout << endl << "global_b: ";
	for (int i = 0; i < global_b.size(); i++)
		cout << global_b[i] << " ";*/

	Boundary_conditions_2();
	Boundary_conditions_1(t);

	au = al;

	/*cout << endl << endl<<"al: ";
	for (int i = 0; i < al.size(); i++)
		cout << al[i] << " ";

	cout << endl <<"di: ";
	for (int i = 0; i < di.size(); i++)
		cout << di[i] << " ";

	cout << endl <<"global_b: ";
	for (int i = 0; i < global_b.size(); i++)
		cout << global_b[i] << " ";*/

	LU_decomposition();
	straight_run();
	reverse();

	for (int i = 0; i < n; i++)
		global_b[i] = di[i] = 0;

	for (int i = 0; i < ia[n]; i++)
		al[i] = 0;
}

