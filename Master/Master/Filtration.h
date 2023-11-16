#pragma once
#include "InitAndFEM.h"
class Filtration : InitAndFEM
{
public:
	void calc_of_phase_flows();

	void calc_time_step();

	void calc_of_phase_volumes_and_saturations();

	void calc_of_new_x();

	void Calc();

	Filtration(string t_grid, string parameters, string grid, string nu_water,
		string pressure, string saturation, string total_oil, string total_water, string x_, string nu_out,
		string tick, string parity, string s_res, string press) // рассчЄт нового состо€ни€ €чеек
	{
		setlocale(LC_ALL, "Russian");
		Reading_grid(grid);
		Building_grid();
		Resize();
		Reading_parameters(parameters);
		Reading_viscosity_of_water(nu_water);
		Init();
		//Portrait();

		fpressure.open(pressure + ".txt");
		fcalc.open(saturation + ".txt");
		ftotal_oil.open(total_oil + ".txt");
		ftotal_water.open(total_water + ".txt");
		fxi.open(x_ + ".txt");
		fnu_out.open(nu_out + ".txt");
		ft_grid.open(t_grid + ".txt");
		fs_res.open(s_res + ".txt");
		ftick.open(tick + ".txt");
		fparity.open(parity + ".txt");
		fpress.open(press + ".txt");
		double time_tick = 0;
		double tick_ = 200;
		int parity_index = 0;

		for (t = a_time; t <= b_time; t += dt)
		{
			int t_help = t;
			/*if (!(t_help/3600==0||(t_help/3600==4)))
			{
				x[1][0] = 0;
			}*/
			/*if (t >= 7200.0)
			{
				x[1][0] = 0;
			}*/
			if (time_tick == tick_)
			{
				time_tick = 0;
			}
			fcalc << endl << "¬рем€: " << t;
			ft_grid << t << endl;
			fpressure << endl << "¬рем€: " << t;
			fxi << endl << "¬рем€: " << t;
			ke = 0;
			double ss;
			FEM(t);
			pp = q;
			Calc();
			if (parity_index % 1 == 0)
			{
				fparity << total_V_out_oil << endl;
			}
			if (time_tick == 0)
			{
				for (size_t e = 1; e < s[1].size(); e++)
				{
					ftick << s[1][e] << " ";
					//ftick << s_res_oil_vec[e] << " ";
				}
				ftick << endl;
			}
			for (size_t e = 0; e < k[0].size(); e++)
			{
				s_transitional[e] = (s[0][e] - 0.0) / (1.0 - s_res_oil_vec[e] - 0.0);

				k[0][e] = 0.05 * pow(s_transitional[e], 1.5);
				k[1][e] = 0.425 * pow(1.0 - s_transitional[e], 2.0);
				if (isnan(k[0][e]) == true)
				{
					cout << "!!";
				}
				if (k[0][e] >= 1.0 or k[1][e] > 1.0)
				{
					cout << "!!";
				}
				if (e % 1 == 0)
				{
					fpress << y_vec[e] << endl;
				}
			}
			time_tick += dt;
			parity_index++;
		}
		fpressure.close();
		fcalc.close();
		ftotal_oil.close();
		ftotal_water.close();
		fxi.close();
		ft_grid.close();
		fnu_out.close();
		fs_res.close();
		ftick.close();
		fparity.close();
		fpress.close();
	}
};

