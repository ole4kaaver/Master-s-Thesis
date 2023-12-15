#include "CalcOverflows.h"

void CalcOverflows::PhaseFlows(Init Object, vector <double> pressure)
{
	Object.phases[0].alfa[0] = 1.0;
	Object.phases[1].alfa[0] = 0.0;
	double mult1, mult2;
	for (int e = 1; e < Object.elements.size() + 1; e++)
	{
		mult1 = Object.phases[0].multiplierToPhasePermeability[e] / Object.phases[0].viscosity[e];
		mult2 = Object.phases[1].multiplierToPhasePermeability[e] / Object.phases[1].viscosity[e];
		Object.phases[0].alfa[e] = mult1 / (mult1 + mult2);
		Object.phases[1].alfa[e] = mult2 / (mult1 + mult2);
	}

	Object.flow[0] = Object.flow0 * Object.crossSectionalArea;
	for (int e = 1; e < Object.elements.size() + 1; e++) // расчёт перетекающего потока смеси через грани кэ
	{
		//Q[e] = Q[0];
		mult1 = Object.phases[0].multiplierToPhasePermeability[e] / Object.phases[0].viscosity[e];
		mult2 = Object.phases[1].multiplierToPhasePermeability[e] / Object.phases[1].viscosity[e];
		Object.flow[e] = Object.crossSectionalArea * Object.permeability * (mult1 + mult2) * (pressure[e] - pressure[e - 1]) 
			/ (Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin);
	}

	for (int e = 0; e < Object.elements.size() + 1; e++) // расчёт перетекающего потока фазы через грани кэ
		for (int m = 0; m < 2; m++)
			Object.phases[m].flow[e] = Object.phases[m].alfa[e] * fabs(Object.flow[e]);
}

void CalcOverflows::TimeStep(Init Object)
{
	double dt1;
	for (int e = 1; e < Object.elements.size() + 1; e++)
	{
		if (Object.phases[1].saturation[e] > Object.saturationOilCrit)
		{
			dt1 = fabs(Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin) * Object.crossSectionalArea * 
				(Object.phases[1].saturation[e] - Object.saturationOilRes) * Object.porosity / Object.phases[1].flow[e];
			if (dt1 < Object.dt)
				Object.dt = dt1;
		}
	}
}

void CalcOverflows::PhaseVolumesAndSaturations(Init Object)
{
	for (int e = 0; e < Object.elements.size() + 1; e++)
		for (int m = 0; m < 2; m++)
			Object.phases[m].volumeOut[e] = Object.phases[m].flow[e] * Object.dt;

	bool flag = false;
	double fraction_of_oil = 0.0;
	for (int e = 1; e < Object.elements.size() + 1; e++)
	{
		Object.phases[0].volumeCur[e] = fabs(Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin) 
			* Object.crossSectionalArea * Object.phases[0].saturation[e] * Object.porosity;
		Object.phases[0].volumeCur[e] += Object.phases[0].volumeOut[e - 1] - Object.phases[0].volumeOut[e];
		Object.phases[1].volumeCur[e] = fabs(Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin)
			* Object.crossSectionalArea * Object.phases[1].saturation[e] * Object.porosity;
		Object.phases[1].volumeCur[e] += Object.phases[1].volumeOut[e - 1] - Object.phases[1].volumeOut[e];

		for (int m = 0; m < 2; m++)
			Object.phases[m].saturation[e] = Object.phases[m].volumeCur[e] 
			/ (fabs(Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin) * Object.porosity * Object.crossSectionalArea);
		

		if (Object.phases[1].volumeOut[e] > fabs(Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin) * 
			Object.porosity * Object.crossSectionalArea * (Object.phases[1].saturation[e] - Object.saturationOilRes))
		{
			double V_out_start = Object.phases[1].volumeOut[e];
			double V_out_mix = Object.phases[1].volumeOut[e] + Object.phases[0].volumeOut[e];

			Object.phases[0].volumeOut[e] += Object.phases[1].volumeOut[e];
			Object.phases[1].volumeOut[e] = fabs(Object.elements[e - 1].xEnd - Object.elements[e - 1].xBegin) * 
				Object.porosity * Object.crossSectionalArea * (Object.phases[1].saturation[e] - Object.saturationOilRes);
			Object.phases[0].volumeOut[e] -= Object.phases[1].volumeOut[e];
			double V_out_mix_res = Object.phases[1].volumeOut[e] + Object.phases[0].volumeOut[e];

			Object.phases[1].saturation[e] = Object.saturationOilRes;
		}
	}
}

void CalcOverflows::СomponentProperties(Init Object)
{
	double s_mobile = 0;
	double mass;
	double portion_oil, portion_surfactant;
	double amount_moll_oil, amount_moll_surfactant;
	//for (int e = 1; e < Object.elements.size() + 1; e++) // количество вещества (вода, полимер, пав, нефть)
	//{
	//	for (int l = 0; l < 3; l++) // вода, полимер, пав
	//		amount_of_substance[l][e] = (density[0] * V_out[0][e - 1] * x[l][e - 1] + density[0] * x[l][e] * (fabs(nodes[e] - nodes[e - 1]) * s_s * F * s[0][e] - V_out[0][e])) / M[l];

	//	mass = V_m[1][e] * density[1];
	//	amount_of_substance[3][e] = mass / M[3];

	//	for (int l = 0; l < 3; l++)
	//		x[l][e] = amount_of_substance[l][e] * M[l] / (amount_of_substance[0][e] * M[0] + amount_of_substance[1][e] * M[1] + amount_of_substance[2][e] * M[2]);

	//	//x[0][e] = 1.0 - x[1][e] - x[2][e];

	//	if (x[2][e] > 0.5e-5)
	//	{
	//		portion_oil = int(amount_of_substance[3][e] / 0.00002);
	//		amount_moll_surfactant = 0.0000005 * portion_oil;
	//		amount_moll_oil = portion_oil * 0.00002;
	//		if (amount_moll_surfactant > amount_of_substance[2][e])
	//		{
	//			portion_surfactant = int(amount_of_substance[2][e] / 0.0000005);
	//			amount_moll_oil = 0.00002 * portion_surfactant;
	//		}

	//		double V_new = amount_moll_oil * M[3] / density[1];
	//		s_res_oil_vec[e] -= V_new / (fabs(nodes[e] - nodes[e - 1]) * F * s_s);
	//		/*if (s_res_oil_vec[e] < 0.4)
	//			int y = 0;*/
	//		if (s_res_oil_vec[e] < 0.0)
	//			s_res_oil_vec[e] = 0;
	//	}
	//	s_mobile += (s[1][e] - s_res_oil_vec[e]) * F * s_s * fabs(nodes[1] - nodes[0]);
	//}
	//fs_res << s_mobile << endl;
	//for (int e = 1; e < n - 1; e++) // вязкость
	//	nu_water[e] = Function_nu_in_point(x[1][e]);
	//fnu_out << nu_water[n - 2] << endl;
}

