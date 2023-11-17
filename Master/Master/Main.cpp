#pragma once
//#include "Filtration.h"
#include "FEM.h"

//int main()
//{
//	Filtration fil = Filtration("t_grid", "parameters",
//		"grid", "nu_water", "pressure", "saturation", "total_oil", "total_water", "x", "nu_out", "tick", "parity", "s_res", "press");
//	return 0;
//}

int main()
{
	Init Object = Init("grid1", "parameters1", "nuWaterPhase1");
	FEM Method = FEM(Object);
	//CalcDeltaV Calc = CalcDeltaV(Method.y);
}