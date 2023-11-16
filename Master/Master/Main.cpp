#pragma once
#include "Filtration.h"
#include "Objects.cpp"
#include "FEM.cpp"

//int main()
//{
//	Filtration fil = Filtration("t_grid", "parameters",
//		"grid", "nu_water", "pressure", "saturation", "total_oil", "total_water", "x", "nu_out", "tick", "parity", "s_res", "press");
//	return 0;
//}

int main()
{
	Init Object = Init("grid", "parameters", "nuWaterPhase");
	FEM Method = FEM(Object.elements);
	Filtration(Method.y);

	return 0;
}