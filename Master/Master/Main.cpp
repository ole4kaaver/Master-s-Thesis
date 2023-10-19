#pragma once
#include "Filtration.h"

int main()
{
	Filtration fil = Filtration("t_grid", "parameters",
		"grid", "nu_water", "pressure", "saturation", "total_oil", "total_water", "x", "nu_out", "tick", "parity", "s_res", "press");
	return 0;
}