#include <vector>
#include <iostream>
#include "FEM.h"
using namespace std;

class CalcOverflows
{
public:
	Init Mixture;

	void PhaseFlows(Init Object, vector <double> pressure);

	void TimeStep(Init Object);

	void PhaseVolumesAndSaturations(Init Object);

	void ÑomponentProperties(Init Object);

	CalcOverflows(Init Object, vector <double> pressure)
	{
		Mixture = Object;
		Mixture.phases[0].alfa.resize(Object.elements.size() + 1);
		Mixture.phases[1].alfa.resize(Object.elements.size() + 1);
		Mixture.flow.resize(Object.elements.size() + 1);
		Mixture.phases[0].flow.resize(Object.elements.size() + 1);
		Mixture.phases[1].flow.resize(Object.elements.size() + 1);
		Mixture.phases[0].volumeOut.resize(Object.elements.size() + 1);
		Mixture.phases[1].volumeOut.resize(Object.elements.size() + 1);
		Mixture.phases[0].volumeCur.resize(Object.elements.size() + 1);
		Mixture.phases[1].volumeCur.resize(Object.elements.size() + 1);

		PhaseFlows(Mixture, pressure);
		TimeStep(Mixture);
		PhaseVolumesAndSaturations(Mixture);
		ÑomponentProperties(Mixture);
	}
};

