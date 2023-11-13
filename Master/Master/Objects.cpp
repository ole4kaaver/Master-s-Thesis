#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

class Component
{
public:
	int number;
	double massFraction;
	double amountOfSubstance;
};

class Phase : Component
{
public:
	int number;
	double saturation;
	double multiplierToPhasePermeability;
	double viscosity;
	double density;
};