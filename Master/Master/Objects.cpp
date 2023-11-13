#include <vector>
#include <fstream>
#include <iostream>
#include <list>
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
	list <Component> components;
};

class FiniteElement : Phase
{
public:
	int number;
	double xBegin, xEnd;
	list <Phase> phases;
};

class Init : FiniteElement
{
public:
	double gridBegin, gridEnd, dischargeRatio;
	int numberOfPartitions;
	double tBegin, tEnd, dt;

	void ReadingGrid(string grid)
	{
		ifstream fGrid;
		fGrid.open(grid + ".txt");
		fGrid >> gridBegin >> gridEnd >> numberOfPartitions >> dischargeRatio;
		fGrid >> tBegin >> tEnd >> dt;
	}

	Init(string grid)
	{
		ReadingGrid(grid);
	}
};