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

	Component(int num, double massFrac, double amountOfSub)
	{
		number = num;
		massFraction = massFrac;
		amountOfSubstance = amountOfSub;
	}
};

class Phase
{
public:
	int number;
	double saturation;
	double multiplierToPhasePermeability;
	double viscosity;
	double density;

	Phase(int num, double satur, double multiplierToPhasePerm, double viscos, double den)
	{
		number = num;
		saturation = satur;
		multiplierToPhasePermeability = multiplierToPhasePerm;
		viscosity = viscos;
		density = den;
	}
};

class FiniteElement
{
public:
	int number;
	double xBegin, xEnd;
	list <Phase> phases;
	list <Component> components;

	FiniteElement(int num, double x1, double x2, list <Phase> phasesCur, list <Component> componentsCur)
	{
		number = num;
		xBegin = x1;
		xEnd = x2;
		phases = phasesCur;
		components = componentsCur;
	}
};

class Init
{
public:
	double gridBegin, gridEnd, dischargeRatio;
	int numberOfPartitions;
	double tBegin, tEnd, dt;
	vector <FiniteElement> elements;

	void ReadingGrid(string grid)
	{
		ifstream fGrid;
		fGrid.open(grid + ".txt");
		fGrid >> gridBegin >> gridEnd >> numberOfPartitions >> dischargeRatio;
		fGrid >> tBegin >> tEnd >> dt;
	}

	void BuildingGrid(vector <FiniteElement> elements)
	{
		double h;
		double length = gridEnd - gridBegin;
		if (dischargeRatio != 1)
		{
			if (dischargeRatio < 0)
				dischargeRatio = 1 / abs(dischargeRatio);
			h = length * (1 - dischargeRatio) / (1 - pow(dischargeRatio, numberOfPartitions));
		}
		else h = length / numberOfPartitions;
	}

	void ReadingParameters(string parameters)
	{
		ifstream fParam;
		fParam.open(parameters + ".txt");
		int count = 0;
		fParam >> count;
		for (int i = 0; i < count; i++)
		{
			
		}


	}

	Init(string grid, string parameters)
	{
		ReadingGrid(grid);
		elements.resize(numberOfPartitions);
		BuildingGrid(elements);
	}
};