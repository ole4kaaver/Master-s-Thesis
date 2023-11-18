#include <vector>
#include <fstream>
#include <iostream>
#include <list>
using namespace std;

class Component
{
public:
	int number = 0;
	vector <double> massFraction;
	vector <double> amountOfSubstance;

	Component(int num, vector <double> massFrac, vector <double> amountOfSub)
	{
		number = num;
		massFraction = massFrac;
		amountOfSubstance = amountOfSub;
	}
	Component()
	{

	}
};

class Phase
{
public:
	int number = 0;
	vector <double> saturation;
	vector <double> multiplierToPhasePermeability;
	vector <double> viscosity;
	vector <double> density;
	vector <double> alfa;
	vector <Component> components;
	vector <double> flow;
	vector <double> volumeCur, volumeIn, volumeOut;

	Phase(int num, vector <double> satur, vector <double> multiplierToPhasePerm, vector <double> viscos,
		vector <double> den, vector <Component> comp)
	{
		number = num;
		saturation = satur;
		multiplierToPhasePermeability = multiplierToPhasePerm;
		viscosity = viscos;
		density = den;
		components = comp;
	}
	Phase()
	{

	}

};

class FiniteElement
{
public:
	int number = 0;
	double xBegin = 0.0, xEnd = 0.0;

	FiniteElement(int num, double x1, double x2)
	{
		number = num;
		xBegin = x1;
		xEnd = x2;
	}
	FiniteElement()
	{

	}
};

class Init
{
public:
	double gridBegin, gridEnd, dischargeRatio;
	int numberOfPartitions;
	double tBegin, tEnd, dt;
	vector <Phase> phases;
	vector <FiniteElement> elements;
	vector <vector<double>> discreteNuFunc;
	vector <double> flow;
	double flow0; // начальный поток
	double permeability, porosity, crossSectionalArea; // структурная проницаемость, структурная пористость, площадь сечения керна
	double saturationOilCrit, saturationOilRes;

	void ReadingGrid(string grid);

	void BuildingGrid();

	void ReadingParameters(string parameters);

	void Reading_viscosity_of_water(string nuWaterPhase);

	Init(string grid, string parameters, string nuWaterPhase)
	{
		ReadingGrid(grid);
		BuildingGrid();
		ReadingParameters(parameters);
		Reading_viscosity_of_water(nuWaterPhase);
	}

	Init()
	{

	}
};