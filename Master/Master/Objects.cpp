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
	Component()
	{

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
	vector <Component> components;

	Phase(int num, double satur, double multiplierToPhasePerm, double viscos, double den, vector <Component> comp)
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
	int number;
	double xBegin, xEnd;

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
	vector<Phase> phases;
	vector <FiniteElement> elements;

	void ReadingGrid(string grid)
	{
		ifstream fGrid;
		fGrid.open(grid + ".txt");
		fGrid >> gridBegin >> gridEnd >> numberOfPartitions >> dischargeRatio;
		fGrid >> tBegin >> tEnd >> dt;
	}

	void BuildingGrid()
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

		// заполнение массива elements
		double coordCur = 0;
		for (int i = 0; i < numberOfPartitions; i++)
		{
			elements[i] = FiniteElement(i, coordCur, coordCur + h);
			coordCur = coordCur + h;
			h = h * dischargeRatio;
		}
	}

	void ReadingParameters(string parameters)
	{
		ifstream fParam;
		fParam.open(parameters + ".txt");
		int countPhase, countComponent = 0;
		int numberPhase;
		double saturationCur;
		double multiplierToPhasePermeabilityCur;
		double viscosityCur;
		double densityCur;

		int numberComponent;
		double massFractionCur;
		double amountOfSubstanceCur;
		vector<Component> components;
		
		 /*чтение кол-ва фаз
		 цикл по фазам (номер фазы, хар-ки, кол-во составляющих компонент, хар-ки компонент)*/
		fParam >> countPhase;
		phases.resize(countPhase);
		
		for (int i = 0; i < countPhase; i++)
		{
			fParam >> numberPhase >> saturationCur >> multiplierToPhasePermeabilityCur >> viscosityCur >> densityCur;
			fParam >> countComponent;
			components.resize(countComponent);
			for (int j = 0; j < countComponent; j++)
			{
				fParam >> numberComponent >> massFractionCur >> amountOfSubstanceCur;
				components[j] = Component(numberComponent, massFractionCur, amountOfSubstanceCur);
			}
			phases[i] = Phase(numberPhase, saturationCur, multiplierToPhasePermeabilityCur, viscosityCur, densityCur, components);
		}
	}

	Init(string grid, string parameters)
	{
		ReadingGrid(grid);
		elements.resize(numberOfPartitions);
		BuildingGrid();
	}
};