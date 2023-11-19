#include "Objects.h"

void Init::ReadingGrid(string grid)
{
	ifstream fGrid;
	fGrid.open(grid + ".txt");
	fGrid >> gridBegin >> gridEnd >> numberOfPartitions >> dischargeRatio;
	fGrid >> tBegin >> tEnd >> dt;
}

void Init::BuildingGrid()
{
	elements.resize(numberOfPartitions);
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

void Init::ReadingParameters(string parameters)
{
	ifstream fParam;
	fParam.open(parameters + ".txt");
	int countPhase, countComponent = 0;
	int numberPhase;
	vector <double> saturationCur;
	vector <double> multiplierToPhasePermeabilityCur;
	vector <double> viscosityCur;
	vector <double> densityCur;

	int numberComponent;
	vector <double> massFractionCur;
	vector <double> amountOfSubstanceCur;
	vector<Component> components;

	saturationCur.resize(elements.size() + 1);
	densityCur.resize(elements.size() + 1);
	massFractionCur.resize(elements.size() + 1);
	amountOfSubstanceCur.resize(elements.size() + 1);
	multiplierToPhasePermeabilityCur.resize(elements.size() + 1);
	viscosityCur.resize(elements.size() + 1);

	// чтение параметров породы
	fParam >> porosity >> permeability;
	fParam >> crossSectionalArea;
	fParam >> flow0;
	fParam >> saturationOilCrit >> saturationOilRes;

	/*чтение кол-ва фаз
	цикл по фазам (номер фазы, хар-ки, кол-во составляющих компонент, хар-ки компонент)*/
	fParam >> countPhase;
	phases.resize(countPhase);

	for (int i = 0; i < countPhase; i++)
	{
		fParam >> numberPhase;
		for (int k = 0; k < elements.size() + 1; k++)
			fParam >> saturationCur[k];
		for (int k = 0; k < elements.size() + 1; k++)
			fParam >> densityCur[k];
		for (int k = 0; k < elements.size() + 1; k++)
			fParam >> multiplierToPhasePermeabilityCur[k];
		for (int k = 0; k < elements.size() + 1; k++)
			fParam >> viscosityCur[k];

		fParam >> countComponent;
		components.resize(countComponent);
		for (int j = 0; j < countComponent; j++)
		{
			fParam >> numberComponent;
			for (int k = 0; k < elements.size() + 1; k++)
				fParam >> massFractionCur[k];
			for (int k = 0; k < elements.size() + 1; k++)
				fParam >> amountOfSubstanceCur[k];
			components[j] = Component(numberComponent, massFractionCur, amountOfSubstanceCur);
		}
		phases[i] = Phase(numberPhase, saturationCur, multiplierToPhasePermeabilityCur, viscosityCur, densityCur, components);
	}
}

void Init::Reading_viscosity_of_water(string nuWaterPhase)
{
	ifstream fNuWaterPhase;

	fNuWaterPhase.open(nuWaterPhase + ".txt");

	int discreteNuFuncSize;
	fNuWaterPhase >> discreteNuFuncSize;
	discreteNuFunc.resize(discreteNuFuncSize);
	for (int i = 0; i < discreteNuFuncSize; i++)
	{
		discreteNuFunc[i].resize(2);
		fNuWaterPhase >> discreteNuFunc[i][0] >> discreteNuFunc[i][1];
	}
}