//Author: Santiago A. Flores Roman
//Description: Performs MC-NVT simulations for simple (monomolecular) species in bulk.
//To do list:
//	Make it more flexible for the implementation of more fluid-fluid potentials.
//	Extend the program to muVT simulations.
//  Implement the calculation of the radial distribution function.
//	Reorganize code to avoid utilization of global variables.

#include <string>
#include <cstddef>
#include <iostream>
#include <climits>

#include "MC.h"

int main(int argc, char** argv){
	srand((unsigned)time(NULL)); //seed
	string inFileName;
	int particle, nSteps, nInitSteps, nEquilSteps, printEvery, obstructed;
	int obstructionLimit=10;
	MC mc;

	system("rm -r stats; mkdir stats");

	inFileName = argv[1];
	mc.ReadInputFile(inFileName);
	mc.ComputeBoxSize();
	mc.PrintParams();
	mc.AssignPossitions();
	nSteps = mc.GetNSteps();
	nInitSteps = mc.GetNInitSteps();
	nEquilSteps = mc.GetNEquilSteps();
	printEvery = mc.GetPrintEvery();
	for (int step=1; step<=nSteps; step++){
		if (step >= nInitSteps) obstructionLimit = INT_MAX;
		particle = mc.SelectParticle();
		mc.MoveParticle(particle);
		mc.Metropolis(particle);
		obstructed = mc.GetObstructed(particle);
		if (obstructed >= obstructionLimit){
			mc.ResetParticle(particle);
			mc.IncrementObstruction();
		}
		if (step >= nEquilSteps){
			mc.ComputeChemicalPotential();
			if (step%printEvery == 0){
				mc.PrintStats(step);
				mc.CreateEXYZ(step);
				mc.CreateLogFile(step);
			}
		}
	}
	return EXIT_SUCCESS;
}

