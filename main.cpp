//Author: Santiago A. Flores Roman
//Description: Performs MC (GCMC and NVT) simulations for simple (monomolecular) species in bulk.
//To do list:
//	Make it more flexible for the implementation of more fluid-fluid potentials.
//  Implement the calculation of the radial distribution function.
//  Extend it to simulations of confined fluids.
//	Reorganize code to avoid utilization of global variables.

#include <string>
#include <cstddef>
#include <iostream>
#include <climits>

#include "MC.h"

int main(int argc, char** argv){
	srand((unsigned)time(NULL)); //seed
	string inFileName;
	int nSets, nEquilSets, nStepsPerSet, printEvery;
	double* moveProbs;
	double rand;
	MC mc;

	inFileName = argv[1];
	mc.ReadInputFile(inFileName);
	mc.PrintParams();
	mc.OutputDirectory(true);
	mc.InitialConfig();
	nSets = mc.GetNSets();
	nEquilSets = mc.GetNEquilSets();
	nStepsPerSet = mc.GetNStepsPerSet();
	printEvery = mc.GetPrintEvery();
	moveProbs = mc.GetMCMoveProbabilities();
	mc.MinimizeEnergy();
	for (int set=1; set<=nSets; set++){
		//Reinitialize MC and Widom statistics every set.
		mc.ResetStats();
		for (int step=1; step<=nStepsPerSet; step++){
			rand = mc.Random();
			if (rand <= moveProbs[0]) mc.MoveParticle(); //Try displacement.
			else if (rand <= moveProbs[1]) mc.ChangeVolume(); //Try volume change.
			else if (rand <= moveProbs[2]) mc.ExchangeParticle(); //Try exchange.
			if (set > nEquilSets){
				mc.ComputeWidom();
				mc.ComputeRDF();
			}
		}
		if (set%printEvery == 0) mc.PrintStats(set);
		if ((set%printEvery == 0) && (set > nEquilSets)){
			mc.ComputeChemicalPotential();
			mc.CreateEXYZ(set);
			mc.CreateLogFile(set);
		}
		if (set <= nEquilSets) mc.AdjustMCMoves();
	}
	mc.PrintRDF();
	return EXIT_SUCCESS;
}

