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
	int particle, nSets, nEquilSets, nStepsPerSet, printEvery, obstructed;
	int obstructionLimit=10;
	double* moveProbs;
	double rand;
	MC mc;

	system("rm -r stats; mkdir stats");

	inFileName = argv[1];
	mc.ReadInputFile(inFileName);
	mc.ComputeBoxSize();
	mc.PrintParams();
	mc.InitialConfig();
	nSets = mc.GetNSets();
	nEquilSets = mc.GetNEquilSets();
	nStepsPerSet = mc.GetNStepsPerSet();
	printEvery = mc.GetPrintEvery();
	moveProbs = mc.GetMCMoveProbabilities();
	for (int set=1; set<=nSets; set++){
		//Reinitialize MC and Widom statistics every set.
		mc.ResetStats();
		mc.ResetWidom();
		if (set >= nEquilSets) obstructionLimit = INT_MAX;
		for (int step=1; step<=nStepsPerSet; step++){
			rand = mc.Random();
			if (rand <= moveProbs[0]){ //Try displacement.
				particle = mc.MoveParticle();
				//obstructed = mc.GetObstructed(particle);
				//if (obstructed >= obstructionLimit){
					//mc.ResetParticle(particle);
					//mc.IncrementObstruction();
				//}
			}else if (rand <= moveProbs[1]){ //Try exchange.
				mc.ExchangeParticle();
			}
			if (set >= nEquilSets){
				if (moveProbs[0] == 1) mc.ComputeWidom(); //If NVT ensemble...
			}
		}
		if (set >= nEquilSets){
			if (moveProbs[0] == 1) mc.ComputeChemicalPotential(); //If NVT ensemble...
		}
		if (set%printEvery == 0){
			mc.PrintStats(set);
			mc.CreateEXYZ(set);
			mc.CreateLogFile(set);
		}
	}
	return EXIT_SUCCESS;
}

