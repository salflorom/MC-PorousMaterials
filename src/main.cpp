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
	int nSets, nEquilSets, nCyclesPerSet, printEvery;
	int nDispAttempts, nVolAttempts, nSwapAttempts;
	int* moves;
	double rand;
	MC mc;

	inFileName = argv[1];
	mc.ReadInputFile(inFileName);
	mc.PrintParams();
	mc.OutputDirectory();
	mc.InitialConfig();
	nSets = mc.GetNSets();
	nEquilSets = mc.GetNEquilSets();
	nCyclesPerSet = mc.GetNCyclesPerSet();
	printEvery = mc.GetPrintEvery();
	moves = mc.GetMCMoves();
	nDispAttempts = moves[0];
	nVolAttempts = moves[1];
	nSwapAttempts = moves[2];
	mc.MinimizeEnergy();
	mc.PrintTrajectory(0);
	mc.PrintLog(0);
	for (int set=1; set<=nSets; set++){
		//Reinitialize MC and Widom statistics every set.
		mc.ResetStats();
		for (int cycle=1; cycle<=nCyclesPerSet; cycle++){
			rand = mc.Random() * (nDispAttempts+nVolAttempts+nSwapAttempts);
			if (rand <= nDispAttempts) mc.MoveParticle(); //Try displacement.
			else if (rand <= nDispAttempts + nVolAttempts) mc.ChangeVolume(); //Try volume change.
			else if (rand <= nDispAttempts + nVolAttempts + nSwapAttempts) mc.SwapParticle(); //Try swap.
			if (set > nEquilSets){
				mc.ComputeWidom();
				mc.ComputeRDF();
			}
		}
		if (set%printEvery == 0) mc.PrintStats(set);
		if ((set%printEvery == 0) && (set > nEquilSets)){
			mc.ComputeChemicalPotential();
			mc.PrintTrajectory(set);
			mc.PrintLog(set);
		}
		if (set <= nEquilSets) mc.AdjustMCMoves();
	}
	mc.PrintRDF();
	return EXIT_SUCCESS;
}

