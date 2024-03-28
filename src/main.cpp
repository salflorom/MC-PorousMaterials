//Author: Santiago A. Flores Roman
//Description: Performs MC (GCMC, NVT, NPT, GEMC-NPT, GEMC-NVT) simulations for simple (monoatomic or single-bead) species.
//To do list:
//	Reorganize code to avoid utilization of global variables.
//	Baptize the program.

#include <string> // string
#include <cstdlib> // EXIT_SUCCESS
#include <iostream> // cout
#include <chrono> // chrono
#include <ctime> // time_t, ctime
#include "Eigen/Dense"

#include "MC.h"

int main(int argc, char** argv){
	string inFileName;
	size_t nSets, nEquilSets, nStepsPerSet, printEvery;
	int nDispAttempts, nVolAttempts, nSwapAttempts;
	size_t currentSet;
	double rand;
	Eigen::Vector3i moves;
	chrono::time_point<chrono::system_clock> startMinimization, endMinimization;
	chrono::time_point<chrono::system_clock> start, startSimulation, endSimulation;
	chrono::duration<double> elapsedMinimization, elapsedSimulation;
	MC mc;

	start = chrono::system_clock::now();
	time_t startTime = chrono::system_clock::to_time_t(start);
	cout << "Starting simulation on " << ctime(&startTime) << endl;

	inFileName = argv[1];
	currentSet = mc.ReadInputFile(inFileName);
	mc.PrintParams();
	mc.OutputDirectory();
	mc.InitialConfig();
	nSets = mc.GetNSets();
	nEquilSets = mc.GetNEquilSets();
	nStepsPerSet = mc.GetNStepsPerSet();
	printEvery = mc.GetPrintEvery();
	moves = mc.GetMCMoves();
	nDispAttempts = moves(0);
	nVolAttempts = moves(1);
	nSwapAttempts = moves(2);

	if (currentSet == 0) mc.PrintTrajectory(0);

	cout << "Starting minimization..." << endl;
	startMinimization = chrono::system_clock::now();
	mc.MinimizeEnergy(currentSet);
	endMinimization = chrono::system_clock::now();
	elapsedMinimization = endMinimization-startMinimization;
	cout << "Minimization finished" << endl;
	cout << "Elapsed time of minimization: " << elapsedMinimization.count() << " s" << endl;
	cout << endl;

	if (currentSet == 0){
		mc.PrintTrajectory(0);
		mc.PrintLog(0);
	}

	startSimulation = chrono::system_clock::now();
	for (size_t set=currentSet+1; set<=nSets; set++){
		//Reinitialize MC and Widom statistics every set.
		mc.ResetStats();
		for (size_t step=1; step<=nStepsPerSet; step++){
			rand = mc.Random() * (nDispAttempts+nVolAttempts+nSwapAttempts);
			if (rand <= nDispAttempts) mc.MoveParticle(); //Try displacement.
			else if (rand <= nDispAttempts + nVolAttempts) mc.ChangeVolume(); //Try volume change.
			else if (rand <= nDispAttempts + nVolAttempts + nSwapAttempts) mc.SwapParticle(); //Try swap.
			mc.CorrectEnergy();
			if (set > nEquilSets){
				mc.ComputeWidom();
				mc.ComputeRDF();
			}
		}
		if (set%printEvery == 0){
			mc.ComputeChemicalPotential();
			mc.PrintStats(set);
			mc.PrintTrajectory(set);
			mc.PrintLog(set);
		}
		mc.AdjustMCMoves(set);
	}
	endSimulation = chrono::system_clock::now();
	elapsedSimulation = endSimulation-startSimulation;

	mc.PrintRDF();

	time_t endTime = chrono::system_clock::to_time_t(endSimulation);
	cout << endl;
	cout << "Simulation finished" << endl;
	cout << "Elapsed time of simulation: " << elapsedSimulation.count() << " s" << endl;
	cout << "Computation finished at " << ctime(&endTime) << endl;
	cout << endl;
	return EXIT_SUCCESS;
}

