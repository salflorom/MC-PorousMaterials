#include <string>
#include <cstddef>
#include <iostream>

#include "MC.h"

int main(int argc, char** argv){
	srand((unsigned)time(NULL)); //seed
	string inFileName;
	int particle, nSteps, nEquilSteps, printEvery, obstructed, obstructionLimit;
	//FILE *file;
	MC mc;

	inFileName = argv[1];
	system("mkdir stats");
	mc.ReadInputFile(inFileName);
	mc.ComputeBoxSize();
	mc.PrintParams();
	system("cd stats");
	//file = fopen("stats/energy.dat", "w");
	//fclose(file);
	mc.AssignPossitions();
	mc.RDF();
	//PrintRDF();
	nSteps = mc.GetNSteps();
	nEquilSteps = mc.GetNEquilSteps();
	printEvery = mc.GetPrintEvery();
	obstructionLimit = mc.GetObstructionLimit();
	for (int step=1; step<=nSteps; step++){
		particle = mc.SelectParticle();
		mc.MoveParticle(particle);
		mc.Metropolis(particle);
		obstructed = mc.GetObstructed(particle);
		if (obstructed >= obstructionLimit){
			mc.ResetParticle(particle);
			mc.IncrementObstruction();
		}
		if (step >= nEquilSteps) mc.RDF();
		if ((step%printEvery == 0) && (step >= nEquilSteps)){
			//PrintRDF();
			mc.PrintStats(step);
			mc.CreateEXYZ(step);
		}
	}
	return EXIT_SUCCESS;
}

