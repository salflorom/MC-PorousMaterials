#include <cstdio>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include "MC-Integrated.h"
#include "FluidFluidPotentials.h"

using namespace std;

extern void ReadInputFile(string inFileName){
	sim.file = fopen(inFileName.c_str(), "r");
	fscanf(sim.file,"productionSteps: %i", &sim.nSteps);
	fscanf(sim.file,"equilibriumSteps: %i", &sim.nEquilSteps);
	fscanf(sim.file,"printEvery: %i", &sim.printEvery);
	fscanf(sim.file,"ChemicalSymbol: %s", &fluid.symbol);
	fscanf(sim.file,"molarMass: %lf", &fluid.molarMass);
	fscanf(sim.file,"sigma: %lf", &fluid.sigma);
	fscanf(sim.file,"epsilon: %lf", &fluid.epsilon);
	fscanf(sim.file,"rcut: %lf", &fluid.rcut);
	fscanf(sim.file,"numberOfParticles: %i", &fluid.nParts);
	fscanf(sim.file,"density: %lf", &fluid.dens);
	fscanf(sim.file,"temperature: %lf", &fluid.temp);
	fclose(sim.file);
}
extern void ComputeBoxSize(void){
	double boxWidth, mass, volume;

	mass = fluid.nParts/na*fluid.molarMass; //g
	volume = mass/fluid.dens; //cm^3
	boxWidth = cbrt(volume); //cm
	boxWidth *= 1e8; //AA
	fluid.boxWidth = boxWidth;
}
extern void PrintParams(void){
	cout << "Production steps: " << sim.nSteps << endl;
	cout << "Print every n steps: " << 	sim.printEvery << endl;
	cout << "Temperature (K): " << fluid.temp << endl;
	cout << "Molar mass (g/mol): " << fluid.molarMass << endl;
	cout << "sigma (nm): " << fluid.sigma << endl;
	cout << "epsilon (K): " << fluid.epsilon << endl;
	cout << "rcut (nm): " << fluid.epsilon << endl;
	cout << "Number of particles: " << fluid.nParts << endl;
	cout << "Density (g/cm^3): " << fluid.dens << endl;
	cout << "Box width (AA): " << fluid.boxWidth << endl;
}
extern void AssignPossitions(void){
	for (int i=1; i<=fluid.nParts; i++){
		part[i].x = Random()*fluid.boxWidth; //AA
		part[i].x = Random()*fluid.boxWidth; //AA
		part[i].z = Random()*fluid.boxWidth; //AA
	}
}
//extern void RDF(void){
	//double dw = stats.binWidth;

	//for (int i=1; i<=fluid.nParts; i++){
		//stats.rdfx[int(part[i].x/dw)+1]++;
		//stats.rdfy[int(part[i].y/dw)+1]++;
		//stats.rdfz[int(part[i].z/dw)+1]++;

	//}
//}
//extern void PrintRDF(Simulation sim, Stats stats){
	//double pot = 0;
	//FILE *outFile;
	//char outFileName[50];
	//double dw = stats.binWidth;

	//sprintf(outFileName, "temp/rdf_%i.dat", sim.step);
	//outFile = fopen(outFileName, "w");
	//for (int i=1; i<=NBINS; i++){
		//fprintf(outFile, "%f %f", (i-0.5)*dw, redfx[i]/())
	//}
//}
extern int SelectParticle(void){
	int index;

	do{
		index = int(Random()*fluid.nParts)+1;
	}while (index > fluid.nParts);
	part[0] = part[index]; //Save old position of selected particle.
	return index;
}
extern void MoveParticle(int index){
	part[index].x += (Random()*fluid.nParts-0.5)*sim.dx; //AA
	part[index].y += (Random()*fluid.nParts-0.5)*sim.dy; //AA
	part[index].z += (Random()*fluid.nParts-0.5)*sim.dz; //AA
}
// PBC for a box with the origin at the lower left vertex.
extern void PBC(int index){
	part[index].x -= int(part[index].x/fluid.boxWidth) * fluid.boxWidth; //AA
	part[index].y -= int(part[index].y/fluid.boxWidth) * fluid.boxWidth; //AA
	part[index].z -= int(part[index].z/fluid.boxWidth) * fluid.boxWidth; //AA
}
extern void Metropolis(int index){
	double expArg;
	expArg = -EnergyChange(index)/(kb*fluid.temp);
	if (exp(expArg) > Random()){
		stats.acceptance++;
		part[index].obstructed = 0;
	}else{
		part[index] = part[0];
		part[index].obstructed++;
		stats.rejection++;
	}
}
extern double EnergyChange(int index){
	double deltaEner;
	deltaEner = LJ_Energy(index);
	return deltaEner;
}
extern double Distance(int i, int j){
	float dist, dx, dy, dz;

	dx = part[j].x - part[i].x;
	dx -= round(dx/fluid.boxWidth) * fluid.boxWidth;
	dy = part[j].y - part[i].y;
	dy -= round(dy/fluid.boxWidth) * fluid.boxWidth;
	dz = part[j].z - part[i].z;
	dz -= round(dz/fluid.boxWidth) * fluid.boxWidth;
	dist = sqrt(dx*dx + dy*dy + dz*dz);
	return dist;
}
extern void ResetParticle(int index){
	part[index].x = Random()*fluid.boxWidth; //AA
	part[index].y = Random()*fluid.boxWidth; //AA
	part[index].z = Random()*fluid.boxWidth; //AA
	part[index].obstructed = 0;
}
extern void CreateEXYZ(int step){
	string outFileName = "stats/trajectory.xyz";
	double tmp = 0.0;

	sim.file = fopen(outFileName.c_str(), "a");
	if (step == 0){sim.file = fopen(outFileName.c_str(), "w");}
	fprintf(sim.file, "%i\n%%PBC\n", fluid.nParts);
	for (int i=1; i<=fluid.nParts; i++){
		fprintf(sim.file, "\t%s\t%lf\t%lf\t%lf\n", fluid.symbol, part[i].x, part[i].y, part[i].z); //AA
	}
	fprintf(sim.file, "\n");
	fprintf(sim.file, "Vector1\t%lf\t%lf\t%lf\n", fluid.boxWidth, tmp, tmp); //AA
	fprintf(sim.file, "Vector2\t%lf\t%lf\t%lf\n", tmp, fluid.boxWidth, tmp); //AA
	fprintf(sim.file, "Vector3\t%lf\t%lf\t%lf\n", tmp, tmp, fluid.boxWidth); //AA
	fprintf(sim.file, "Offset\t%lf\t%lf\t%lf\n", tmp, tmp, tmp); //AA
}
extern void PrintStats(int step){
	cout << "Step: " << step;
	cout << " Acceptance rate: " << stats.acceptance*1./step;
	cout << " Rejection rate: " << stats.rejection*1./step;
	cout << " Obstruction rate: " << stats.obstruction*1./step;
	cout << endl;
}
int main(int argc, char** argv){
	srand((unsigned)time(NULL)); //seed
	string inFileName;
	int particle;

	inFileName = argv[1];
	system("mkdir stats");
	ReadInputFile(inFileName);
	ComputeBoxSize();
	stats.binWidth = fluid.boxWidth/NBINS;
	PrintParams();
	system("cd stats");
	//sim.file = fopen("stats/energy.dat", "w");
	//fclose(sim.file);
	AssignPossitions();
	CreateEXYZ(0);
	RDF();
	//PrintRDF();
	for (int step=1; step<=sim.nSteps; step++){
		particle = SelectParticle();
		MoveParticle(particle);
		PBC(particle);
		Metropolis(particle);
		if (part[particle].obstructed >= sim.obstructionLimit){
			ResetParticle(particle);
			stats.obstruction++;
		}
		if (step >= sim.nEquilSteps) RDF();
		if (step%sim.printEvery == 0){
			//PrintRDF();
			PrintStats(step);
			CreateEXYZ(step);
		}
	}
	return 0;
}
// EoS
