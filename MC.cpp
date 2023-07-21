//Author: Santiago A. Flores Roman

//#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <regex>
#include <cctype>

#include "operations.h"
#include "MC.h"

using namespace std;

MC::MC(void){
	int i;
	stats.acceptance = stats.rejection = stats.obstruction = 0;
	stats.nDisplacements = stats.widomInsertions = stats.insertionParam = 0;
	fluid.mu = 0;
	sim.dx = 1, sim.dy = 1, sim.dz = 1;
	sim.nInitSteps = sim.nEquilSteps = sim.nSteps = 0;
	for (i=0; i<=NBINS; i++){
		stats.rdfx[i] = 0;
		stats.rdfy[i] = 0;
		stats.rdfz[i] = 0;
	}
	for (i=0; i<3; i++) {sim.PBC[i] = true;}
	for (i=0; i<MAXPART; i++) {part[i].obstructed = 0;}
	fluid.epsilon = fluid.sigma = fluid.rcut = 0;
	fluid.nParts = fluid.dens = fluid.temp = fluid.molarMass = 0;
	fluid.initialE = fluid.finalE = 0;
	fluid.name = fluid.vdwPot = "";
}
void MC::ReadInputFile(string inFileName){
	string line, command, value;
	Operations ops;
	ifstream inFile;
	regex match("^\\s*([a-z0-9_]+)\\s+?([0-9_\\.a-z@,\\s\\-+]+);.*",regex_constants::icase);
	smatch params;

	inFile.open(inFileName);
	if (!inFile.is_open()){ //Check if file was opened successfully.
		cout << "Error opening file" << endl;
		exit(EXIT_FAILURE);
	}
	while (getline(inFile, line)) {
		if (regex_search(line, params, match)){
			command = ops.LowerCase(params[1]);
			value = ops.LowerCase(params[2]);
		}
		//cout << command << " " << value << endl;
		if (command == "productionsteps") {sim.nSteps = stod(value);}
		else if (command == "initializationsteps") {sim.nInitSteps = stod(value);}
		else if (command == "equilibriumsteps") {sim.nEquilSteps = stod(value);}
		else if (command == "printevery") {sim.printEvery = stod(value);}
		else if (command == "temperature") {fluid.temp = stod(value);}
		else if (command == "name") {fluid.name = value;}
		else if (command == "numberofparticles") {fluid.nParts = stod(value);}
		else if (command == "density") {fluid.dens = stod(value);}
		else if (command == "molarmass") {fluid.molarMass = stod(value);}
		else if (command == "vdwpotential") {fluid.vdwPot = value;}
		else if (command == "sigma") {fluid.sigma = stod(value);}
		else if (command == "epsilon") {fluid.epsilon = stod(value);}
		else if (command == "rcut") {fluid.rcut = stod(value);}
	}
	inFile.close();
}
void MC::ComputeBoxSize(void){
	double boxWidth, mass, volume;
	int nPerSide=round(cbrt(fluid.nParts));

	fluid.nParts = nPerSide*nPerSide*nPerSide;
	mass = fluid.nParts/na*fluid.molarMass; //g
	volume = mass/fluid.dens; //cm^3
	boxWidth = cbrt(volume); //cm
	boxWidth *= 1e8; //AA
	fluid.boxWidth = boxWidth;
	//stats.binWidth = boxWidth/NBINS;
	stats.binWidth = fluid.rcut/NBINS;
	sim.dx = sim.dy = sim.dz = fluid.sigma; //AA. Step size set to sigma.
}
void MC::PrintParams(void){
	cout << "Production steps: " << sim.nSteps << endl;
	cout << "Initialization steps: " << sim.nInitSteps << endl;
	cout << "Equilibration steps: " << sim.nEquilSteps << endl;
	cout << "Print every n steps: " << 	sim.printEvery << endl;
	cout << "Temperature (K): " << fluid.temp << endl;
	cout << "Fluid name: " << fluid.name << endl;
	cout << "Molar mass (g/mol): " << fluid.molarMass << endl;
	cout << "VdW potential: " << fluid.vdwPot << endl;
	cout << "sigma (AA): " << fluid.sigma << endl;
	cout << "epsilon (K): " << fluid.epsilon << endl;
	cout << "rcut (AA): " << fluid.rcut << endl;
	cout << "Number of particles: " << fluid.nParts << endl;
	cout << "Density (g/cm^3): " << fluid.dens << endl;
	cout << "Box width (AA): " << fluid.boxWidth << endl;
	cout << endl;
}
void MC::AssignPossitions(void){
	int nPerSide=round(cbrt(fluid.nParts));
	double interMolSpace=fluid.boxWidth/nPerSide;

	for (int k=0; k<nPerSide; k++){
		for (int j=0; j<nPerSide; j++){
			for (int i=0; i<nPerSide; i++){
				part[1+i+j*nPerSide+k*nPerSide*nPerSide].x = (1 + 2*i + ((j+k)%2))*0.5*interMolSpace;
				part[1+i+j*nPerSide+k*nPerSide*nPerSide].y = (1 + sqrt(3) * (j+(k%2)/3.))*0.5*interMolSpace;
				part[1+i+j*nPerSide+k*nPerSide*nPerSide].z = (1 + 2*sqrt(6)/3. * k)*0.5*interMolSpace;
			}
		}
	}
}
int MC::SelectParticle(void){
	int index;

	do{
		index = int(Random()*fluid.nParts)+1;
	}while (index > fluid.nParts);
	part[0] = part[index]; //Save old position of selected particle.
	fluid.initialE = EnergyOfParticle(index);
	return index;
}
void MC::MoveParticle(int index){
	part[index].x += (Random()-0.5)*sim.dx; //AA
	part[index].y += (Random()-0.5)*sim.dy; //AA
	part[index].z += (Random()-0.5)*sim.dz; //AA
	PBC(index);
	fluid.finalE = EnergyOfParticle(index);
	stats.nDisplacements++;
}
// PBC for a box with the origin at the lower left vertex.
void MC::PBC(int index){
	part[index].x -= floor(part[index].x/fluid.boxWidth) * fluid.boxWidth; //AA
	part[index].y -= floor(part[index].y/fluid.boxWidth) * fluid.boxWidth; //AA
	part[index].z -= floor(part[index].z/fluid.boxWidth) * fluid.boxWidth; //AA
}
void MC::Metropolis(int index){
	double argument;

	//cout << fluid.finalE-fluid.initialE << endl;
	argument = -(fluid.finalE-fluid.initialE)/fluid.temp;
	if (exp(argument) > Random()){
		stats.acceptance++;
		part[index].obstructed = 0;
	}else{
		part[index] = part[0];
		part[index].obstructed++;
		stats.rejection++;
	}
}
double MC::EnergyOfParticle(int index){
	double energy=0;

	//Fluid-Fluid energy
	if (fluid.vdwPot == "lj"){ //Lennard-Jones potential.
		for (int i=1; i<=fluid.nParts; i++) energy += LJ_Energy(index,i);
		energy *= 0.5;
	}else if (fluid.vdwPot == "eam_ga"){ //EAM potential for Ga.
		energy += EAMGA_Energy(index);
	}
	return energy;
}
double MC::SystemEnergy(void){
	double energy=0;

	for (int i=1; i<=fluid.nParts; i++) energy += EnergyOfParticle(i);
	return energy;
}
void MC::ComputeChemicalPotential(void){
	double thermalWL, partMass, rho, muIdeal, muExcess, argument;

	stats.widomInsertions++;
	partMass = fluid.molarMass/na*1e3; // kg/particle
	rho = fluid.dens/fluid.molarMass*na*1e6; // particle/m^3
	thermalWL = plank/sqrt(2*pi*partMass*kb*fluid.temp); //m^3
	muIdeal = fluid.temp*log(rho*thermalWL*thermalWL*thermalWL); //K
	// Virtual particle.
	part[MAXPART-1].x = Random()*fluid.boxWidth;
	part[MAXPART-1].y = Random()*fluid.boxWidth;
	part[MAXPART-1].z = Random()*fluid.boxWidth;
	argument = exp(-EnergyOfParticle(MAXPART-1)/fluid.temp);
	stats.insertionParam += argument/stats.widomInsertions;
	muExcess = -fluid.temp*log(stats.insertionParam); //K
	fluid.mu = muIdeal+muExcess; //K
}
void MC::ResetParticle(int index){
	part[index].x = Random()*fluid.boxWidth; //AA
	part[index].y = Random()*fluid.boxWidth; //AA
	part[index].z = Random()*fluid.boxWidth; //AA
}
void MC::RDF(void){
	for (int i=1; i<=int(fluid.nParts); i++){
		stats.rdfx[int(part[i].x/stats.binWidth)+1]++;
		stats.rdfy[int(part[i].y/stats.binWidth)+1]++;
		stats.rdfz[int(part[i].z/stats.binWidth)+1]++;

	}
}
void MC::PrintRDF(int step){
	Operations ops;
	ofstream rdfFile;
	string outFileName="stats/rdf.dat";
	double dw=stats.binWidth;
	double density=fluid.nParts/NBINS;
	double relativeStep=step+sim.printEvery;

	if (ops.FileExists(outFileName)) rdfFile.open(outFileName, ios::app);
	else rdfFile.open(outFileName);
	rdfFile << NBINS << "\nrx(AA)\tg(rx)\try(AA)\tg(ry)\trz(AA)\tg(rz)" << endl;
	for (int i=1; i<=NBINS; i++){
		rdfFile << (i-0.5)*dw << "\t" << stats.rdfx[i]/(relativeStep*density) << "\t";
		rdfFile << (i-0.5)*dw << "\t" << stats.rdfy[i]/(relativeStep*density) << "\t";
		rdfFile << (i-0.5)*dw << "\t" << stats.rdfz[i]/(relativeStep*density) << endl;
	}
	rdfFile.close();
}
void MC::PrintStats(int step){
	cout << "Step: " << step;
	cout << "\t\tNum. particles: " << fluid.nParts;
	cout << "\t\tEnergy/particle: " << SystemEnergy()/fluid.nParts;
	cout << "\t\tAcceptance rate: " << stats.acceptance*1./step;
	cout << "\t\tRejection rate: " << stats.rejection*1./step;
	cout << "\t\tObstruction rate: " << stats.obstruction*1./step;
	cout << endl;
}
void MC::CreateEXYZ(int step){
	Operations ops;
	ofstream trajFile;
	string outFileName="stats/trajectory.exyz";
	double tmp=0.0;

	if (ops.FileExists(outFileName)) trajFile.open(outFileName, ios::app);
	else trajFile.open(outFileName);
	trajFile << fluid.nParts << "\n";
	trajFile << "Lattice=\" " << fluid.boxWidth << " " << tmp << " " << tmp << " ";
	trajFile << tmp << " " << fluid.boxWidth << " " << tmp << " ";
	trajFile << tmp << " " << tmp << " " << fluid.boxWidth << "\" ";
	trajFile << "Properties=species:S:1:pos:R:3 Time=" << 1.*step << "\n";
	for (int i=1; i<=fluid.nParts; i++){
		trajFile << fluid.name << "\t";
		trajFile << part[i].x << "\t" << part[i].y << "\t" << part[i].z << "\n"; //AA
	}
	trajFile.close();
}
void MC::CreateLogFile(int step){
	Operations ops;
	ofstream logFile;
	string outFileName="stats/simulation.log";
	double energyPerPart;

	if (!ops.FileExists(outFileName)){
		logFile.open(outFileName);
		logFile << "Step\tN parts\tBoxWidth [AA]\tE/part [K]\tmu [K]\n";
	}else{
		energyPerPart = SystemEnergy()/fluid.nParts;
		logFile.open(outFileName, ios::app);
		logFile << step << "\t";
		logFile << fluid.nParts << "\t";
		logFile << fluid.boxWidth << "\t";
		logFile << energyPerPart << "\t";
		logFile	<< fluid.mu << "\n";
	}
	logFile.close();
}
double MC::NeighDistance(int i, int j){
	float dist, dx, dy, dz;
	float boxWidth=fluid.boxWidth;

	dx = part[j].x - part[i].x; //AA
	dx -= round(dx/boxWidth) * boxWidth; //AA
	dy = part[j].y - part[i].y; //AA.
	dy -= round(dy/boxWidth) * boxWidth; //AA
	dz = part[j].z - part[i].z; //AA.
	dz -= round(dz/boxWidth) * boxWidth; //AA
	dist = sqrt(dx*dx + dy*dy + dz*dz); //AA
	return dist; //AA
}
double MC::LJ_Energy(int i, int j){
	double rij;
	double eps=fluid.epsilon, sig=fluid.sigma, rcut=fluid.rcut;
	double energy=0;

	if (i != j) {
		rij = NeighDistance(i, j);
		if (rij < rcut) energy += 4*eps*(pow(sig/rij,12)-pow(sig/rij,6));
	}
	return energy; //K
}
double MC::EAMGA_Energy(int index){
	double rij, energy;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	for (int j=1; j<=fluid.nParts; j++){
		if (j != index){
			rij = NeighDistance(index, j);
			if (rij < fluid.rcut){
				rho += eDens(rij);
				phiLow += PairPot(rij);
			}
		}
	}
	energy = EmbPot(rho) + phiLow;
	energy *= eVToK;
	return energy;
}
double MC::StepUnit(double radius, double leftLim, double rightLim){
	if (leftLim <= radius && radius < rightLim) return 1.;
	else return 0.;
}
double MC::EmbPot(double rho){
	double rhoIntervals[7] = {1.00000,  0.92000,  0.870000,  0.800000,  0.750000,  0.650000,  1.400000};
	double aValues[7]      = {0.00000, -1.91235, -1.904030, -1.897380, -1.883520, -1.852620, -1.822820};
	double bValues[7]      = {0.00000,  0.00000, -0.208000, -0.058000, -0.338000, -0.898000,  0.302000};
	double cValues[7]      = {0.00000,  1.30000, -1.500000,  2.000000,  5.600000, -6.000000,  2.000000};
	double phi = 0;

	if (rhoIntervals[1] <= rho && rho <= rhoIntervals[6]) phi = aValues[1] + cValues[1]*(rho-rhoIntervals[0])*(rho-rhoIntervals[0]);
	for (int i=2; i<6; i++){
		if (rhoIntervals[i] <= rho && rho <= rhoIntervals[i-1]){
			phi = aValues[i] + bValues[i]*(rho-rhoIntervals[i-1]) + cValues[i]*(rho-rhoIntervals[i-1])*(rho-rhoIntervals[i-1]);
			break;
		}
	}
	if (rho <= rhoIntervals[5]){
		phi = (aValues[6] + bValues[6]*(rho-rhoIntervals[5]) + cValues[6]*(rho-rhoIntervals[5])*(rho-rhoIntervals[5])) \
* (2*rho/rhoIntervals[5]-(rho/rhoIntervals[5])*(rho/rhoIntervals[5]));
	}
	return phi;
}
double MC::eDens(double radius){
	double pValues[3] = {0, 2.24450, 1.2};

	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::PairPot(double radius){
	double rIntervals[7] = {0.0, 2.15, 2.75, 3.35, 4.00, 6.50, 8.30};
	double aValues[9][6] = {{0.0, -0.65052509307861e-01, -0.15576396882534e+00, -0.13794735074043e+00, 0.13303710147738e-01,  0.00000000000000e+00},
							{0.0, -0.32728102803230e+00, -0.16365580260754e-01,  0.78778542578220e-01, 0.59769893996418e-02,  0.00000000000000e+00},
							{0.0,  0.51590444127493e+01,  0.20955204046244e+00, -0.83622260891495e-01, 0.57411338894840e-01, -0.60454444423660e-02},
							{0.0,  0.90195221829217e+02, -0.97550604734748e+00, -0.44410858010987e+01, 0.19517888219051e+00, -0.13258585494287e+00},
							{0.0,  0.72322004859499e+03, -0.11625479189815e+02, -0.36415106938231e+02, 0.32162310059276e+00, -0.34988482891053e+00},
							{0.0,  0.27788989409594e+04, -0.58549935696765e+02, -0.13414583419234e+03, 0.30195698240893e+00, -0.45183606796559e+00},
							{0.0,  0.56037895713613e+04, -0.15186293377510e+03, -0.25239146992011e+03, 0.14850603977640e+00, -0.31733856650298e+00},
							{0.0,  0.57428084950480e+04, -0.19622924502226e+03, -0.23858760191913e+03, 0.36233874262589e-01, -0.11493645479281e+00},
							{0.0,  0.23685488320885e+04, -0.98789413798382e+02, -0.90270667293646e+02, 0.34984220138018e-02, -0.16768950999376e-01}};
	double phi=0;

	if (rIntervals[1] < radius && radius <= rIntervals[6]){
		for (int i=1; i<6; i++){
			for (int m=0; m<9; m++){
				phi += aValues[m][i] * pow((radius-rIntervals[i+1]),m) * StepUnit(radius, rIntervals[i], rIntervals[i+1]);
			}
		}
	}else if (rIntervals[0] < radius && radius <= rIntervals[1]){
		phi = 0.619588 - 51.86268*(2.15-radius) + 27.8*(exp(1.96*(2.15-radius))-1);
	}
	return phi;
}

