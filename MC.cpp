#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <regex>
#include <cctype>
#include <vector>

#include "MC.h"
#include "operations.h"
#include "FluidFluidPotentials.h"

using namespace std;

MC::MC(void){
	stats.acceptance = stats.rejection = stats.obstruction = 0;
	stats.nDisplacements = 0;
	sim.obstructionLimit = 10;
	sim.dx = 1, sim.dy = 1, sim.dz = 1;
	for (int i=0; i<=NBINS; i++){
		stats.rdfx[i] = 0, stats.rdfy[i] = 0, stats.rdfz[i] = 0;
	}
	for (int i=0; i<3; i++) {sim.PBC[i] = true;}
	for (int i=0; i<MAXPART; i++) {part[i].obstructed = 0;}
	fluid.epsilon = fluid.sigma = fluid.rcut = 0;
	fluid.nParts = fluid.dens = fluid.temp = fluid.molarMass = 0;
	fluid.initialE = fluid.finalE = 0;
	fluid.symbol = "";
}
void MC::ReadInputFile(string inFileName){
	int i;
	string line, command, value;
	Operations ops;
	ifstream inFile;
	regex match("^ *([a-z0-9_]+) +?([0-9\\.a-z@, \\-+]+);.*",regex_constants::icase);
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
		else if (command == "equilibriumsteps") {sim.nEquilSteps = stod(value);}
		else if (command == "printevery") {sim.printEvery = stod(value);}
		else if (command == "numberofparticles") {fluid.nParts = stod(value);}
		else if (command == "density") {fluid.dens = stod(value);}
		else if (command == "temperature") {fluid.temp = stod(value);}
		else if (command == "chemicalsymbol") {fluid.symbol = value;}
		else if (command == "molarmass") {fluid.molarMass = stod(value);}
		else if (command == "sigma") {fluid.sigma = stod(value);}
		else if (command == "epsilon") {fluid.epsilon = stod(value);}
		else if (command == "rcut") {fluid.rcut = stod(value);}
	}
	inFile.close();
}
void MC::ComputeBoxSize(void){
	double boxWidth, mass, volume;

	mass = fluid.nParts/na*fluid.molarMass; //g
	volume = mass/fluid.dens; //cm^3
	boxWidth = cbrt(volume); //cm
	boxWidth *= 1e8; //AA
	fluid.boxWidth = boxWidth;
	stats.binWidth = boxWidth/NBINS;
}
void MC::PrintParams(void){
	cout << "Production steps: " << sim.nSteps << endl;
	cout << "Equilibration steps: " << sim.nEquilSteps << endl;
	cout << "Print every n steps: " << 	sim.printEvery << endl;
	cout << "Temperature (K): " << fluid.temp << endl;
	cout << "Molar mass (g/mol): " << fluid.molarMass << endl;
	cout << "sigma (nm): " << fluid.sigma << endl;
	cout << "epsilon (K): " << fluid.epsilon << endl;
	cout << "rcut (nm): " << fluid.rcut << endl;
	cout << "Number of particles: " << fluid.nParts << endl;
	cout << "Density (g/cm^3): " << fluid.dens << endl;
	cout << "Box width (AA): " << fluid.boxWidth << endl;
	cout << endl;
}
void MC::AssignPossitions(void){
	for (int i=1; i<=fluid.nParts; i++){
		part[i].x = Random()*fluid.boxWidth; //AA
		part[i].y = Random()*fluid.boxWidth; //AA
		part[i].z = Random()*fluid.boxWidth; //AA
	}
}
void MC::CreateEXYZ(int step){
	string outFileName = "stats/trajectory.exyz";
	ofstream trajFile;
	double tmp = 0.0;

	if (step == sim.nEquilSteps){trajFile.open(outFileName);}
	else{trajFile.open(outFileName, ios::app);}
	trajFile << fluid.nParts << "\n";
	trajFile << "Lattice=\" " << fluid.boxWidth << " " << tmp << " " << tmp << " ";
	trajFile << tmp << " " << fluid.boxWidth << " " << tmp << " ";
	trajFile << tmp << " " << tmp << " " << fluid.boxWidth << "\" ";
	trajFile << "Properties=species:S:1:pos:R:3 Time=" << 1.*step << "\n";
	for (int i=1; i<=fluid.nParts; i++){
		trajFile << fluid.symbol << "\t";
		trajFile << part[i].x << "\t" << part[i].y << "\t" << part[i].z << "\n"; //AA
	}
}
void MC::RDF(void){
	double dw = stats.binWidth;

	for (int i=1; i<=fluid.nParts; i++){
		stats.rdfx[int(part[i].x/dw)+1]++;
		stats.rdfy[int(part[i].y/dw)+1]++;
		stats.rdfz[int(part[i].z/dw)+1]++;

	}
}
//void MC::PrintRDF(Simulation sim, Stats stats){
	//double MC::pot = 0;
	//FILE *outFile;

	//char outFileName[50];
	//double MC::dw = stats.binWidth;

	//sprintf(outFileName, "temp/rdf_%i.dat", sim.step);
	//outFile = fopen(outFileName, "w");
	//for (int i=1; i<=NBINS; i++){
		//fprintf(outFile, "%f %f", (i-0.5)*dw, redfx[i]/())
	//}
//}
int MC::SelectParticle(void){
	int index;

	do{
		index = int(Random()*fluid.nParts)+1;
	}while (index > fluid.nParts);
	part[0] = part[index]; //Save old position of selected particle.
	fluid.initialE = EnergyOfParticle(0);
	return index;
}
void MC::MoveParticle(int index){
	part[index].x += (Random()-0.5)*sim.dx; //AA
	part[index].y += (Random()-0.5)*sim.dy; //AA
	part[index].z += (Random()-0.5)*sim.dz; //AA
	PBC(index);
	fluid.finalE = EnergyOfParticle(index);
	stats.nDisplacements++;
	AdjustStepSize();
}
void MC::AdjustStepSize(){
	//double ratio;

	//ratio = (1.*stats.acceptance)/(1.*stats.nDisplacements);
	//if (ratio < 0.5) {
		//sim.dx /= 1.1;
		//sim.dy /= 1.1;
		//sim.dz /= 1.1;
	//}else if ((ratio > 0.6) && (ratio <= fluid.boxWidth*0.5/1.1)){
		//sim.dx *= 1.1;
		//sim.dy *= 1.1;
		//sim.dz *= 1.1;
	//}
	sim.dx = sim.dy = sim.dz = fluid.sigma;
	//cout << "New step distance for a ratio "<< ratio << " :" << sim.dx << endl;
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
	double totalE=0;
	unordered_map<string,double> params=GetFluidParams();
	vector<vector<double>> positions=GetPositions();

	totalE = LJ_Energy(index, params, positions);
	return totalE;
}
void MC::ResetParticle(int index){
	part[index].x = Random()*fluid.boxWidth; //AA
	part[index].y = Random()*fluid.boxWidth; //AA
	part[index].z = Random()*fluid.boxWidth; //AA
	part[index].obstructed = 0;
}
void MC::PrintStats(int step){
	cout << "Step: " << step;
	cout << " Acceptance rate: " << stats.acceptance*1./step;
	cout << " Rejection rate: " << stats.rejection*1./step;
	cout << " Obstruction rate: " << stats.obstruction*1./step;
	cout << endl;
}
unordered_map<string,double> MC::GetFluidParams(void){
	unordered_map<string,double> fluidParams;
	unordered_map<string,vector<double>> particles;

	fluidParams["epsilon[K]"] = fluid.epsilon;
	fluidParams["sigma[AA]"] = fluid.sigma;
	fluidParams["rcut[AA]"] = fluid.rcut;
	fluidParams["molarMass[g/mol]"] = fluid.molarMass;
	fluidParams["temperature[K]"] = fluid.temp;
	fluidParams["density[g/mol]"] = fluid.dens;
	fluidParams["numOfParticles"] = (double)fluid.nParts;
	fluidParams["boxWidth[AA]"] = fluid.boxWidth;
	return fluidParams;
}
vector<vector<double>> MC::GetPositions(void){
	vector<vector<double>> positions(fluid.nParts+1, vector<double>(3));

	for (int i=1; i<=fluid.nParts; i++){
		positions[i][0] = part[i].x;
		positions[i][1] = part[i].y;
		positions[i][2] = part[i].z;
	}
	return positions;
}

