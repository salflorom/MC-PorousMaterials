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
	stats.acceptance = stats.rejection = stats.nDisplacements = 0;
	stats.acceptanceVol = stats.rejectionVol = stats.nVolChanges = 0;
	stats.insertion = stats.deletion = stats.nExchanges = 0;
	stats.obstruction = 0;
	stats.widomInsertions = 0;
	stats.widom = 0.;
	fluid.mu = 0.;
	fluid.extPress = 0.;
	sim.dx = 1., sim.dy = 1., sim.dz = 1., sim.dv = 1.;
	sim.nEquilSets = sim.nSets = sim.nStepsPerSet = 0;
	sim.displaceProb = 1.;
	sim.volumeProb = sim.exchangeProb = 0.;
	for (i=0; i<=NBINS; i++) stats.rdf[i] = 0.;
	for (i=0; i<3; i++) {sim.PBC[i] = true;}
	//for (i=0; i<MAXPART; i++) {part[i].obstructed = 0;}
	fluid.epsilon = fluid.sigma = fluid.rcut = 0.;
	fluid.nParts = 0;
	fluid.dens = fluid.temp = fluid.molarMass = 0.;
	fluid.name = fluid.vdwPot = sim.geometry = "";
}
void MC::ResetStats(void){
	stats.acceptance = stats.rejection = stats.nDisplacements = 0;
	stats.acceptanceVol = stats.rejectionVol = 0;
	stats.insertion = stats.deletion = 0;
	stats.widomInsertions = 0;
	stats.widom = 0.;
	stats.obstruction = 0;
	for (int bin=1; bin<=NBINS; bin++) stats.rdf[bin] = 0.;
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
		if (command == "productionsets") {sim.nSets = stod(value);}
		else if (command == "equilibriumsets") {sim.nEquilSets = stod(value);}
		else if (command == "stepsperset") {sim.nStepsPerSet = stod(value);}
		else if (command == "printeverynsets") {sim.printEvery = stod(value);}
		else if (command == "temperature") {fluid.temp = stod(value);}
		else if (command == "externalpressure") {fluid.extPress = stod(value);  // Pa
		}
		else if (command == "geometry") {sim.geometry = value;}
		else if (command == "name") {fluid.name = value;}
		else if (command == "numberofparticles") {fluid.nParts = stoi(value);}
		else if (command == "density") {fluid.dens = stod(value);}
		else if (command == "molarmass") {fluid.molarMass = stod(value);}
		else if (command == "vdwpotential") {fluid.vdwPot = value;}
		else if (command == "sigma") {fluid.sigma = stod(value);}
		else if (command == "epsilon") {fluid.epsilon = stod(value);}
		else if (command == "rcut") {fluid.rcut = stod(value);}
		else if (command == "chemicalpotential") {fluid.mu = stod(value);}
		else if (command == "displacementprobability") {sim.displaceProb = stod(value);}
		else if (command == "changevolumeprobability") {sim.volumeProb = stod(value);}
		else if (command == "exchangeprobability") {sim.exchangeProb = stod(value);}

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
	stats.binWidth = fluid.rcut/NBINS;
	sim.dx = sim.dy = sim.dz = fluid.sigma; //AA. Step size set to sigma.
	sim.dv = fluid.sigma*fluid.sigma*fluid.sigma; //AA^3. Volume step size set to sigma^3.
}
double* MC::GetMCMoveProbabilities(void){
	static double cumulativeMoveProbs[3]; //Size equals num. of MC moves.
	double totalProb;

	totalProb = sim.displaceProb + sim.exchangeProb + sim.volumeProb;
	cumulativeMoveProbs[0] = sim.displaceProb;
	cumulativeMoveProbs[1] = cumulativeMoveProbs[0] + sim.volumeProb;
	cumulativeMoveProbs[2] = cumulativeMoveProbs[1] + sim.exchangeProb;
	for (int i=0; i<3; i++) cumulativeMoveProbs[i] /= totalProb; //Normalized cumulative prob.
	return cumulativeMoveProbs;
}
void MC::PrintParams(void){
	cout << "Production sets: " << sim.nSets << endl;
	cout << "Equilibration sets: " << sim.nEquilSets << endl;
	cout << "Steps per set: " << sim.nStepsPerSet << endl;
	cout << "Print every n sets: " << 	sim.printEvery << endl;
	cout << "Temperature (K): " << fluid.temp << endl;
	if (sim.volumeProb > 0) cout << "External pressure (Pa): " << fluid.extPress << endl;
	if (sim.exchangeProb > 0) cout << "mu (K): " << fluid.mu << endl;
	cout << "Fluid name: " << fluid.name << endl;
	cout << "Molar mass (g/mol): " << fluid.molarMass << endl;
	cout << "VdW potential: " << fluid.vdwPot << endl;
	cout << "sigma (AA): " << fluid.sigma << endl;
	cout << "epsilon (K): " << fluid.epsilon << endl;
	cout << "rcut (AA): " << fluid.rcut << endl;
	cout << "Number of particles: " << fluid.nParts << endl;
	cout << "Density (g/cm^3): " << fluid.dens << endl;
	cout << "Initial box width (AA): " << fluid.boxWidth << endl;
	cout << "Displacement probability: " << sim.displaceProb << endl;
	cout << "Change volume probability: " << sim.volumeProb << endl;
	cout << "Exchange probability: " << sim.exchangeProb << endl;
	cout << endl;
}
void MC::InitialConfig(void){
	for (int i=1; i<=fluid.nParts; i++){
		part[i].x = Random()*fluid.boxWidth;
		part[i].y = Random()*fluid.boxWidth;
		part[i].z = Random()*fluid.boxWidth;
	}
}
void MC::MoveParticle(void){
	int index;
	double initialEnergy, finalEnergy, argument;

	// Select particle.
	do{
		index = int(Random()*fluid.nParts)+1;
	}while (index > fluid.nParts);
	part[0] = part[index]; //Save old position of selected particle.
	initialEnergy = EnergyOfParticle(index);
	// Move particle.
	part[index].x += (Random()-0.5)*sim.dx; //AA
	part[index].y += (Random()-0.5)*sim.dy; //AA
	part[index].z += (Random()-0.5)*sim.dz; //AA
	PBC(index);
	finalEnergy = EnergyOfParticle(index);
	stats.nDisplacements++;
	// Metropolis.
	argument = -(finalEnergy-initialEnergy)/fluid.temp;
	if (exp(argument) > Random()){
		stats.acceptance++;
		part[index].obstructed = 0;
	}else{
		part[index] = part[0];
		part[index].obstructed++;
		stats.rejection++;
	}
}
// PBC for a box with the origin at the lower left vertex.
void MC::PBC(int index){
	part[index].x -= floor(part[index].x/fluid.boxWidth) * fluid.boxWidth; //AA
	part[index].y -= floor(part[index].y/fluid.boxWidth) * fluid.boxWidth; //AA
	part[index].z -= floor(part[index].z/fluid.boxWidth) * fluid.boxWidth; //AA
}
void MC::ChangeVolume(void){
	double initEnergy, finalEnergy, deltaEnergy;
	double initVol, finalVol, logFinalVol, deltaVol;
	double initBoxWidth, finalBoxWidth;
	double argument;
	double extPress;

	stats.nVolChanges++;
	extPress = fluid.extPress/kb*1e-30; // K/AA^3
	initEnergy = SystemEnergy();
	initBoxWidth = fluid.boxWidth;
	initVol = fluid.boxWidth*fluid.boxWidth*fluid.boxWidth;
	logFinalVol = log(initVol) + (Random()-0.5)*sim.dv; //Perform volume change step.
	finalVol = exp(logFinalVol);
	deltaVol = finalVol-initVol;
	finalBoxWidth = cbrt(finalVol);
	for (int i=1; i<=fluid.nParts; i++){ //Rescale center of mass.
		part[i].x *= finalBoxWidth/initBoxWidth;
		part[i].y *= finalBoxWidth/initBoxWidth;
		part[i].z *= finalBoxWidth/initBoxWidth;
	}
	fluid.boxWidth = finalBoxWidth;
	finalEnergy = SystemEnergy();
	deltaEnergy = finalEnergy-initEnergy;
	argument = -(deltaEnergy + extPress*deltaVol - (fluid.nParts+1)*log(finalVol/initVol)*fluid.temp)/fluid.temp;
	if ((Random() > exp(argument)) || (finalBoxWidth < 2*fluid.rcut)){ //Rejected change-volume trial.
		for (int i=1; i<=fluid.nParts; i++){ //Rescale center of mass to initial.
			part[i].x *= initBoxWidth/finalBoxWidth;
			part[i].y *= initBoxWidth/finalBoxWidth;
			part[i].z *= initBoxWidth/finalBoxWidth;
		}
		fluid.boxWidth = initBoxWidth;
		stats.rejectionVol++;
	}else stats.acceptanceVol++;
}
void MC::ExchangeParticle(void){
	double partMass, thermalWL, zz, energy, argument, volume;
	double randNum;
	int index;

	partMass = fluid.molarMass/na*1e3; // kg/particle
	//rho = fluid.dens/fluid.molarMass*na*1e6; // particle/m^3
	thermalWL = plank/sqrt(2*pi*partMass*kb*fluid.temp); // m^3
	zz=exp(fluid.mu/fluid.temp)/(thermalWL*thermalWL*thermalWL); // m^-3
	volume = fluid.boxWidth*fluid.boxWidth*fluid.boxWidth*1e-30; //m^3
	randNum = Random();
	if ((randNum < 0.5 ) && (fluid.nParts > 0)){ // Try removing particle.
		// Select particle.
		do{
			index = int(Random()*fluid.nParts)+1;
		}while (index > fluid.nParts);
		energy = EnergyOfParticle(index); //K
		argument = fluid.nParts*exp(energy/fluid.temp)/(zz*volume); // Motropolis (for removing particle).
		if (Random() < argument){ // Accepted: Remove particle.
			part[index] = part[fluid.nParts];
			fluid.nParts--;
			stats.insertion++;
			stats.nExchanges++;
		}
	}else if (randNum > 0.5){ // Try inserting particle.
		// Insert particle at random position.
		index = fluid.nParts+1;
		part[index].x = Random()*fluid.boxWidth; //AA
		part[index].y = Random()*fluid.boxWidth; //AA
		part[index].z = Random()*fluid.boxWidth; //AA
		energy = EnergyOfParticle(index); //K
		argument = zz*volume*exp(-energy/fluid.temp)/(fluid.nParts+1); // Metropolis (for inserting particle).
		if (Random() < argument){
			fluid.nParts++; // Accepted: Insert particle.
			stats.deletion++;
			stats.nExchanges++;
		}
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
void MC::ComputeWidom(void){
	double volume;
	if (sim.exchangeProb == 0){
		// Insert virtual particle.
		part[MAXPART-1].x = Random()*fluid.boxWidth;
		part[MAXPART-1].y = Random()*fluid.boxWidth;
		part[MAXPART-1].z = Random()*fluid.boxWidth;
		stats.widomInsertions++;
		if (sim.volumeProb > 0){
			volume = fluid.boxWidth*fluid.boxWidth*fluid.boxWidth;
			stats.widom += volume*exp(-EnergyOfParticle(MAXPART-1)/fluid.temp);
		}else stats.widom += exp(-EnergyOfParticle(MAXPART-1)/fluid.temp);
	}
}
void MC::ComputeChemicalPotential(void){
	double thermalWL, partMass, volume, muIdeal, muExcess, insertionParam;
	double extPress;

	if (sim.exchangeProb == 0){
		partMass = fluid.molarMass/na*1e-3; // kg/particle
		thermalWL = plank/sqrt(2*pi*partMass*kb*fluid.temp)*1e10; //AA
		volume = fluid.boxWidth*fluid.boxWidth*fluid.boxWidth; //AA^3
		insertionParam = stats.widom/stats.widomInsertions;
		if (sim.volumeProb > 0){
			extPress = fluid.extPress/kb*1e-30; // K/AA^3
			muIdeal = fluid.temp*log(thermalWL*thermalWL*thermalWL*extPress/fluid.temp); //K
			muExcess = -fluid.temp*log(insertionParam*extPress/(fluid.temp*(fluid.nParts+1))); //K
		}else{
			muIdeal = fluid.temp*log(thermalWL*thermalWL*thermalWL*(fluid.nParts+1)/volume); //K
			muExcess = -fluid.temp*log(insertionParam); //K
		}
		fluid.mu = muIdeal+muExcess; //K
	}
}
void MC::ResetParticle(int index){
	part[index].x = Random()*fluid.boxWidth; //AA
	part[index].y = Random()*fluid.boxWidth; //AA
	part[index].z = Random()*fluid.boxWidth; //AA
}
void MC::ComputeRDF(void){
	int bin;
	double deltaR;
	Operations ops;

	deltaR = fluid.rcut/NBINS;
	for (int i=1; i<=int(fluid.nParts)-1; i++){
		for (int j=i+1; j<=int(fluid.nParts); j++){
			bin = int(NeighDistance(i, j)/deltaR)+1;
			if (bin <= NBINS) stats.rdf[bin] += 2; // Takes into account both i->j and j->i.
		}
	}
}

// Check PrintRDF!!
void MC::PrintRDF(int set){
	double constant, deltaR, lowR, highR, rho, dn_ideal, gofr[NBINS+1];
	Operations ops;
	ofstream rdfFile;
	string outFileName="stats/rdf.dat";

	deltaR = fluid.rcut/NBINS;
	rho = fluid.nParts/ops.Pow(fluid.boxWidth,3);
	constant = (4*pi*rho/3);
	for (int bin=1; bin<=NBINS; bin++){
		gofr[bin] = stats.rdf[bin]/(1.*fluid.nParts*sim.nStepsPerSet);
		lowR = bin*deltaR;
		highR = lowR + deltaR;
		dn_ideal = constant*(ops.Pow(highR,3)-ops.Pow(lowR,3));
		gofr[bin] /= dn_ideal;
	}
	//Write into the file.
	if (ops.FileExists(outFileName)) rdfFile.open(outFileName, ios::app);
	else rdfFile.open(outFileName);
	rdfFile << "nBins; set: "<< NBINS << "; " << set << endl;
	rdfFile << "r(AA)\tg(r)" << endl;
	for (int bin=1; bin<=NBINS; bin++){
		rdfFile << (bin+0.5)*deltaR << "\t" << gofr[bin] << endl;
	}
	rdfFile.close();
}
void MC::PrintStats(int set){
	cout << "Set: " << set << endl;
	if (stats.nDisplacements > 0){
		cout << "AcceptDispRatio; RegectDispRatio:\t";
		cout << stats.acceptance*1./stats.nDisplacements << "; ";
		cout << stats.rejection*1./stats.nDisplacements << endl;
	}
	if (sim.volumeProb > 0){
		cout << "AcceptVolRatio; RejectVolRatio:\t";
		cout << stats.acceptanceVol*1./stats.nVolChanges << "; ";
		cout << stats.rejectionVol*1./stats.nVolChanges << endl;
	}
	if (sim.exchangeProb > 0){
		cout << "InsertRatio; DeletRatio:\t";
		cout << stats.insertion*1./(stats.nExchanges+1.) << "; ";
		cout << stats.deletion*1./(stats.nExchanges+1.) << endl;
	}
	if (fluid.nParts > 0){
		cout << "NumPartcles; Energy/Particle:\t";
		cout << fluid.nParts << "; ";
		cout << SystemEnergy()/fluid.nParts << endl;
	}
	cout << endl;
}
void MC::CreateEXYZ(int set){
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
	trajFile << "Properties=species:S:1:pos:R:3 Time=" << 1.*set << "\n";
	for (int i=1; i<=fluid.nParts; i++){
		trajFile << fluid.name << "\t";
		trajFile << part[i].x << "\t" << part[i].y << "\t" << part[i].z << "\n"; //AA
	}
	trajFile.close();
}
void MC::CreateLogFile(int set){
	Operations ops;
	ofstream logFile;
	string outFileName="stats/simulation.log";
	double energyPerPart, density, volume;

	density = fluid.nParts/(fluid.boxWidth*fluid.boxWidth*fluid.boxWidth);
	density *= fluid.molarMass/na*1e24; // g/cm^3
	volume = fluid.boxWidth*fluid.boxWidth*fluid.boxWidth; //AA
	if (!ops.FileExists(outFileName)){
		logFile.open(outFileName);
		logFile << "Set\tNParts\tDens[g/cm^3]\t";
		logFile << "BoxWidth[AA]\tVolume[AA^3]\tE/part[K]\tmu[K]\n";
	}else{
		energyPerPart = SystemEnergy()/fluid.nParts;
		logFile.open(outFileName, ios::app);
		logFile << set << "\t";
		logFile << fluid.nParts << "\t";
		logFile << density << "\t";
		logFile << fluid.boxWidth << "\t";
		logFile << volume << "\t";
		logFile << energyPerPart << "\t";
		logFile	<< fluid.mu << "\n";
	}
	logFile.close();
}
double MC::NeighDistance(int i, int j){
	double dist, dx, dy, dz;
	double boxWidth=fluid.boxWidth;

	dx = part[j].x - part[i].x; //AA
	dy = part[j].y - part[i].y; //AA.
	dz = part[j].z - part[i].z; //AA.
	dx -= round(dx/boxWidth) * boxWidth; //AA
	dy -= round(dy/boxWidth) * boxWidth; //AA
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
		phi = (aValues[6] + bValues[6]*(rho-rhoIntervals[5]) + cValues[6]*(rho-rhoIntervals[5])*(rho-rhoIntervals[5])) * (2*rho/rhoIntervals[5]-(rho/rhoIntervals[5])*(rho/rhoIntervals[5]));
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

