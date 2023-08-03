//Author: Santiago A. Flores Roman

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <regex>
#include <cctype>
#include <iomanip>
#include <climits>

#include "operations.h"
#include "MC.h"

using namespace std;

MC::MC(void){
	int i;
	stats.acceptance = stats.rejection = stats.nDisplacements = 0;
	stats.acceptanceVol = stats.rejectionVol = stats.nVolChanges = 0;
	stats.acceptInsertion = stats.rejectInsertion = 0;
	stats.acceptDeletion = stats.rejectDeletion = stats.nExchanges = 0;
	stats.widomInsertions = 0;
	stats.widom = 0.;
	sim.dx = 1., sim.dy = 1., sim.dz = 1., sim.dv = 1.;
	sim.nEquilSets = sim.nSets = sim.nStepsPerSet = 0;
	sim.displaceProb = 1.;
	sim.volumeProb = sim.exchangeProb = 0.;
	sim.systemEnergy = 0.;
	for (i=0; i<=NBINS; i++) stats.rdf[i] = 0.;
	for (i=0; i<3; i++) {pore.PBC[i] = true;}
	fluid.epsilon = fluid.sigma = fluid.rcut = 0.;
	fluid.nParts = 0;
	fluid.temp = fluid.molarMass = fluid.dens = 0.;
	fluid.name = fluid.vdwPot = "";
	fluid.ffEnergy = fluid.mu = 0.;
	fluid.extPress = 0.;
	pore.sfEnergy = 0.;
	for (i=0; i<3; i++) pore.boxWidth[i] = 0;
	pore.name = pore.sfPot = "";
	pore.geometry = "bulk";
}
void MC::ResetStats(void){
	stats.acceptance = stats.rejection = stats.nDisplacements = 0;
	stats.acceptanceVol = stats.rejectionVol = 0;
	stats.acceptInsertion = stats.rejectInsertion = 0;
	stats.acceptDeletion = stats.rejectDeletion = stats.nExchanges = 0;
	stats.widomInsertions = 0;
	stats.widom = 0.;
	for (int bin=1; bin<=NBINS; bin++) stats.rdf[bin] = 0.;
}
void MC::ReadInputFile(string inFileName){
	double volume;
	double density = 0;
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
		if (command == "productionsets") sim.nSets = stod(value);
		else if (command == "equilibriumsets") sim.nEquilSets = stod(value);
		else if (command == "stepsperset") sim.nStepsPerSet = stod(value);
		else if (command == "printeverynsets") sim.printEvery = stod(value);
		else if (command == "temperature") fluid.temp = stod(value); // K
		else if (command == "externalpressure") fluid.extPress = stod(value);  // Pa
		else if (command == "chemicalpotential") fluid.mu = stod(value); // K
		else if (command == "name") fluid.name = value;
		else if (command == "numberofparticles") fluid.nParts = stoi(value);
		else if (command == "molarmass") fluid.molarMass = stod(value); // g/mol
		else if (command == "vdwpotential") fluid.vdwPot = value;
		else if (command == "sigma") fluid.sigma = stod(value); // AA
		else if (command == "epsilon") fluid.epsilon = stod(value); // K
		else if (command == "rcut") fluid.rcut = stod(value); // AA
		else if (command == "frameworkname") pore.name = value;
		else if (command == "geometry") pore.geometry = value;
		else if (command == "displacementprobability") sim.displaceProb = stod(value);
		else if (command == "changevolumeprobability") sim.volumeProb = stod(value);
		else if (command == "exchangeprobability") sim.exchangeProb = stod(value);
		else if (command == "density") density = stod(value); // g/cm^3
		else if (command == "poresize") pore.boxWidth[2] = stod(value); // AA. boxWidth[2] (along z axis) will always keep the pore size.
	}
	inFile.close();

	if (density > 0){
		pore.boxWidth[0] = pore.boxWidth[1] = 2*fluid.rcut; // Temporal box widths.
		density *= na/fluid.molarMass*1e-24; // AA^-3
		volume = fluid.nParts*1./density; // AA^3
		pore.boxWidth[2] = ComputeBoxWidth(volume); // AA
	}
	// Get remaining box widths according to the geometry given.
	if (pore.geometry == "sphere") {pore.boxWidth[0] = pore.boxWidth[1] = pore.boxWidth[2];
	}else if (pore.geometry == "cylinder"){
		pore.boxWidth[0] = 2*fluid.rcut; // PBC along x axis.
		pore.boxWidth[1] = pore.boxWidth[2];
	}else if (pore.geometry == "slit") {pore.boxWidth[0] = pore.boxWidth[1] = 2*fluid.rcut; // PBC along xy plane.
	}else pore.boxWidth[0] = pore.boxWidth[1] = pore.boxWidth[2]; // Bulk phase.
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
	double density;

	density = fluid.nParts/ComputeVolume(); //AA^-3
	density *= fluid.molarMass/na*1e24; // g/cm^3
	cout << "Production sets: " << sim.nSets << endl;
	cout << "Equilibration sets: " << sim.nEquilSets << endl;
	cout << "Steps per set: " << sim.nStepsPerSet << endl;
	cout << "Print every n sets: " << 	sim.printEvery << endl;
	cout << "Temperature (K): " << fluid.temp << endl;
	if (sim.volumeProb > 0) cout << "External pressure (Pa): " << fluid.extPress << endl;
	if (sim.exchangeProb > 0) cout << "Chemical potential (K): " << fluid.mu << endl;
	cout << "Fluid name: " << fluid.name << endl;
	cout << "Molar mass (g/mol): " << fluid.molarMass << endl;
	cout << "VdW potential: " << fluid.vdwPot << endl;
	cout << "sigma (AA): " << fluid.sigma << endl;
	cout << "epsilon (K): " << fluid.epsilon << endl;
	cout << "rcut (AA): " << fluid.rcut << endl;
	cout << "Number of particles: " << fluid.nParts << endl;
	cout << "Initial density (g/cm^3): " << density << endl;
	if (pore.name != "") cout << "Framework name: " << pore.name << endl;
	cout << "Geometry: " << pore.geometry << endl;
	cout << "Initial box width or pore size (AA): " << pore.boxWidth[2] << endl;
	cout << "Displacement probability: " << sim.displaceProb << endl;
	cout << "Change volume probability: " << sim.volumeProb << endl;
	cout << "Exchange probability: " << sim.exchangeProb << endl;
	cout << endl;

	if ((pore.boxWidth[2] < 2*fluid.rcut) && (pore.geometry == "bulk")){
		cout << "Error: Box length must be larger than twice the cut-off radius.";
		cout << endl;
		exit(EXIT_FAILURE);
	}
}
void MC::InsertParticle(int index){
	double radius, theta, phi; // Spherical geometry.
	double rho; // Cylindrical geometry.

	if (pore.geometry == "sphere"){
		radius = Random()*pore.boxWidth[2]*0.5;
		theta = Random()*pi;
		phi = Random()*2*pi;
		part[index].x = radius*sin(theta)*cos(phi);
		part[index].y = radius*sin(theta)*sin(phi);
		part[index].z = radius*cos(theta);
		// Frame is at (0,0,0). Moving part. to box center.
		part[index].x += 0.5*pore.boxWidth[0];
		part[index].y += 0.5*pore.boxWidth[1];
		part[index].z += 0.5*pore.boxWidth[2];
	}else if (pore.geometry == "cylinder"){
		rho = Random()*pore.boxWidth[2]*0.5;
		phi = Random()*2*pi;
		part[index].x = (Random()-0.5)*pore.boxWidth[0];
		part[index].y = rho*sin(phi);
		part[index].z = rho*cos(phi);
		// Frame is at (0,0,0). Moving part. to box center.
		part[index].x += 0.5*pore.boxWidth[0];
		part[index].y += 0.5*pore.boxWidth[1];
		part[index].z += 0.5*pore.boxWidth[2];
	}else if (pore.geometry == "slit"){
		part[index].x = Random()*pore.boxWidth[0];
		part[index].y = Random()*pore.boxWidth[1];
		part[index].z = Random()*pore.boxWidth[2];
	}else{
		part[index].x = Random()*pore.boxWidth[0];
		part[index].y = Random()*pore.boxWidth[1];
		part[index].z = Random()*pore.boxWidth[2];
	}
}
void MC::InitialConfig(void){
	for (int i=1; i<=fluid.nParts; i++) InsertParticle(i);
	if (pore.geometry == "sphere"){
		pore.PBC[0] = false;
		pore.PBC[1] = false;
		pore.PBC[2] = false;
	}else if (pore.geometry == "cylinder"){
		pore.PBC[0] = true;
		pore.PBC[1] = false;
		pore.PBC[2] = false;
	}else if (pore.geometry == "slit"){
		pore.PBC[0] = true;
		pore.PBC[1] = true;
		pore.PBC[2] = false;
	}else if (pore.geometry == "bulk"){
		pore.PBC[0] = true;
		pore.PBC[1] = true;
		pore.PBC[2] = true;
	}
}
void MC::MinimizeEnergy(void){
	int initMoves = 200;

	sim.systemEnergy = SystemEnergy()[0];
	cout << "Minimizing initial configuration" << endl;
	cout << "Energy/part. before minimization: ";
	cout << sim.systemEnergy/fluid.nParts << " K"<< endl;
	for (int i=1; i<=initMoves; i++){
		for (int j=1; j<=fluid.nParts; j++) MoveParticle();
	}
	sim.systemEnergy = SystemEnergy()[0];
	fluid.ffEnergy = SystemEnergy()[1];
	pore.sfEnergy = SystemEnergy()[2];
	cout << "Energy/part. after minimization (";
	cout << fluid.nParts*initMoves << " moves): ";
	cout << sim.systemEnergy/fluid.nParts << " K" << endl;
	cout << endl;
}
void MC::MoveParticle(void){
	int index;
	double initialEnergy, finalEnergy, deltaEnergy, argument;
	double initialffEnergy, finalffEnergy;
	double initialsfEnergy, finalsfEnergy;

	stats.nDisplacements++;
	// Select particle.
	index = int(Random()*fluid.nParts)+1;
	part[0] = part[index]; //Save old position of selected particle.
	EnergyOfParticle(index);
	initialEnergy = sim.deltaParticleEnergy;
	initialffEnergy = fluid.deltaffEnergy;
	initialsfEnergy = pore.deltasfEnergy;
	// Move particle.
	part[index].x += (Random()-0.5)*sim.dx; //AA
	part[index].y += (Random()-0.5)*sim.dy; //AA
	part[index].z += (Random()-0.5)*sim.dz; //AA
	PBC(index);
	EnergyOfParticle(index);
	finalEnergy = sim.deltaParticleEnergy;
	finalffEnergy = fluid.deltaffEnergy;
	finalsfEnergy = pore.deltasfEnergy;
	// Metropolis.
	deltaEnergy = finalEnergy-initialEnergy;
	argument = -deltaEnergy/fluid.temp;
	if (exp(argument) > Random()){
		stats.acceptance++;
		sim.systemEnergy += deltaEnergy;
		fluid.ffEnergy += finalffEnergy-initialffEnergy;
		pore.sfEnergy += finalsfEnergy-initialsfEnergy;
	}
	else{
		part[index] = part[0];
		stats.rejection++;
	}
}
// PBC for a box with the origin at the lower left vertex.
void MC::PBC(int index){
	if (pore.PBC[0]) part[index].x -= floor(part[index].x/pore.boxWidth[0]) * pore.boxWidth[0]; //AA
	if (pore.PBC[1]) part[index].y -= floor(part[index].y/pore.boxWidth[1]) * pore.boxWidth[1]; //AA
	if (pore.PBC[2]) part[index].z -= floor(part[index].z/pore.boxWidth[2]) * pore.boxWidth[2]; //AA
}
double MC::ComputeVolume(void){
	Operations ops;

	if (pore.geometry == "sphere") return (4*pi/3.)*ops.Pow(pore.boxWidth[2]*0.5,3);
	else if (pore.geometry == "cylinder") return pi*ops.Pow(pore.boxWidth[2]*0.5,2)*pore.boxWidth[0];
	else if (pore.geometry == "slit") return pore.boxWidth[0]*pore.boxWidth[1]*pore.boxWidth[2];
	else return ops.Pow(pore.boxWidth[2],3);
}
// Computes the length of the box along the z axis.
double MC::ComputeBoxWidth(double volume){
	if (pore.geometry == "sphere") return 2*cbrt(3*volume/(4*pi));
	else if (pore.geometry == "cylinder") return 2*sqrt(volume/(pi*pore.boxWidth[0]));
	else if (pore.geometry == "slit") return volume/(pore.boxWidth[0]*pore.boxWidth[1]);
	else return cbrt(volume); //Bulk phase.
}
void MC::RescaleCenterOfMass(double initBoxWidth, double finalBoxWidth, bool rescale){
	int i;

	if (rescale){
		if (pore.geometry == "sphere"){
			for (i=1; i<=fluid.nParts; i++){
				part[i].x *= finalBoxWidth/initBoxWidth;
				part[i].y *= finalBoxWidth/initBoxWidth;
				part[i].z *= finalBoxWidth/initBoxWidth;
			}
		}else if (pore.geometry == "cylinder"){
			for (i=1; i<=fluid.nParts; i++){
				part[i].y *= finalBoxWidth/initBoxWidth;
				part[i].z *= finalBoxWidth/initBoxWidth;
			}
		}else if (pore.geometry == "slit"){
			for (i=1; i<=fluid.nParts; i++) part[i].z *= finalBoxWidth/initBoxWidth;
		}else{ //Bulk phase.
			for (i=1; i<=fluid.nParts; i++){
				part[i].x *= finalBoxWidth/initBoxWidth;
				part[i].y *= finalBoxWidth/initBoxWidth;
				part[i].z *= finalBoxWidth/initBoxWidth;
			}
		}
	}else{
		if (pore.geometry == "sphere"){
			for (i=1; i<=fluid.nParts; i++){
				part[i].x *= initBoxWidth/finalBoxWidth;
				part[i].y *= initBoxWidth/finalBoxWidth;
				part[i].z *= initBoxWidth/finalBoxWidth;
			}
		}else if (pore.geometry == "cylinder"){
			for (i=1; i<=fluid.nParts; i++){
				part[i].y *= initBoxWidth/finalBoxWidth;
				part[i].z *= initBoxWidth/finalBoxWidth;
			}
		}else if (pore.geometry == "slit"){
			for (i=1; i<=fluid.nParts; i++) part[i].z *= initBoxWidth/finalBoxWidth;
		}else{ //Bulk phase.
			for (i=1; i<=fluid.nParts; i++){
				part[i].x *= initBoxWidth/finalBoxWidth;
				part[i].y *= initBoxWidth/finalBoxWidth;
				part[i].z *= initBoxWidth/finalBoxWidth;
			}
		}
	}
}
void MC::ChangeVolume(void){
	double initEnergy, finalEnergy, deltaEnergy;
	double finalffEnergy, finalsfEnergy;
	double initVol, finalVol, logFinalVol, deltaVol;
	double initBoxWidth, finalBoxWidth;
	double argument;
	double extPress;

	stats.nVolChanges++;
	extPress = fluid.extPress/kb*1e-30; // K/AA^3
	// Record current config.
	initEnergy = sim.systemEnergy;
	initBoxWidth = pore.boxWidth[2];
	initVol = ComputeVolume();
	logFinalVol = log(initVol) + (Random()-0.5)*sim.dv; //Perform volume change step.
	finalVol = exp(logFinalVol);
	deltaVol = finalVol-initVol;
	finalBoxWidth = ComputeBoxWidth(finalVol);
	RescaleCenterOfMass(initBoxWidth, finalBoxWidth, true); //Rescale to trial config.
	// Set box sizes according to trial volume.
	if (pore.geometry == "sphere") pore.boxWidth[0] = pore.boxWidth[1] = pore.boxWidth[2] = finalBoxWidth;
	else if (pore.geometry == "cylinder") pore.boxWidth[1] = pore.boxWidth[2] = finalBoxWidth;
	else if (pore.geometry == "slit") pore.boxWidth[2] = finalBoxWidth;
	else pore.boxWidth[0] = pore.boxWidth[1] = pore.boxWidth[2] = finalBoxWidth;
	// Compute trial energy.
	finalEnergy = SystemEnergy()[0];
	finalffEnergy = SystemEnergy()[1];
	finalsfEnergy = SystemEnergy()[2];
	deltaEnergy = finalEnergy-initEnergy;
	argument = -(deltaEnergy + extPress*deltaVol - (fluid.nParts+1)*log(finalVol/initVol)*fluid.temp)/fluid.temp;
	if (Random() > exp(argument)){ //Rejected trial move.
		RescaleCenterOfMass(initBoxWidth, finalBoxWidth, false); //Rescale to previous config.
		// Redo box sizes.
		if (pore.geometry == "sphere") pore.boxWidth[0] = pore.boxWidth[1] = pore.boxWidth[2] = initBoxWidth;
		else if (pore.geometry == "cylinder") pore.boxWidth[1] = pore.boxWidth[2] = initBoxWidth;
		else if (pore.geometry == "slit") pore.boxWidth[2] = initBoxWidth;
		else pore.boxWidth[0] = pore.boxWidth[1] = pore.boxWidth[2] = initBoxWidth;
		stats.rejectionVol++;
	}else{ // Accepted trial move. Keep trial config.
		stats.acceptanceVol++;
		sim.systemEnergy = finalEnergy;
		fluid.ffEnergy = finalffEnergy;
		pore.sfEnergy = finalsfEnergy;
	}
}
void MC::ExchangeParticle(void){
	Operations ops;
	int index;
	double particleMass, thermalWL, zz, energy, volume;
	double argument=0;

	particleMass = fluid.molarMass/na*1e3; // kg/particle
	thermalWL = plank/sqrt(2*pi*particleMass*kb*fluid.temp)*1e10; // AA
	zz=exp(fluid.mu/fluid.temp)/ops.Pow(thermalWL,3); // AA^-3
	volume = ComputeVolume(); //AA^3
	if (Random() < 0.5) { // Try inserting particle.
		stats.nExchanges++;
		// Insert particle at random position.
		index = fluid.nParts+1;
		InsertParticle(index);
		EnergyOfParticle(index); //K
		energy = sim.deltaParticleEnergy;
		argument = zz*volume*exp(-energy/fluid.temp)/(fluid.nParts+1); // Acceptance criterion (for inserting particle).
		if (Random() < argument){
			stats.acceptInsertion++; // Accepted: Insert particle.
			fluid.nParts++;
			sim.systemEnergy += energy;
			fluid.ffEnergy += fluid.deltaffEnergy;
			pore.sfEnergy += pore.deltasfEnergy;
		}else stats.rejectInsertion++;
	}else if (fluid.nParts > 0){ // Try removing particle (only if there are particles in the box).
		stats.nExchanges++;
		// Select random particle.
		index = int(Random()*fluid.nParts)+1;
		EnergyOfParticle(index); //K
		energy = sim.deltaParticleEnergy;
		argument = fluid.nParts*exp(energy/fluid.temp)/(zz*volume); // Acceptance criterion (for removing particle).
		if (Random() < argument){ // Accepted: Remove particle.
			stats.acceptDeletion++;
			for (int i=index; i<fluid.nParts; i++) part[i] = part[i+1]; //Fill in the gap left by removal.
			fluid.nParts--;
			sim.systemEnergy -= energy;
			fluid.ffEnergy -= fluid.deltaffEnergy;
			pore.sfEnergy -= pore.deltasfEnergy;
		}else stats.rejectDeletion++;
	}
}
void MC::EnergyOfParticle(int index){
	Operations ops;
	double xPos, yPos, zPos, sqrPoreRadius;
	double sqrDistFromBoxCenter = 0.;

	sim.deltaParticleEnergy = fluid.deltaffEnergy = pore.deltasfEnergy = 0.;
	// Particle must not be out of boundaries.
	sqrPoreRadius = pore.boxWidth[2]*pore.boxWidth[2]*0.25;
	xPos = part[index].x - 0.5*pore.boxWidth[0]; //Move positions virtually to center of box.
	yPos = part[index].y - 0.5*pore.boxWidth[1];
	zPos = part[index].z - 0.5*pore.boxWidth[2];
	if (pore.geometry == "sphere") sqrDistFromBoxCenter = xPos*xPos + yPos*yPos + zPos*zPos;
	else if (pore.geometry == "cylinder") sqrDistFromBoxCenter = yPos*yPos + zPos*zPos;
	else if (pore.geometry == "slit") sqrDistFromBoxCenter = zPos*zPos;
	if (sqrDistFromBoxCenter >= sqrPoreRadius) pore.deltasfEnergy = INT_MAX;
	//Solid-Fluid energy
	//	To implement ...
	//Fluid-Fluid energy
	if (fluid.vdwPot == "lj"){ //Lennard-Jones potential.
		for (int i=1; i<=fluid.nParts; i++) fluid.deltaffEnergy += LJ_Energy(index,i);
		fluid.deltaffEnergy *= 0.5;
	}else if (fluid.vdwPot == "eam_ga"){ //EAM potential for Ga.
		fluid.deltaffEnergy += EAMGA_Energy(index);
	}
	sim.deltaParticleEnergy = fluid.deltaffEnergy + pore.deltasfEnergy;
}
double* MC::SystemEnergy(void){
	double energy=0, ffEnergy=0, sfEnergy=0;
	static double energies[3];

	for (int i=1; i<=fluid.nParts; i++){
		EnergyOfParticle(i);
		energy += sim.deltaParticleEnergy;
		ffEnergy += fluid.deltaffEnergy;
		sfEnergy += pore.deltasfEnergy;
	}
	energies[0] = energy;
	energies[1] = ffEnergy;
	energies[2] = sfEnergy;
	return energies;
}
void MC::ComputeWidom(void){
	Operations ops;
	double volume, energy;

	if (sim.exchangeProb == 0){
		InsertParticle(MAXPART-1); // Insert virtual particle.
		stats.widomInsertions++;
		EnergyOfParticle(MAXPART-1);
		energy = sim.deltaParticleEnergy;
		if (sim.volumeProb > 0){
			volume = ComputeVolume();
			stats.widom += volume*exp(-energy/fluid.temp);
		}else stats.widom += exp(-energy/fluid.temp);
	}
}
void MC::ComputeChemicalPotential(void){
	Operations ops;
	double thermalWL, partMass, volume, muIdeal, muExcess, insertionParam;
	double extPress;

	if (sim.exchangeProb == 0){
		partMass = fluid.molarMass/na*1e-3; // kg/particle
		thermalWL = plank/sqrt(2*pi*partMass*kb*fluid.temp)*1e10; //AA
		volume = ComputeVolume(); //AA^3
		insertionParam = stats.widom/stats.widomInsertions;
		if (sim.volumeProb > 0){
			extPress = fluid.extPress/kb*1e-30; // K/AA^3
			muIdeal = fluid.temp*log(ops.Pow(thermalWL,3)*extPress/fluid.temp); //K
			muExcess = -fluid.temp*log(insertionParam*extPress/(fluid.temp*(fluid.nParts+1))); //K
		}else{
			muIdeal = fluid.temp*log(ops.Pow(thermalWL,3)*(fluid.nParts+1)/volume); //K
			muExcess = -fluid.temp*log(insertionParam); //K
		}
		fluid.mu = muIdeal+muExcess; //K
	}
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
void MC::PrintRDF(int set){
	double constant, deltaR, lowR, highR, rho, dn_ideal, gofr[NBINS+1];
	Operations ops;
	ofstream rdfFile;
	string outFileName="stats/rdf.dat";

	deltaR = fluid.rcut/NBINS;
	rho = fluid.nParts/ComputeVolume();
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
	cout << fixed;
	cout << setprecision(5);
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
		cout << "AcceptInsertRatio; RejectInsertRatio:\t";
		cout << stats.acceptInsertion*1./stats.nExchanges << "; ";
		cout << stats.rejectInsertion*1./stats.nExchanges << endl;
		cout << "AcceptDeletRAtio; RejectDeletRatio:\t";
		cout << stats.acceptDeletion*1./stats.nExchanges << "; ";
		cout << stats.rejectDeletion*1./stats.nExchanges << endl;
	}
	if (fluid.nParts > 0){
		cout << "NumPartcles; Energy/Particle:\t";
		cout << fluid.nParts << "; ";
		cout << sim.systemEnergy/fluid.nParts << endl;
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
	trajFile << "Lattice=\" " << pore.boxWidth[0] << " " << tmp << " " << tmp << " ";
	trajFile << tmp << " " << pore.boxWidth[1] << " " << tmp << " ";
	trajFile << tmp << " " << tmp << " " << pore.boxWidth[2] << "\" ";
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
	double density, volume;

	volume = ComputeVolume(); //AA^3
	density = fluid.nParts/volume; // AA^-3
	density *= fluid.molarMass/na*1e24; // g/cm^3
	if (!ops.FileExists(outFileName)){
		logFile.open(outFileName);
		logFile << "Set\tNParts\tDens[g/cm^3]\t";
		logFile << "boxWidth[AA]\tVolume[AA^3]\t";
		logFile << "ffE/part[K]\tsfE/part[K]\tE/part[K]\tmu[K]\n";
	}else{
		logFile.open(outFileName, ios::app);
		logFile << fixed;
		logFile << setprecision(5);
		logFile << set << "\t";
		logFile << fluid.nParts << "\t";
		logFile << density << "\t";
		logFile << pore.boxWidth[2] << "\t";
		logFile << volume << "\t";
		logFile << fluid.ffEnergy/fluid.nParts << "\t";
		logFile << pore.sfEnergy/fluid.nParts << "\t";
		logFile << sim.systemEnergy/fluid.nParts << "\t";
		logFile	<< fluid.mu << "\n";
	}
	logFile.close();
}
double MC::NeighDistance(int i, int j){
	double dist, dx, dy, dz;

	dx = part[j].x - part[i].x; //AA
	dy = part[j].y - part[i].y; //AA.
	dz = part[j].z - part[i].z; //AA.
	if (pore.PBC[0]) dx -= round(dx/pore.boxWidth[0]) * pore.boxWidth[0]; //AA
	if (pore.PBC[1]) dy -= round(dy/pore.boxWidth[1]) * pore.boxWidth[1]; //AA
	if (pore.PBC[2]) dz -= round(dz/pore.boxWidth[2]) * pore.boxWidth[2]; //AA
	dist = sqrt(dx*dx + dy*dy + dz*dz); //AA
	return dist; //AA
}
// ---------- Fluid-Fluid potentials ---------- //
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
// EAM Ga potential vvv //
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
// EAM Ga potential ^^^ //
// ---------- Fluid-Fluid potentials ---------- //

