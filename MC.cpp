//Author: Santiago A. Flores202020202020 Roman

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cctype>
#include <iomanip>
#include <climits>
#include <filesystem>

#include "tools.h"
#include "MC.h"

using namespace std;

MC::MC(void){
	int i, j, k;
	ostringstream boxName;

	stats.acceptance = stats.rejection = stats.nDisplacements = 0;
	stats.acceptanceVol = stats.rejectionVol = stats.nVolChanges = 0;
	stats.acceptSwap = stats.rejectSwap = stats.nSwaps = 0;
	stats.widomInsertions = 0;
	for (i=0; i<MAXBOX; i++){
		for (j=0; j<MAXSPECIES; j++) stats.widom[i][j] = 0.;
	}
	for (i=0; i<=NBINS; i++) stats.rdf[i] = 0.;
	sim.rdf[0] = sim.rdf[1] = -1;
	sim.dr = 1;
	sim.dv = 1e-1;
	sim.nEquilSets = sim.nSets = sim.nCyclesPerSet = 0.;
	sim.nVolAttempts = sim.nSwapAttempts = 0.;
	thermoSys.temp = thermoSys.volume = thermoSys.press = -1.;
	thermoSys.nBoxes = 0;
	thermoSys.nSpecies = 0;
	thermoSys.nParts = 0;
	for (i=0; i<MAXBOX; i++){
		for (j=0; j<3; j++){
			box[i].width[j] = 0;
			box[i].PBC[j] = true;
		}
		for (j=0; j<MAXSPECIES; j++){
			box[i].fluid[j].nParts = 0;
			box[i].fluid[j].vdwPot[0] = "";
			box[i].fluid[j].epsilon[0] = box[i].fluid[j].sigma[0] = 0.;
			for (k=0; k<MAXPART; k++){
				box[i].fluid[j].particle[k].x = 0.;
				box[i].fluid[j].particle[k].y = 0.;
				box[i].fluid[j].particle[k].z = 0.;
				box[i].fluid[j].particle[k].energy = 0.;
			}
		}
		boxName << "Box" << i; box[i].name = boxName.str();
		box[i].geometry = "bulk";
		box[i].nLayersPerWall = 0;
		box[i].deltaLayers = 0.;
		box[i].solidDens = 0.;
		box[i].nParts = 0.;
		box[i].vdwPot = 0.;
		box[i].manyBodyE = 0.;
		box[i].pairPotE = 0.;
		box[i].boxE = 0.;
		box[i].energy = 0.;
	}
	for (i=0; i<MAXSPECIES; i++){
		fluid[i].name = "";
		fluid[i].molarMass = 0.;
		fluid[i].mu = 0.;
		for (j=0; j<MAXSPECIES; j++){
			fluid[i].vdwPot[j] = "";
			fluid[i].epsilon[j] = fluid[i].sigma[j] = fluid[i].rcut[j] = 0.;
		}
	}
}
double MC::ComputeVolume(Box ithBox){
	Tools tls;
	string geometry = ithBox.geometry;
	double xLength=ithBox.width[0], yLength=ithBox.width[1], zLength=ithBox.width[2];

	if (geometry == "sphere") return (pi/6.)*tls.Pow(zLength,3);
	else if (geometry == "cylinder") return (pi/4.)*tls.Pow(zLength,2)*xLength;
	else if (geometry == "slit") return xLength*yLength*zLength;
	else return tls.Pow(zLength,3); // Bulk phase.
}
void MC::ReadInputFile(string inFileName){
	int i, j, ithBox, ithSpecies, jthSpecies;
	double maxRcut;
	string* commands;
	string line;
	Tools tls;
	ifstream inFile;

	inFile.open(inFileName);
	if (!inFile.is_open()){ //Check if file was opened successfully.
		cout << "Error opening file" << endl;
		exit(EXIT_FAILURE);
	}
	while (getline(inFile, line)) {
		commands = tls.SplitString(line, ' ');
		for (i=0; i<4; i++) commands[i] = tls.LowerCase(commands[i]);
		// Simulation parameters.
		if (commands[0] == "productionsets") sim.nSets = stod(commands[1]);
		else if (commands[0] == "equilibriumsets") sim.nEquilSets = stod(commands[1]);
		else if (commands[0] == "cyclesperset") sim.nCyclesPerSet = stod(commands[1]);
		else if (commands[0] == "printeverynsets") sim.printEvery = stod(commands[1]);
		else if (commands[0] == "nvolumeattempts") sim.nVolAttempts = stoi(commands[1]);
		else if (commands[0] == "nswapattempts") sim.nSwapAttempts = stoi(commands[1]);
		else if (commands[0] == "externaltemperature") thermoSys.temp = stod(commands[1]); // K
		else if (commands[0] == "externalpressure") thermoSys.press = stod(commands[1]); // Pa
		// Species parameters.
		else if (commands[0] == "fluidname"){
			thermoSys.nSpecies++;
			fluid[thermoSys.nSpecies-1].name = commands[1];
		}else if (commands[0] == "molarmass") {fluid[thermoSys.nSpecies-1].molarMass = stod(commands[1]); // g/mol
		}else if (commands[0] == "chemicalpotential") {fluid[thermoSys.nSpecies-1].mu = stod(commands[1]);} // K
		// Box parameters.
		else if (commands[0] == "boxname"){
			thermoSys.nBoxes++;
			box[thermoSys.nBoxes-1].name = commands[1];
		}else if (commands[0] == "geometry") {box[thermoSys.nBoxes-1].geometry = commands[1];
		}else if (commands[0] == "surfacedensity") {box[thermoSys.nBoxes-1].solidDens = stod(commands[1]); // AA^-2
		}else if (commands[0] == "nlayersperwall") {box[thermoSys.nBoxes-1].nLayersPerWall = stoi(commands[1]);
		}else if (commands[0] == "distancebetweenlayers") {box[thermoSys.nBoxes-1].deltaLayers = stod(commands[1]);
		}else if (commands[0] == "size") {box[thermoSys.nBoxes-1].width[2] = stod(commands[1]);} // AA. width[2] (z axis) always keeps for pore size.
		// Interaction parameters.
		else if (commands[2] == "vdwpotential"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].vdwPot[jthSpecies] = commands[3];
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].vdwPot[ithSpecies] = commands[3];
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].vdwPot[0] = commands[3];
		}else if (commands[2] == "sigma"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].sigma[jthSpecies] = stod(commands[3]); // AA
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].sigma[ithSpecies] = stod(commands[3]); // AA
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].sigma[0] = stod(commands[3]); // AA
		}else if (commands[2] == "epsilon"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].epsilon[jthSpecies] = stod(commands[3]); // K
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].epsilon[ithSpecies] = stod(commands[3]); // K
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].epsilon[0] = stod(commands[3]); // AA
		}else if (commands[2] == "rcut"){
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[ithSpecies].rcut[jthSpecies] = stod(commands[3]); // AA
			if (ithSpecies >= 0 && jthSpecies >= 0) fluid[jthSpecies].rcut[ithSpecies] = stod(commands[3]); // AA
		}else if (commands[2] == "numberofparticles"){
			ithBox = tls.FindIndex(box, thermoSys.nBoxes, commands[0]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			if (ithBox >= 0 && jthSpecies >= 0) box[ithBox].fluid[jthSpecies].nParts = stoi(commands[3]);
		}
		// Sampling parameters.
		else if (commands[0] == "samplerdf"){
			ithSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[1]);
			jthSpecies = tls.FindIndex(fluid, thermoSys.nSpecies, commands[2]);
			sim.rdf[0] = ithSpecies;
			sim.rdf[1] = jthSpecies;
		}
		for (i=0; i<4; i++) commands[i] = "";
	}
	inFile.close();
	// Get remaining box widths and initial volumes according to the geometries given.
	maxRcut = 0.;
	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++) if (maxRcut < fluid[i].rcut[j]) maxRcut = fluid[i].rcut[j];
	}
	for (i=0; i<thermoSys.nBoxes; i++){
		if (box[i].geometry == "sphere") {box[i].width[0] = box[i].width[1] = box[i].width[2];
		}else if (box[i].geometry == "cylinder"){
			box[i].width[0] = 2*maxRcut; // PBC along x axis.
			box[i].width[1] = box[i].width[2];
		}else if (box[i].geometry == "slit") {box[i].width[0] = box[i].width[1] = 2*maxRcut; // PBC along xy plane.
		}else box[i].width[0] = box[i].width[1] = box[i].width[2]; // Bulk phase.
		box[i].volume = ComputeVolume(box[i]);
		thermoSys.volume += box[i].volume;
	}
	// Get total number of particles and density per box.
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j<thermoSys.nSpecies; j++) box[i].nParts += box[i].fluid[j].nParts;
		thermoSys.nParts += box[i].nParts;
	}
	// Set PBC restrictions.
	for (i=0; i<thermoSys.nBoxes; i++){
		if (box[i].geometry == "sphere"){
			box[i].PBC[0] = false;
			box[i].PBC[1] = false;
			box[i].PBC[2] = false;
		}else if (box[i].geometry == "cylinder"){
			box[i].PBC[0] = true;
			box[i].PBC[1] = false;
			box[i].PBC[2] = false;
		}else if (box[i].geometry == "slit"){
			box[i].PBC[0] = true;
			box[i].PBC[1] = true;
			box[i].PBC[2] = false;
		}else if (box[i].geometry == "bulk"){
			box[i].PBC[0] = true;
			box[i].PBC[1] = true;
			box[i].PBC[2] = true;
		}
	}
	// Set step size.
	//sim.dr = fluid[0].sigma[0]; //AA
}
void MC::PrintParams(void){
	int i, j, ithSpecies, jthSpecies;
	double maxRcut = 0.;

	cout << "Simulation parameters:" << endl;
	cout << "\tEquilibration sets: " << sim.nEquilSets << endl;
	cout << "\tProduction sets: " << sim.nSets << endl;
	cout << "\tCycles per set: " << sim.nCyclesPerSet << endl;
	cout << "\tPrint every n sets: " << sim.printEvery << endl;
	cout << "\tNum. of change volume attempts: " << sim.nVolAttempts << endl;
	cout << "\tNum. of swap attempts: " << sim.nSwapAttempts << endl;
	cout << endl;

	cout << "System parameters:" << endl;
	cout << "\tNum. of particles: " << thermoSys.nParts << endl;
	cout << "\tVolume (AA^3): " << thermoSys.volume << endl;
	cout << "\tExternal temperature (K): " << thermoSys.temp << endl;
	if (thermoSys.press >= 0) cout << "\tExternal pressure (Pa): " << thermoSys.press << endl;
	cout << endl;

	cout << "Fluid parameters:" << endl;
	for (i=0; i<thermoSys.nSpecies; i++){
		cout << "Species " << i << ": " << fluid[i].name << endl;
		cout << "\tMolar mass (g/mol): " << fluid[i].molarMass << endl;
		if (fluid[i].mu > 0) cout << "\tChemical potential (K): " << fluid[i].mu << endl;
	}
	cout << endl;

	cout << "Fluid-fluid interaction parameters: " << endl;
	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++){
			cout << fluid[i].name << "-" << fluid[j].name << endl;
			cout << "\tVdW potential: " << fluid[i].vdwPot[j] << endl;
			cout << "\tsigma (AA): " << fluid[i].sigma[j] << endl;
			cout << "\tepsilon (K): " << fluid[i].epsilon[j] << endl;
			cout << "\trcut (AA): " << fluid[i].rcut[j] << endl;
		}
	}
	cout << endl;

	cout << "Box parameters:" << endl;
	for (i=0; i<thermoSys.nBoxes; i++){
		cout << "Box " << i << ": " << box[i].name << endl;
		cout << "\tGeometry: " << box[i].geometry << endl;
		cout << "\tInitial width (AA): " << box[i].width[2] << endl;
		if (box[i].solidDens > 0) cout << "\tSurface density (AA^-2): " << box[i].solidDens << endl;
		if (box[i].nLayersPerWall > 1) cout << "\tLayers per wall: " << box[i].nLayersPerWall << endl;
		if (box[i].deltaLayers > 0) cout << "\tDistance between layers (AA): " << box[i].deltaLayers << endl;
		cout << "\tInitial number of particles: " << box[i].nParts << endl;
	}
	cout << endl;

	cout << "Box-fluid interaction parameters: " << endl;
	for (i=0; i<thermoSys.nBoxes; i++){
		for (j=0; j<thermoSys.nSpecies; j++){
			cout << box[i].name << "-" << fluid[j].name << endl;
			if (box[i].fluid[j].vdwPot[0] != "") cout << "\tVdW potential: " << box[i].fluid[j].vdwPot[0] << endl;
			if (box[i].fluid[j].sigma[0] > 0) cout << "\tsigma (AA): " << box[i].fluid[j].sigma[0] << endl;
			if (box[i].fluid[j].epsilon[0] > 0) cout << "\tepsilon (K): " << box[i].fluid[j].epsilon[0] << endl;
			cout << "\tInitial number of particles: " << box[i].fluid[j].nParts << endl;
		}
	}
	cout << endl;

	if (sim.rdf[0] > -1 && sim.rdf[1] > -1){
		ithSpecies = sim.rdf[0];
		jthSpecies = sim.rdf[1];
		cout << "Sample RDF: " << fluid[ithSpecies].name << "-" << fluid[jthSpecies].name << endl;
	}
	cout << endl;

	for (i=0; i<thermoSys.nSpecies; i++){
		for (j=i; j<thermoSys.nSpecies; j++){
			if (maxRcut < fluid[i].rcut[j]) maxRcut = fluid[i].rcut[j];
		}
	}
	for (i=0; i<thermoSys.nBoxes; i++){
		if ((box[i].width[2] < 2*maxRcut) && (box[i].geometry == "bulk")){
			cout << "Error: Box size must be larger than twice the cut-off radius." << endl;
			cout << "Box: " << box[i].name << endl;
			cout << endl;
			exit(EXIT_FAILURE);
		}
	}
}
void MC::OutputDirectory(void){
	Tools tls;
	ostringstream outDirName, boxDirName, command;

	outDirName << "output_" << thermoSys.nParts << "_" << thermoSys.temp << "K";
	if (! std::filesystem::is_directory(outDirName.str())){
		command << "mkdir " << outDirName.str();
		system(command.str().c_str());
		command.str(""); command.clear();
	}
	for (int i=0; i<thermoSys.nBoxes; i++){
		boxDirName << outDirName.str() << "/" << box[i].name;
		if (! filesystem::is_directory(boxDirName.str())){
			command << "mkdir " << boxDirName.str();
			system(command.str().c_str());
			command.str(""); command.clear();
		}
		boxDirName.str(""); boxDirName.clear();
	}
}
void MC::InsertParticle(int ithBox, int ithSpecies, int index){
	double radius, theta, phi; // Spherical geometry.
	double rho; // Cylindrical geometry.

	if (box[ithBox].geometry == "sphere"){
		radius = Random()*box[ithBox].width[2]*0.5;
		theta = Random()*pi;
		phi = Random()*2*pi;
		box[ithBox].fluid[ithSpecies].particle[index].x = radius*sin(theta)*cos(phi);
		box[ithBox].fluid[ithSpecies].particle[index].y = radius*sin(theta)*sin(phi);
		box[ithBox].fluid[ithSpecies].particle[index].z = radius*cos(theta);
		// Frame is at (0,0,0). Moving part. to box center.
		box[ithBox].fluid[ithSpecies].particle[index].x += 0.5*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y += 0.5*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z += 0.5*box[ithBox].width[2];
	}else if (box[ithBox].geometry == "cylinder"){
		rho = Random()*box[ithBox].width[2]*0.5;
		phi = Random()*2*pi;
		box[ithBox].fluid[ithSpecies].particle[index].x = (Random()-0.5)*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y = rho*sin(phi);
		box[ithBox].fluid[ithSpecies].particle[index].z = rho*cos(phi);
		// Frame is at (0,0,0). Moving part. to box center.
		box[ithBox].fluid[ithSpecies].particle[index].x += 0.5*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y += 0.5*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z += 0.5*box[ithBox].width[2];
	}else if (box[ithBox].geometry == "slit"){
		box[ithBox].fluid[ithSpecies].particle[index].x = Random()*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y = Random()*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z = Random()*box[ithBox].width[2];
	}else{
		box[ithBox].fluid[ithSpecies].particle[index].x = Random()*box[ithBox].width[0];
		box[ithBox].fluid[ithSpecies].particle[index].y = Random()*box[ithBox].width[1];
		box[ithBox].fluid[ithSpecies].particle[index].z = Random()*box[ithBox].width[2];
	}
}
void MC::InitialConfig(void){
	for (int i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++){
			for (int k=1; k<=box[i].fluid[j].nParts; k++) InsertParticle(i, j, k);
		}
	}
	cout << "Initial configuration set." << endl;
	cout << endl;
}
int* MC::GetMCMoves(void){
	static int moves[3];

	moves[0] = thermoSys.nParts;
	moves[1] = sim.nVolAttempts;
	moves[2] = sim.nSwapAttempts;
	return moves;
}
void MC::EnergyOfParticle(int ithBox, int ithSpecies, int index){
	Tools tls;
	double* tmp;
	double xPos, yPos, zPos, sqrPoreRadius;
	double sqrDistFromBoxCenter = 0.;

	box[ithBox].fluid[ithSpecies].particle[index].manyBodyE = 0.;
	box[ithBox].fluid[ithSpecies].particle[index].pairPotE = 0.;
	box[ithBox].fluid[ithSpecies].particle[index].boxE = 0.;

	sqrPoreRadius = box[ithBox].width[2]*box[ithBox].width[2]*0.25;
	xPos = box[ithBox].fluid[ithSpecies].particle[index].x;
	yPos = box[ithBox].fluid[ithSpecies].particle[index].y;
	zPos = box[ithBox].fluid[ithSpecies].particle[index].z;
	//Virtually, move positions to center of box.
	xPos -= 0.5*box[ithBox].width[0];
	yPos -= 0.5*box[ithBox].width[1];
	zPos -= 0.5*box[ithBox].width[2];
	if (box[ithBox].geometry == "sphere") sqrDistFromBoxCenter = xPos*xPos + yPos*yPos + zPos*zPos;
	else if (box[ithBox].geometry == "cylinder") sqrDistFromBoxCenter = yPos*yPos + zPos*zPos;
	else if (box[ithBox].geometry == "slit") sqrDistFromBoxCenter = zPos*zPos;
	// Particle must not be out of boundaries.
	if ((sqrDistFromBoxCenter > sqrPoreRadius) || (sqrDistFromBoxCenter < 0)){
		box[ithBox].fluid[ithSpecies].particle[index].boxE = INT_MAX;
	}
	//Particle-box energy
	if (box[ithBox].vdwPot == "lj" && box[ithBox].geometry == "sphere"){
		//box[ithBox].fluid[ithSpecies].particle[index].boxE += SphericalLJ(ithBox, ithPart);
	}
	else if (box[ithBox].vdwPot == "lj" && box[ithBox].geometry == "cylinder"){
		//box[ithBox].fluid[ithSpecies].particle[index].boxE += CylindricalLJ(ithBox, ithPart);
	}else if (box[ithBox].vdwPot == "lj" && box[ithBox].geometry == "slit"){
		//box[ithBox].fluid[ithSpecies].particle[index].boxE += SlitLJ(ithBox, ithPart);
	}
	//Fluid-Fluid energy
	for (int jthSpecies=0; jthSpecies<thermoSys.nSpecies; jthSpecies++){
		if (fluid[ithSpecies].vdwPot[jthSpecies] == "lj"){ //Lennard-Jones potential.
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += LJ_Pot(ithBox, ithSpecies, jthSpecies, index);
		}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "hs"){
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += HardSphere_Pot(ithBox, ithSpecies, jthSpecies, index);
		//}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "eam_ga"){ //EAM potential for Ga.
			//tmp = EAMGA_Energy(index);
			//box[ithBox].fluid[ithSpecies].particle[index].manyBodyE += tmp[0];
			//box[ithBox].fluid[ithSpecies].particle[index].pairPotE += tmp[1];
		}
	}
	box[ithBox].fluid[ithSpecies].particle[index].energy = box[ithBox].fluid[ithSpecies].particle[index].pairPotE;
	box[ithBox].fluid[ithSpecies].particle[index].energy += box[ithBox].fluid[ithSpecies].particle[index].manyBodyE;
	box[ithBox].fluid[ithSpecies].particle[index].energy += box[ithBox].fluid[ithSpecies].particle[index].boxE;
}
void MC::BoxEnergy(int ithBox){
	double manyBodyE, pairPotE, boxE;

	box[ithBox].manyBodyE = box[ithBox].pairPotE = box[ithBox].boxE = 0.;
	box[ithBox].energy = 0.;
	for (int i=0; i<thermoSys.nSpecies; i++){
		manyBodyE = pairPotE = boxE = 0.;
		for (int j=1; j<=box[ithBox].fluid[i].nParts; j++){
			EnergyOfParticle(ithBox, i, j);
			manyBodyE += box[ithBox].fluid[i].particle[j].manyBodyE;
			pairPotE += box[ithBox].fluid[i].particle[j].pairPotE;
			boxE += box[ithBox].fluid[i].particle[j].boxE;
		}
		box[ithBox].manyBodyE += manyBodyE;
		box[ithBox].pairPotE += 0.5*pairPotE;
		box[ithBox].boxE += boxE;
		box[ithBox].energy += manyBodyE + 0.5*pairPotE + boxE;
	}
}
// PBC for a box with the origin at the lower left vertex.
void MC::PBC(int ithBox, Particle& part){
	if (box[ithBox].PBC[0]) part.x -= floor(part.x/box[ithBox].width[0]) * box[ithBox].width[0]; //AA
	if (box[ithBox].PBC[1]) part.y -= floor(part.y/box[ithBox].width[1]) * box[ithBox].width[1]; //AA
	if (box[ithBox].PBC[2]) part.z -= floor(part.z/box[ithBox].width[2]) * box[ithBox].width[2]; //AA
}
void MC::MoveParticle(void){
	Particle ithPart;
	int i, ithBox, ithSpecies, index, tmp;
	double rand;
	double oldEnergy, newEnergy, deltaEnergy, arg;
	double oldPairPotE, oldManyBodyE, newPairPotE, newManyBodyE;
	double oldBoxEnergy, newBoxEnergy;
	double deltaManyBodyE=0, deltaPairPotE=0, deltaBoxEnergy=0;

	stats.nDisplacements++;
	// Select particle.
	rand = int(Random()*thermoSys.nParts);
	tmp = 0;
	for (i=0; i<thermoSys.nBoxes; i++){
		tmp += box[i].nParts;
		if (rand < tmp){
			ithBox = i;
			break;
		}
	}
	rand = int(Random()*box[ithBox].nParts);
	tmp = 0;
	for (i=0; i<thermoSys.nSpecies; i++){
		tmp += box[ithBox].fluid[i].nParts;
		if (rand < tmp){
			ithSpecies = i;
			break;
		}
	}
	index = int(Random()*box[ithBox].fluid[ithSpecies].nParts)+1;
	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	box[ithBox].fluid[ithSpecies].particle[0] = ithPart; //Save old position of selected particle.
	EnergyOfParticle(ithBox, ithSpecies, index);
	oldManyBodyE = box[ithBox].fluid[ithSpecies].particle[index].manyBodyE;
	oldPairPotE = box[ithBox].fluid[ithSpecies].particle[index].pairPotE;
	oldBoxEnergy = box[ithBox].fluid[ithSpecies].particle[index].boxE;
	oldEnergy = box[ithBox].fluid[ithSpecies].particle[index].energy;
	// Move particle.
	ithPart.x += (Random()-0.5)*sim.dr; //AA
	ithPart.y += (Random()-0.5)*sim.dr; //AA
	ithPart.z += (Random()-0.5)*sim.dr; //AA
	PBC(ithBox, ithPart);
	box[ithBox].fluid[ithSpecies].particle[index] = ithPart;
	// Metropolis.
	EnergyOfParticle(ithBox, ithSpecies, index);
	newManyBodyE = box[ithBox].fluid[ithSpecies].particle[index].manyBodyE;
	newPairPotE = box[ithBox].fluid[ithSpecies].particle[index].pairPotE;
	newBoxEnergy = box[ithBox].fluid[ithSpecies].particle[index].boxE;
	newEnergy = box[ithBox].fluid[ithSpecies].particle[index].energy;
	deltaEnergy = newEnergy-oldEnergy;
	arg = deltaEnergy/thermoSys.temp;
	if (Random() < exp(-arg)){
		stats.acceptance++;
		deltaManyBodyE = newManyBodyE - oldManyBodyE;
		deltaPairPotE = newPairPotE - oldPairPotE;
		deltaBoxEnergy = newBoxEnergy - oldBoxEnergy;
		box[ithBox].manyBodyE += deltaManyBodyE;
		box[ithBox].pairPotE += 0.5*deltaPairPotE;
		box[ithBox].boxE += deltaBoxEnergy;
		box[ithBox].energy += deltaManyBodyE + 0.5*deltaPairPotE + deltaBoxEnergy;
	}else{
		stats.rejection++;
		ithPart = box[ithBox].fluid[ithSpecies].particle[0];
		box[ithBox].fluid[ithSpecies].particle[index] = ithPart;
	}
}
void MC::MinimizeEnergy(void){
	int initMoves = 200;

	for (int i=0; i<thermoSys.nBoxes; i++){
		if (box[i].nParts > 0){
			cout << "Minimizing initial configuration of box \"" << box[i].name << "\"..."<< endl;
			BoxEnergy(i);
			cout << "\tEnergy/part. of box \"" << box[i].name << "\" before minimization: ";
			cout << box[i].energy/box[i].nParts << " K" << endl;
			for (int j=0; j<box[i].nParts; j++){
				for (int l=0; l<initMoves; l++) MoveParticle();
			}
			BoxEnergy(i);
			cout << "\tEnergy/part. of box \"" << box[i].name << "\" after minimization (";
			cout << box[i].nParts*initMoves << " moves): ";
			cout << box[i].energy/box[i].nParts << " K" << endl;
			cout << endl;
		}
	}
}
void MC::PrintStats(int set){
	cout << fixed;
	cout << setprecision(5);
	if (set < sim.nEquilSets) cout << "Equilibrium set: " << set << endl;
	else cout << "Set: " << set << endl;
	if (stats.nDisplacements > 0){
		cout << "AcceptDispRatio; RegectDispRatio:\t";
		cout << stats.acceptance*1./stats.nDisplacements << "; ";
		cout << stats.rejection*1./stats.nDisplacements << endl;
		if (set < sim.nEquilSets) cout << "Step size: " << sim.dr << endl;
	}
	if (sim.nVolAttempts > 0){
		cout << "AcceptVolRatio; RejectVolRatio:\t";
		cout << stats.acceptanceVol*1./stats.nVolChanges << "; ";
		cout << stats.rejectionVol*1./stats.nVolChanges << endl;
		if (set < sim.nEquilSets) cout << "Volume step size: " << sim.dv << endl;
	}
	//if (sim.exchangeProb > 0){
		//cout << "AcceptInsertRatio; RejectInsertRatio:\t";
		//cout << stats.acceptInsertion*1./stats.nExchanges << "; ";
		//cout << stats.rejectInsertion*1./stats.nExchanges << endl;
		//cout << "AcceptDeletRAtio; RejectDeletRatio:\t";
		//cout << stats.acceptDeletion*1./stats.nExchanges << "; ";
		//cout << stats.rejectDeletion*1./stats.nExchanges << endl;
	//}
	for (int i=0; i<thermoSys.nBoxes; i++){
		if (box[i].nParts > 0){
			cout << "Box " << box[i].name << ":" << endl;
			cout << "\tNumParticles; BoxSize; Energy/Particle:\t";
			cout << box[i].nParts << "; ";
			cout << box[i].width[2] << "; ";
			cout << box[i].energy/box[i].nParts << endl;
		}
	}
	cout << endl;
}
void MC::ResetStats(void){
	stats.acceptance = stats.rejection = 0;
	stats.nDisplacements = 1;
	stats.acceptanceVol = stats.rejectionVol = 0;
	stats.nVolChanges = 1;
	stats.acceptSwap = stats.rejectSwap = 0;
	stats.nSwaps = 1;
	stats.widomInsertions = 0;
	for (int i=0; i<thermoSys.nBoxes; i++){
		for (int j=0; j<thermoSys.nSpecies; j++) stats.widom[i][j] = 0.;
	}
}
void MC::AdjustMCMoves(void){
	double ratio;

	ratio = stats.acceptance*1./(1.*stats.nDisplacements);
	if (ratio < 0.2) sim.dr /= 1.1;
	else sim.dr *= 1.1;
	//if (sim.nVolAttempts > 0){
		//ratio = stats.acceptanceVol*1./(1.*stats.nVolChanges);
		//if (ratio < 0.3) sim.dv /= 1.1;
		//else sim.dv *= 1.1;
	//}
}
// Computes the length of the box along the z axis.
double MC::ComputeBoxWidth(Box ithBox, double volume){
	cout << ""; // To avoid Core dump.
	if (ithBox.geometry == "sphere") return 2*cbrt(3*volume/(4*pi));
	else if (ithBox.geometry == "cylinder") return 2*sqrt(volume/(pi*ithBox.width[0]));
	else if (ithBox.geometry == "slit") return volume/(ithBox.width[0]*ithBox.width[1]);
	else return cbrt(volume); //Bulk phase.
}
void MC::RescaleCenterOfMass(Box oldBox, Box& newBox){
	Tools tls;
	int i, j;
	double ratio = newBox.width[2]/oldBox.width[2];

	cout << ""; // To avoid core dump.
	for (i=0; i<thermoSys.nSpecies; i++){
		//if (pore.geometry == "sphere"){
			//for (i=1; i<=fluid.nParts; i++){
				//particle[i].x *= ratio;
				//particle[i].y *= ratio;
				//particle[i].z *= ratio;
			//}
		//}else if (pore.geometry == "cylinder"){
			//for (i=1; i<=fluid.nParts; i++){
				//particle[i].y *= ratio;
				//particle[i].z *= ratio;
			//}
		//}else if (pore.geometry == "slit"){
			//for (i=1; i<=fluid.nParts; i++) particle[i].z *= ratio;
		//}else{ //Bulk phase.
		for (j=1; j<=newBox.fluid[i].nParts; j++){
			newBox.fluid[i].particle[j].x *= ratio;
			newBox.fluid[i].particle[j].y *= ratio;
			newBox.fluid[i].particle[j].z *= ratio;
		}
	}
	//sim.dr *= ratio;
	//sim.dv *= ratio;
}
// Performs MC trial move: change of volume.
// Note: It applies either for Gibbs or NPT ensemble.
// Note 2: For NPT ensemble, it only works in box 1.
// Note 3: For Gibbs ensemble, it only works with 2 boxes.
void MC::ChangeVolume(void){
	double deltaE0, deltaE1;
	double logNewVol, deltaVol0, deltaVol1;
	double arg0, arg1;
	double extPress;
	Box oldBox0, oldBox1;

	stats.nVolChanges++;
	if (thermoSys.nBoxes == 1){
		extPress = thermoSys.press/kb*1e-30; // K/AA^3
		// Record current config.
		oldBox0 = box[0]; // Record old information.
		logNewVol = log(box[0].volume) + (Random()-0.5)*sim.dv; //Perform volume change step.
		box[0].volume = exp(logNewVol);
		RescaleCenterOfMass(oldBox0, box[0]); //Rescale to trial config.
		// Set box sizes according to trial volume.
		//if (pore.geometry == "sphere") pore.width[0] = pore.width[1] = pore.width[2] = finalBoxWidth;
		//else if (pore.geometry == "cylinder") pore.width[1] = pore.width[2] = finalBoxWidth;
		//else if (pore.geometry == "slit") pore.width[2] = finalBoxWidth;
		//else pore.width[0] = pore.width[1] = pore.width[2] = finalBoxWidth;
		box[0].width[2] = ComputeBoxWidth(box[0], box[0].volume);
		box[0].width[0] = box[0].width[1] = box[0].width[2];
		// Compute trial energy.
		BoxEnergy(0);
		deltaE0 = box[0].energy-oldBox0.energy;
		deltaVol0 = box[0].volume-oldBox0.volume;
		arg0 = -(deltaE0 + extPress*deltaVol0 - (box[0].nParts+1)*log(box[0].volume/oldBox0.volume)*thermoSys.temp);
		if (Random() > exp(arg0/thermoSys.temp)){ //Rejected trial move.
			stats.rejectionVol++;
			box[0] = oldBox0;
		}else{ // Accepted trial move. Keep trial config.
			stats.acceptanceVol++;
			thermoSys.volume = box[0].volume;
		}
	}else{
		// Record current configr
		oldBox0 = box[0]; // Record old information.
		oldBox1 = box[1];
		logNewVol = log(box[0].volume/box[1].volume) + (Random()-0.5)*sim.dv; //Perform volume change step.
		box[0].volume = thermoSys.volume*exp(logNewVol)/(1+exp(logNewVol));
		box[1].volume = thermoSys.volume-box[0].volume;
		RescaleCenterOfMass(oldBox0, box[0]); //Rescale to trial config.
		RescaleCenterOfMass(oldBox1, box[1]); //Rescale to trial config.
		// Set box sizes according to trial volume.
		box[0].width[2] = ComputeBoxWidth(box[0], box[0].volume);
		box[0].width[0] = box[0].width[1] = box[0].width[2];
		box[1].width[2] = ComputeBoxWidth(box[1], box[1].volume);
		box[1].width[0] = box[1].width[1] = box[1].width[2];
		// Compute trial energy.
		BoxEnergy(0);
		deltaE0 = box[0].energy-oldBox0.energy;
		deltaVol0 = box[0].volume-oldBox0.volume;
		arg0 = -(deltaE0 + extPress*deltaVol0 - (box[0].nParts+1)*log(box[0].volume/oldBox0.volume)*thermoSys.temp);
		BoxEnergy(1);
		deltaE1 = box[1].energy-oldBox1.energy;
		deltaVol1 = box[1].volume-oldBox1.volume;
		arg1 = -(deltaE1 + extPress*deltaVol1 - (box[1].nParts+1)*log(box[1].volume/oldBox1.volume)*thermoSys.temp);
		if (Random() > exp((arg0+arg1)/thermoSys.temp)){ //Rejected trial move.
			stats.rejectionVol++;
			box[0] = oldBox0;
			box[1] = oldBox1;
		}else stats.acceptanceVol++; // Accepted trial move. Keep trial config.
	}
}
// Performs MC trial move: change of volume.
// Note: It applies either for Gibbs or muVT ensemble.
// Note 2: For muVT ensemble, it only works in box 1.
// Note 3: For Gibbs ensemble, it only works with 2 boxes.
void MC::SwapParticle(void){
	Tools tls;
	int tmp, ithSpecies, inIndex, outIndex, inBox, outBox;
	double particleMass, thermalWL, zz;
	double inEnergy, outEnergy, deltaE, pairPotE, manyBodyE, boxE;
	double rand;
	double arg=0, temp=thermoSys.temp;

	if (thermoSys.nBoxes == 1){
		rand = int(Random()*box[0].nParts);
		tmp = 0;
		for (int i=0; i<thermoSys.nSpecies; i++){
			tmp += box[0].fluid[i].nParts;
			if (rand < tmp){
				ithSpecies = i;
				break;
			}
		}
		particleMass = fluid[ithSpecies].molarMass/na*1e-3; // kg/particle
		thermalWL = planck/sqrt(2.0*pi*particleMass*kb*thermoSys.temp)*1e10; // AA
		zz = exp(fluid[ithSpecies].mu/thermoSys.temp)/tls.Pow(thermalWL,3); // AA^-3
		if (Random() < 0.5) { // Try inserting particle.
			stats.nSwaps++;
			// Insert particle at random position.
			inIndex = box[0].fluid[ithSpecies].nParts+1;
			InsertParticle(0, ithSpecies, inIndex);
			EnergyOfParticle(0, ithSpecies, inIndex);
			inEnergy = box[0].fluid[ithSpecies].particle[inIndex].energy; // K
			arg = zz*box[0].volume*exp(-inEnergy/thermoSys.temp)/(box[0].nParts+1); // Acceptance criterion (for inserting particle).
			if (Random() < arg){
				stats.acceptSwap++; // Accepted: Insert particle.
				thermoSys.nParts++;
				box[0].nParts++;
				box[0].fluid[ithSpecies].nParts++;
				manyBodyE = box[0].fluid[ithSpecies].particle[inIndex].manyBodyE;
				pairPotE = 0.5*box[0].fluid[ithSpecies].particle[inIndex].pairPotE;
				boxE = box[0].fluid[ithSpecies].particle[inIndex].boxE;
				box[0].manyBodyE += manyBodyE; // K
				box[0].pairPotE += pairPotE; // K
				box[0].boxE += boxE; // K
				box[0].energy += manyBodyE + pairPotE + boxE;
			}else stats.rejectSwap++;
		}else if (box[0].nParts > 0){ // Try removing particle (only if there are particles in the box).
			stats.nSwaps++;
			// Select random particle.
			do{
				outIndex = int(Random()*box[0].fluid[ithSpecies].nParts)+1;
			}while(outIndex < box[0].fluid[ithSpecies].nParts);
			EnergyOfParticle(0, ithSpecies, outIndex); //K
			outEnergy = box[0].fluid[ithSpecies].particle[outIndex].energy; // K
			arg = box[0].nParts*exp(outEnergy/thermoSys.temp)/(zz*box[0].volume); // Acceptance criterion (for removing particle).
			if (Random() < arg){ // Accepted: Remove particle.
				stats.acceptSwap++;
				box[0].fluid[ithSpecies].particle[outIndex] = box[0].fluid[ithSpecies].particle[box[0].fluid[ithSpecies].nParts];
				thermoSys.nParts--;
				box[0].nParts--;
				box[0].fluid[ithSpecies].nParts--;
				manyBodyE = box[0].fluid[ithSpecies].particle[outIndex].manyBodyE;
				pairPotE = 0.5*box[0].fluid[ithSpecies].particle[outIndex].pairPotE;
				boxE = box[0].fluid[ithSpecies].particle[outIndex].boxE;
				box[0].manyBodyE -= manyBodyE; // K
				box[0].pairPotE -= pairPotE; // K
				box[0].boxE -= boxE; // K
				box[0].energy -= manyBodyE + pairPotE + boxE;
			}else stats.rejectSwap++;
		}
	}else{
		stats.nSwaps++;
		if (Random() < 0.5) {inBox = 0; outBox = 1;
		}else {inBox = 1; outBox = 0;}
		rand = int(Random()*box[inBox].nParts);
		tmp = 0;
		for (int i=0; i<thermoSys.nSpecies; i++){
			tmp += box[inBox].fluid[i].nParts;
			if (rand < tmp){
				ithSpecies = i;
				break;
			}
		}
		inIndex = box[inBox].fluid[ithSpecies].nParts+1;
		InsertParticle(inBox, ithSpecies, inIndex);
		EnergyOfParticle(inBox, ithSpecies, inIndex);
		inEnergy = box[inBox].fluid[ithSpecies].particle[inIndex].energy; // K
		stats.widom[inBox][ithSpecies] += box[inBox].volume*exp(-inEnergy/thermoSys.temp)/(box[inBox].nParts+1); // Compute mu through Widom.
		stats.widomInsertions++;
		if (box[outBox].nParts == 0) return;
		do{
			outIndex = int(Random()*box[outBox].fluid[ithSpecies].nParts)+1;
		}while(outIndex < box[outBox].fluid[ithSpecies].nParts);
		EnergyOfParticle(outBox, ithSpecies, outIndex); //K
		outEnergy = box[outBox].fluid[ithSpecies].particle[outIndex].energy; // K
		deltaE = inEnergy-outEnergy;
		arg = exp(-(deltaE+log(box[outBox].volume*(box[inBox].nParts+1)/(box[inBox].volume*box[outBox].nParts))*temp)/temp); // Acceptance criterion (for removing particle).
		if (Random() < arg){
			stats.acceptSwap++; // Accepted: Insert particle.
			box[inBox].nParts++;
			box[inBox].fluid[ithSpecies].nParts++;
			manyBodyE = box[inBox].fluid[ithSpecies].particle[inIndex].manyBodyE;
			pairPotE = 0.5*box[inBox].fluid[ithSpecies].particle[inIndex].pairPotE;
			boxE = box[inBox].fluid[ithSpecies].particle[inIndex].boxE;
			box[inBox].manyBodyE += manyBodyE; // K
			box[inBox].pairPotE += pairPotE; // K
			box[inBox].boxE += boxE; // K
			box[inBox].energy += manyBodyE + pairPotE + boxE;
			box[outBox].fluid[ithSpecies].particle[outIndex] = box[outBox].fluid[ithSpecies].particle[box[outBox].fluid[ithSpecies].nParts];
			box[outBox].nParts--;
			box[outBox].fluid[ithSpecies].nParts--;
			manyBodyE = box[outBox].fluid[ithSpecies].particle[outIndex].manyBodyE;
			pairPotE = 0.5*box[outBox].fluid[ithSpecies].particle[outIndex].pairPotE;
			boxE = box[outBox].fluid[ithSpecies].particle[outIndex].boxE;
			box[outBox].manyBodyE -= manyBodyE; // K
			box[outBox].pairPotE -= pairPotE; // K
			box[outBox].boxE -= boxE; // K
			box[outBox].energy -= manyBodyE + pairPotE + boxE;
		}else stats.rejectSwap++;
	}
}
// Performs Widom insertions to compute mu in each box.
// Note: It only works for simple fluids.
// Note 2: The Widom parameter is computed randomly in any box unless the chosen box is emtpy.
void MC::ComputeWidom(void){
	Tools tls;
	double volume;
	double energy=0.;

	if (thermoSys.nBoxes == 1 && sim.nSwapAttempts == 0){
		InsertParticle(0, 0, MAXPART-1); // Insert virtual particle.
		stats.widomInsertions++;
		EnergyOfParticle(0, 0, MAXPART-1);
		energy = box[0].fluid[0].particle[MAXPART-1].energy;
		if (sim.nVolAttempts > 0){
			stats.widom[0][0] += box[0].volume*exp(-energy/thermoSys.temp)/(box[0].fluid[0].nParts+1);
		}else stats.widom[0][0] += exp(-energy/thermoSys.temp);
		stats.widom[0][0] += exp(-energy/thermoSys.temp);
	}
}
// Computes mu in each box from the Widom insearsions performed in the ComputeWidom method.
// Note: It only works for simple fluids.
// Note 2: the chemical potential is printed in all boxes.
void MC::ComputeChemicalPotential(void){
	Tools tls;
	double thermalWL, particleMass, volume, muIdeal, muExcess, insertionParam;
	double extPress;

	if (thermoSys.nBoxes == 1 && sim.nSwapAttempts == 0){
		particleMass = fluid[0].molarMass/na*1e-3; // kg/particle
		thermalWL = planck/sqrt(2.0*pi*particleMass*kb*thermoSys.temp)*1e10; //AA
		volume = box[0].volume; //AA^3
		insertionParam = stats.widom[0][0]/stats.widomInsertions;
		if (sim.nVolAttempts > 0){
			extPress = thermoSys.press/kb*1e-30; // K/AA^3
			muIdeal = thermoSys.temp*log(tls.Pow(thermalWL,3)*extPress/thermoSys.temp); //K
			muExcess = -thermoSys.temp*log(insertionParam*extPress/thermoSys.temp); //K
		}else{
			muIdeal = thermoSys.temp*log(tls.Pow(thermalWL,3)*(box[0].fluid[0].nParts+1)/volume); //K
			muExcess = -thermoSys.temp*log(insertionParam); //K
		}
	}else{
		muIdeal = thermoSys.temp*log(tls.Pow(thermalWL,3)); //K
		for (int i=0; i<thermoSys.nSpecies; i++){
			insertionParam = stats.widom[0][i]/stats.widomInsertions;
			muExcess = -thermoSys.temp*log(insertionParam); //K
			box[0].fluid[i].muEx = muExcess; //K
			box[0].fluid[i].mu = muIdeal+muExcess; //K
			insertionParam = stats.widom[1][i]/stats.widomInsertions;
			muExcess = -thermoSys.temp*log(insertionParam); //K
			box[1].fluid[i].muEx = muExcess; //K
			box[1].fluid[i].mu = muIdeal+muExcess; //K
		}
	}
}
// Computes the radial distribution function every MC cycle.
// Note: Due to the expensive PC cost, it only computes the RDF in box 0 for the chosen pair of species.
void MC::ComputeRDF(void){
	int bin, ithSpecies, jthSpecies;
	double deltaR;
	Tools tls;
	Particle ithPart, jthPart;

	if (sim.rdf[0] > -1 && sim.rdf[1] > -1){
		ithSpecies = sim.rdf[0];
		jthSpecies = sim.rdf[1];
		deltaR = fluid[ithSpecies].rcut[jthSpecies]/(1.*NBINS);
		if (ithSpecies == jthSpecies){
			for (int i=1; i<=int(box[0].fluid[ithSpecies].nParts)-1; i++){
				for (int j=i+1; j<=int(box[0].fluid[ithSpecies].nParts); j++){
					ithPart = box[0].fluid[ithSpecies].particle[i];
					jthPart = box[0].fluid[ithSpecies].particle[j];
					bin = int(NeighDistance(0, ithPart, jthPart)/deltaR)+1;
					if (bin <= NBINS) stats.rdf[bin] += 2; // Takes into account both i->j and j->i.
				}
			}
		}else{
			for (int i=1; i<=int(box[0].fluid[ithSpecies].nParts); i++){
				for (int j=1; j<=int(box[0].fluid[jthSpecies].nParts); j++){
					ithPart = box[0].fluid[ithSpecies].particle[i];
					jthPart = box[0].fluid[jthSpecies].particle[j];
					bin = int(NeighDistance(0, ithPart, jthPart)/deltaR)+1;
					if (bin <= NBINS) stats.rdf[bin] += 2; // Takes into account both i->j and j->i.
				}
			}
		}
	}
}
// Prints the radial distribution function into a file.
// Note: RDF is printed after simulation finishes.
void MC::PrintRDF(void){
	double constant, deltaR, lowR, highR, rho, dn_ideal, gofr[NBINS+1];
	int ithSpecies, jthSpecies, nParts;
	Tools tls;
	ofstream rdfFile;
	ostringstream outDirName, outFileName;

	if (sim.rdf[0] > -1 && sim.rdf[1] > -1){
		ithSpecies = sim.rdf[0];
		jthSpecies = sim.rdf[1];
		deltaR = fluid[ithSpecies].rcut[jthSpecies]/NBINS;
		if (ithSpecies == jthSpecies) nParts = box[0].fluid[ithSpecies].nParts;
		else nParts = box[0].fluid[ithSpecies].nParts + box[0].fluid[jthSpecies].nParts;
		rho = nParts/box[0].volume;
		constant = (4*pi*rho/3);
		for (int bin=1; bin<=NBINS; bin++){
			gofr[bin] = stats.rdf[bin]/(1.*nParts*sim.nCyclesPerSet*sim.nSets);
			lowR = bin*deltaR;
			highR = lowR + deltaR;
			dn_ideal = constant*(tls.Pow(highR,3)-tls.Pow(lowR,3));
			gofr[bin] /= dn_ideal;
		}
		//Write into the file.
		outDirName << "output_" << thermoSys.nParts << "_" << thermoSys.temp << "K";
		outFileName << outDirName.str() << box[0].name << "/rdf_";
		outFileName << fluid[ithSpecies].name << "-" << fluid[jthSpecies].name << ".dat";
		rdfFile.open(outFileName.str());
		rdfFile << "nBins: "<< NBINS << endl;
		rdfFile << "r[AA]\tg(r)" << endl;
		for (int bin=1; bin<=NBINS; bin++) rdfFile << (bin+0.5)*deltaR << "\t" << gofr[bin] << endl;
		rdfFile.close();
	}
}
void MC::PrintTrajectory(int set){
	Tools tls;
	ofstream trajFile;
	Particle part;
	ostringstream outDirName, outFileName;
	double tmp=0.0;

	outDirName << "output_" << thermoSys.nParts << "_" << thermoSys.temp << "K";
	for (int i=0; i<thermoSys.nBoxes; i++){
		outFileName << outDirName.str() << "/" << box[i].name << "/trajectory.exyz";
		if (set == 0) trajFile.open(outFileName.str());
		else trajFile.open(outFileName.str(), ios::app);
		trajFile << box[i].nParts << "\n";
		trajFile << "Lattice=\" " << box[i].width[0] << " " << tmp << " " << tmp << " ";
		trajFile << tmp << " " << box[i].width[1] << " " << tmp << " ";
		trajFile << tmp << " " << tmp << " " << box[i].width[2] << "\" ";
		trajFile << "Properties=species:S:1:pos:R:3 Time=" << 1.*set << "\n";
		for (int j=0; j<thermoSys.nSpecies; j++){
			for (int k=1; k<=box[i].fluid[j].nParts; k++){
				part = box[i].fluid[j].particle[k];
				trajFile << fluid[j].name << "\t";
				trajFile << part.x << "\t" << part.y << "\t" << part.z << "\n"; //AA
			}
		}
		trajFile.close();
		outFileName.str(""); outFileName.clear();
		trajFile.clear();
	}
}
void MC::PrintLog(int set){
	Tools tls;
	ofstream logFile;
	ostringstream outDirName, outBoxName, outFileName;
	double density, volume, ffEnergy, energy;

	outDirName << "output_" << thermoSys.nParts << "_" << thermoSys.temp << "K";
	for (int i=0; i<thermoSys.nBoxes; i++){
		outFileName << outDirName.str() << "/" << box[i].name << "/simulation.log";
		volume = box[i].volume; //AA^3
		ffEnergy = box[i].manyBodyE + 0.5*box[i].pairPotE;
		energy = box[i].manyBodyE + 0.5*box[i].pairPotE + box[i].boxE;
		if (set == 0){
			logFile.open(outFileName.str());
			logFile << "Set\tTemp[K]\twidth[AA]\tVolume[AA^3]\t";
			logFile << "ffE/particle[K]\tsfE/particle[K]\tE/particle[K]\t";
			for (int j=0; j<thermoSys.nSpecies; j++){
				logFile << "\tNParts_" << fluid[j].name << "\t";
				logFile << "Dens[g/cm^3]_" << fluid[j].name << "\t";
				logFile << "muEx[K]_" << fluid[j].name << "\t";
				logFile << "mu[K]_" << fluid[j].name << "\n";
			}
		}else{
			logFile.open(outFileName.str(), ios::app);
			logFile << fixed;
			logFile << setprecision(5);
			logFile << set << "\t";
			logFile << thermoSys.temp << "\t" << box[i].width[2] << "\t" << volume << "\t";
			logFile << ffEnergy/box[i].nParts << "\t";
			logFile << box[i].boxE/box[i].nParts << "\t";
			logFile << energy/box[i].nParts << "\t";
			for (int j=0; j<thermoSys.nSpecies; j++){
				density = box[i].fluid[j].nParts/volume; // AA^-3
				density *= fluid[j].molarMass/na*1e24; // g/cm^3
				logFile << box[i].nParts << "\t";
				logFile << density << "\t"; // g/cm^3
				logFile	<< box[i].fluid[j].muEx << "\t"; // K
				logFile	<< box[i].fluid[j].mu << "\n"; // K
			}
		}
		logFile.close();
		outFileName.str(""); outFileName.clear();
		logFile.clear();
	}
}

double MC::NeighDistance(int ithBox, Particle ithPart, Particle jthPart){
	double dist, dx, dy, dz;

	dx = ithPart.x - jthPart.x; //AA
	dy = ithPart.y - jthPart.y; //AA
	dz = ithPart.z - jthPart.z; //AA
	if (box[ithBox].PBC[0]) dx -= round(dx/box[ithBox].width[0]) * box[ithBox].width[0]; //AA
	if (box[ithBox].PBC[1]) dy -= round(dy/box[ithBox].width[1]) * box[ithBox].width[1]; //AA
	if (box[ithBox].PBC[2]) dz -= round(dz/box[ithBox].width[2]) * box[ithBox].width[2]; //AA
	dist = sqrt(dx*dx + dy*dy + dz*dz); //AA
	return dist; //AA
}
// ---------- Fluid-Fluid potentials ---------- //
double MC::HardSphere_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Tools tls;
	Particle ithPart, jthPart;
	double rij, sig;
	double energy=0;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	sig = fluid[ithSpecies].sigma[jthSpecies]; // AA
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		rij = NeighDistance(ithBox, ithPart, jthPart);
		if (rij <= sig) energy += INT_MAX;
	}
	return energy; //K
}
double MC::LJ_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Tools tls;
	Particle ithPart, jthPart;
	double rij, eps, sig, rcut;
	double energy=0;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	eps = fluid[ithSpecies].epsilon[jthSpecies]; // K
	sig = fluid[ithSpecies].sigma[jthSpecies]; // AA
	rcut = fluid[ithSpecies].rcut[jthSpecies]; // AA
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= rcut) energy += 4.*eps*(tls.Pow(sig/rij,12)-tls.Pow(sig/rij,6));
		}
	}
	return energy; //K
}
// EAM Ga potential vvv //
// Paper:
// Belashchenko, D.K., 2012.
// Computer Simulation of the Properties of Liquid Metals: Gallium, Lead, and Bismuth.
// Russ. J. Phys. Chem. A, 86, pp.779-790.
//double* MC::EAMGA_Energy(int index){
	//static double energy[2];
	//double rij;
	//double rho=0, phiLow=0;
	//double eVToK = (801088317./5.0e27)/kb;

	//for (int j=1; j<=fluid.nParts; j++){
		//if (j != index){
			//rij = NeighDistance(index, j);
			//if (rij < fluid.rcut){
				//rho += eDens(rij);
				//phiLow += PairPot(rij);
			//}
		//}
	//}
	//energy[0] = EmbPot(rho)*eVToK;
	//energy[1] = phiLow*eVToK;
	//return energy;
//}
//double MC::StepUnit(double radius, double leftLim, double rightLim){
	//if (leftLim <= radius && radius < rightLim) return 1.;
	//else return 0.;
//}
//double MC::EmbPot(double rho){
	//double rhoIntervals[7] = {1.00000,  0.92000,  0.870000,  0.800000,  0.750000,  0.650000,  1.400000};
	//double aValues[7]      = {0.00000, -1.91235, -1.904030, -1.897380, -1.883520, -1.852620, -1.822820};
	//double bValues[7]      = {0.00000,  0.00000, -0.208000, -0.058000, -0.338000, -0.898000,  0.302000};
	//double cValues[7]      = {0.00000,  1.30000, -1.500000,  2.000000,  5.600000, -6.000000,  2.000000};
	//double phi = 0;

	//if (rhoIntervals[1] <= rho && rho <= rhoIntervals[6]) phi = aValues[1] + cValues[1]*(rho-rhoIntervals[0])*(rho-rhoIntervals[0]);
	//for (int i=2; i<6; i++){
		//if (rhoIntervals[i] <= rho && rho <= rhoIntervals[i-1]){
			//phi = aValues[i] + bValues[i]*(rho-rhoIntervals[i-1]) + cValues[i]*(rho-rhoIntervals[i-1])*(rho-rhoIntervals[i-1]);
			//break;
		//}
	//}
	//if (rho <= rhoIntervals[5]){
		//phi = (aValues[6] + bValues[6]*(rho-rhoIntervals[5]) + cValues[6]*(rho-rhoIntervals[5])*(rho-rhoIntervals[5])) * (2*rho/rhoIntervals[5]-(rho/rhoIntervals[5])*(rho/rhoIntervals[5]));
	//}
	//return phi;
//}
//double MC::eDens(double radius){
	//double pValues[3] = {0, 2.24450, 1.2};

	//return pValues[1]*exp(-pValues[2]*radius);
//}
//double MC::PairPot(double radius){
	//double rIntervals[7] = {0.0, 2.15, 2.75, 3.35, 4.00, 6.50, 8.30};
	//double aValues[9][6] = {{0.0, -0.65052509307861e-01, -0.15576396882534e+00, -0.13794735074043e+00, 0.13303710147738e-01,  0.00000000000000e+00},
							//{0.0, -0.32728102803230e+00, -0.16365580260754e-01,  0.78778542578220e-01, 0.59769893996418e-02,  0.00000000000000e+00},
							//{0.0,  0.51590444127493e+01,  0.20955204046244e+00, -0.83622260891495e-01, 0.57411338894840e-01, -0.60454444423660e-02},
							//{0.0,  0.90195221829217e+02, -0.97550604734748e+00, -0.44410858010987e+01, 0.19517888219051e+00, -0.13258585494287e+00},
							//{0.0,  0.72322004859499e+03, -0.11625479189815e+02, -0.36415106938231e+02, 0.32162310059276e+00, -0.34988482891053e+00},
							//{0.0,  0.27788989409594e+04, -0.58549935696765e+02, -0.13414583419234e+03, 0.30195698240893e+00, -0.45183606796559e+00},
							//{0.0,  0.56037895713613e+04, -0.15186293377510e+03, -0.25239146992011e+03, 0.14850603977640e+00, -0.31733856650298e+00},
							//{0.0,  0.57428084950480e+04, -0.19622924502226e+03, -0.23858760191913e+03, 0.36233874262589e-01, -0.11493645479281e+00},
							//{0.0,  0.23685488320885e+04, -0.98789413798382e+02, -0.90270667293646e+02, 0.34984220138018e-02, -0.16768950999376e-01}};
	//double phi=0;

	//if (rIntervals[1] < radius && radius <= rIntervals[6]){
		//for (int i=1; i<6; i++){
			//for (int m=0; m<9; m++){
				//phi += aValues[m][i] * pow((radius-rIntervals[i+1]),m) * StepUnit(radius, rIntervals[i], rIntervals[i+1]);
			//}
		//}
	//}else if (rIntervals[0] < radius && radius <= rIntervals[1]){
		//phi = 0.619588 - 51.86268*(2.15-radius) + 27.8*(exp(1.96*(2.15-radius))-1);
	//}
	//return phi;
//}
// EAM Ga potential ^^^ //
// ---------- Fluid-Fluid potentials ---------- //

// ---------- Solid-Fluid potentials ---------- //
// Paper:
// Steele, W.A., 1973.
// The Physical Interaction of Gases with Crystalline Solids: I. Gas-Solid Energies and Properties of Isolated Adsorbed Atoms.
// Surf. Sci., 36(1), pp.317-352.
// Steele 10-4-3 potential extended by Jason for multiple layers.
// Domain of potential: z in [0,H], where H is the pore size.
//double MC::SlitLJ(int index){
	//Tools tls;
	//double z, t1, t2, t3, t4;
	//int nLayers = pore.nLayersPerWall;
	//double sizeR = pore.width[2]*0.5; // AA
	//double dens = pore.sfDensity; // AA^-2
	//double eps = pore.sfEpsilon; // K
	//double sig = pore.sfSigma; // AA
	//double delta = pore.deltaLayers; // AA
	//double usf = 0;

	//z = particle[index].z-sizeR;
	//for (int i=0; i<nLayers; i++){
		//t1 = tls.Pow(sig/(sizeR+i*delta+z),10);
		//t2 = tls.Pow(sig/(sizeR+i*delta-z),10);
		//t3 = tls.Pow(sig/(sizeR+i*delta+z),4);
		//t4 = tls.Pow(sig/(sizeR+i*delta-z),4);
		//usf += 0.2*(t1+t2)-0.5*(t3+t4);
	//}
	//usf *= 4.*pi*eps*dens*sig*sig; // K
	//return usf;
//}
// Paper:
// Tjatjopoulos et al., 1988.
// Molecule-Micropore Interaction Potentials.
// J. Phys. Chem., 92(13), pp.4006-4007.
// Source for hypergeometric function:
// Press, W.H., 2007.
// Numerical recipes 3rd edition: The art of scientific computing.
// Cambridge university press.
//double MC::CylindricalLJ(int index){
	//Tools tls;
	//double yPos, zPos;
	//double x, u1, u2, usf, rho;
	//double sizeR = pore.width[2]*0.5; // AA
	//double dens = pore.sfDensity; // AA^-2
	//double eps = pore.sfEpsilon; // K
	//double sig = pore.sfSigma; // AA

	//yPos = particle[index].y - 0.5*pore.width[1];
	//zPos = particle[index].z - 0.5*pore.width[2];
	//rho = sqrt(tls.Pow(yPos,2) + tls.Pow(zPos,2));
	//// Reducing units/
	//dens *= sig*sig;
	//sizeR /= sig;
	//rho /= sig;
	//x = sizeR-rho;
	//u1 = 63. /32. * 1./(tls.Pow(x,10)*tls.Pow(2.0-x/sizeR,10)) * hypgeo(-4.5,-4.5,1.,tls.Pow(1-x/sizeR,2));
	//u2 = 3. / (tls.Pow(x,4)*tls.Pow(2.0-x/sizeR,4)) * hypgeo(-1.5,-1.5,1.,tls.Pow(1.0-x/sizeR,2));
	//usf = (u1-u2)*pi*pi*dens*eps; // K
	//return usf;
//}
// Paper:
// Baksh, M.S.A. and Yang, R.T., 1991.
// Model for Spherical Cavity Radii and Potential Functions of Sorbates in Zeolites.
// AIChE J., 37(6), pp.923-930.
//double MC::SphericalLJ(int index){
	//Tools tls;
	//double xPos, yPos, zPos;
	//double r, x, usf;
	//double sizeR = pore.width[2]*0.5;
	//double dens = pore.sfDensity;
	//double eps = pore.sfEpsilon;
	//double sig = pore.sfSigma;
	//double u1=0, u2=0;

	//xPos = particle[index].x - 0.5*pore.width[0];
	//yPos = particle[index].y - 0.5*pore.width[1];
	//zPos = particle[index].z - 0.5*pore.width[2];
	//r = sqrt(tls.Pow(xPos,2) + tls.Pow(yPos,2) + tls.Pow(zPos,2));
	//// Reducing units/
	//dens *= sig*sig;
	//sizeR /= sig;
	//r /= sig;
	//x = sizeR-r;
	//for (int i=0; i<10; i++){
		//u1 += 1. / tls.Pow(sizeR,i) / tls.Pow(x,10-i);
		//u1 += tls.Pow(-1.0,i) / tls.Pow(sizeR,i) / tls.Pow(x-2.0*sizeR,10-i);
		//if (i < 4){
			//u2 += 1. / tls.Pow(sizeR,i) / tls.Pow(x,4-i);
			//u2 += tls.Pow(-1.0,i) / tls.Pow(sizeR,i) / tls.Pow(x-2.0*sizeR,4-i);
		//}
	//}
	//usf = (2./5.) * u1 - u2;
	//usf *= 2.*pi*eps*dens; //K
	//return usf;
//}
// ---------- Solid-Fluid potentials ---------- //

