//Author: Santiago A. Flores Roman

#include <cstring> // strlen
#include <iostream> // cout
#include <climits> // INT_MAX
#include <cmath> // isnan

#include "MC.h"
#include "tools.h"

using namespace std;

void MC::MinimizeEnergy(void){
	int initMoves = 200;

	for (int i=0; i<thermoSys.nBoxes; i++){
		if (box[i].nParts > 0){
			cout << "Minimizing initial configuration of box \"" << box[i].name << "\"..."<< endl;
			BoxEnergy(i);
			cout << "\tEnergy/part. of box \"" << box[i].name << "\" before minimization: ";
			cout << box[i].energy/box[i].nParts << " K" << endl;
			for (int j=0; j<initMoves; j++){
				for (int l=0; l<box[i].nParts; l++) MoveParticle();
				//AdjustMCMoves();
			}
			BoxEnergy(i);
			cout << "\tEnergy/part. of box \"" << box[i].name << "\" after minimization (";
			cout << box[i].nParts*initMoves << " moves): ";
			cout << box[i].energy/box[i].nParts << " K" << endl;
		}
	}
}
void MC::CorrectEnergy(void){
	double relEnergyChange;
	for (int i=0; i<thermoSys.nBoxes; i++){
		relEnergyChange = box[i].energy/box[i].oldEnergy;
		//if (relEnergyChange > 1 || isnan(relEnergyChange)) BoxEnergy(i);
		if (relEnergyChange > 1) BoxEnergy(i);
	}
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
	if (box[ithBox].fluid[ithSpecies].vdwPot[0] == "lj" && box[ithBox].geometry == "sphere"){
		box[ithBox].fluid[ithSpecies].particle[index].boxE += SphericalLJ_Pot(ithBox, ithSpecies, index);
	}else if (box[ithBox].fluid[ithSpecies].vdwPot[0] == "lj" && box[ithBox].geometry == "cylinder"){
		box[ithBox].fluid[ithSpecies].particle[index].boxE += CylindricalLJ_Pot(ithBox, ithSpecies, index);
	}else if (box[ithBox].fluid[ithSpecies].vdwPot[0] == "lj" && box[ithBox].geometry == "slit"){
		box[ithBox].fluid[ithSpecies].particle[index].boxE += SlitLJ_Pot(ithBox, ithSpecies, index);
	}
	//Fluid-Fluid energy
	for (int jthSpecies=0; jthSpecies<thermoSys.nSpecies; jthSpecies++){
		if (fluid[ithSpecies].vdwPot[jthSpecies] == "lj"){ //Lennard-Jones potential.
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += LJ126_Pot(ithBox, ithSpecies, jthSpecies, index);
		}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "hs"){
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += HardSphere_Pot(ithBox, ithSpecies, jthSpecies, index);
		}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "eam_na"){ //EAM potential for Na.
			tmp = EAMNa_Pot(ithBox, ithSpecies, jthSpecies, index);
			box[ithBox].fluid[ithSpecies].particle[index].manyBodyE += tmp[0];
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += tmp[1];
		}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "eam_k"){ //EAM potential for K.
			tmp = EAMK_Pot(ithBox, ithSpecies, jthSpecies, index);
			box[ithBox].fluid[ithSpecies].particle[index].manyBodyE += tmp[0];
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += tmp[1];
		}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "eam_rb"){ //EAM potential for Rb.
			tmp = EAMRb_Pot(ithBox, ithSpecies, jthSpecies, index);
			box[ithBox].fluid[ithSpecies].particle[index].manyBodyE += tmp[0];
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += tmp[1];
		}else if (fluid[ithSpecies].vdwPot[jthSpecies] == "eam_ga"){ //EAM potential for Ga.
			tmp = EAMGa_Pot(ithBox, ithSpecies, jthSpecies, index);
			box[ithBox].fluid[ithSpecies].particle[index].manyBodyE += tmp[0];
			box[ithBox].fluid[ithSpecies].particle[index].pairPotE += tmp[1];
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

