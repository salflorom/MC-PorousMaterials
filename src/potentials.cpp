//Author: Santiago A. Flores Roman
// Note:
// This file contains all fluid-fluid and solid-fluid potentials.
// To add a new potential:
// 1. Define the potential in this file.
// 2. Declare the potential in MC.h.
// 3. Add it in energy.cpp -> EnergyOfParticle(int, int, int),
//    in its respective section: Particle-box enregy or Fluid-Fluid energy.
//

#include <cstring> // strlen
#include <cmath> // pi, exp, round
#include <climits> // INT_MAX

#include "NumRecipes/hypgeo.h" // hypgeo
#include "MC.h"
#include "tools.h"

using namespace std;


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
double MC::StepUnit(double xValue, double leftLim, double rightLim){
	if (leftLim <= xValue && xValue < rightLim) return 1.;
	else return 0.;
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
double MC::LJ126_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
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
// EAM Na potential //
// Important: This potential is still under test. It is not reliable yet!
// Paper 1:
// Belashchenko, D.K., 2012.
// Electron Contribution to Energy of Alkali Metals in the Scheme of an Embedded Atom Model.
// High Temp., 50(3), pp.331-339.
// Paper 2:
// Belashchenko, D.K., 2009. 
// Application of the embedded atom model to liquid metals: Liquid sodium. 
// High Temp., 47, pp.494-507.
// Cut-off radius: 10.78 Angstrom.
double* MC::EAMNa_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Particle ithPart, jthPart;
	static double energy[2];
	double rij;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= 10.78){
				rho += EAMNa_eDens(rij);
				phiLow += EAMNa_PairPot(rij);
			}
		}
	}
	energy[0] = EAMNa_EmbPot(rho)*eVToK;
	energy[1] = phiLow*eVToK;
	return energy;
}
double MC::EAMNa_EmbPot(double rho){
	Tools tls;
	double rhoIntervals[10] = {1.00000,  0.90000,  0.800000,  0.700000,  0.620000,  0.280000,  0.050000,  0.030000,  1.400000,  0.000000};
	double aValues[10]      = {0.00000, -0.33140, -0.329781, -0.328243, -0.331155, -0.331845, -0.216727,  0.194462,  0.253218, -0.305496};
	double bValues[10]      = {0.00000,  0.00000, -0.032380,  0.001620,  0.056620, -0.039380, -0.637780, -2.937780, -2.937780,  0.129520};
	double cValues[10]      = {0.00000,  0.16190, -0.170000, -0.275000,  0.600000,  0.880000,  5.000000,  0.000000,  0.000000, -0.150000};
	double mValue=1.15, phi = 0;

	if (rhoIntervals[1] <= rho && rho < rhoIntervals[8]) phi = aValues[1] + cValues[1]*tls.Pow(rho-rhoIntervals[0],2);
	for (int i=2; i<8; i++){
		if (rhoIntervals[i] <= rho && rho < rhoIntervals[i-1]){
			phi = aValues[i] + bValues[i]*(rho-rhoIntervals[i-1]) + cValues[i]*tls.Pow(rho-rhoIntervals[i-1],2);
			break;
		}
	}
	if (rho < rhoIntervals[7]){
		phi = aValues[8] + bValues[8]*(rho-rhoIntervals[7])+cValues[8]*tls.Pow(rho-rhoIntervals[7],2);
		phi *= 2.*rho/rhoIntervals[5]-(rho/tls.Pow(rhoIntervals[5],2));
	}
	if (rho >= rhoIntervals[8]) phi = aValues[9] + bValues[9]*(rho-rhoIntervals[8]) + cValues[9]*pow(rho-rhoIntervals[8],mValue);
	return phi;
}
double MC::EAMNa_eDens(double radius){
	double pValues[3] = {0, 3.4418, 1.0245};
	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::EAMNa_PairPot(double radius){
	Tools tls;
	double phi;
	double rIntervals[12]= {0.0, 2.55, 2.80, 2.95, 3.45, 3.95, 4.45, 4.95, 5.45, 5.95, 7.45, 10.78};
	double bValues[12][7]= {{ 0.00000000000,  0.0000000000,  00.00000000000,  000.00000000000,  0000.00000000000,  0000.00000000000,  00.00000000000},
							{ 0.00000000000,  0.0000000000,  00.00000000000,  000.00000000000,  0000.00000000000,  0000.00000000000,  00.00000000000},
							{ 0.35805506000, -2.8231320000,  12.57403700000,  324.38852000000,  1675.13890000000,  2599.46790000000,  00.00000000000},
							{ 0.12708218000, -0.7885621800,  01.46133970000, -023.64693200000,  0000.00000000000,  0000.00000000000,  00.00000000000},
							{-0.11093583000, -0.3013956200, -00.57444694000, -007.67105920000, -0029.51319800000, -0053.35320300000, -35.24244200000},
							{-0.18380286000, -0.0312040760,  00.66662912000,  004.86362490000,  0020.31048900000,  0037.53878100000,  25.74349800000},
							{-0.17446597000,  0.0949913110,  00.61158912000,  003.82919210000,  0012.55137100000,  0018.51867900000,  10.12214300000},
							{-0.13020295000,  0.0879272820, -00.56735449000, -005.27961900000, -0019.75187800000, -0033.86296700000, -21.89686400000},
							{-0.07368651000,  0.1322006400,  00.75321670000,  006.82070260000,  0026.09383900000,  0045.49814000000,  29.74844400000},
							{-0.02636867600,  0.0793784780, -00.07650342800, -000.64246319000, -0003.66128030000, -0008.76413680000, -07.28703330000},
							{ 0.02854006900, -0.0019879458,  00.00569042080,  000.08022363800,  0000.09918809300,  0000.06169853000,  00.01461238100},
							{ 0.73629497e-4,  0.0000000000, -00.40257317e-2, -000.52993510e-2, -0000.16744300e-2, -0000.42268470e-3, -00.65802071e-4}};

	if (rIntervals[0] <= radius && radius <= rIntervals[1]) phi = 0.786149*exp(1.2*(2.55-radius));
	else{
		for (int m=2; m<12; m++){
			if (rIntervals[m-1] < radius && radius <= rIntervals[m]){
				phi = 0;
				for (int n=0; n<7; n++) phi += bValues[m][n]*tls.Pow(radius-rIntervals[m],n);
				break;
			}
		}
	}
	return phi;
}
// EAM K potential //
// Paper:
// Belashchenko, D.K., 2012.
// Electron Contribution to Energy of Alkali Metals in the Scheme of an Embedded Atom Model.
// High Temp., 50(3), pp.331-339.
// Cut-off radius: 9.57 Angstrom.
double* MC::EAMK_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Particle ithPart, jthPart;
	static double energy[2];
	double rij;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= 9.57){
				rho += EAMK_eDens(rij);
				phiLow += EAMK_PairPot(rij);
			}
		}
	}
	energy[0] = EAMK_EmbPot(rho)*eVToK;
	energy[1] = phiLow*eVToK;
	return energy;
}
double MC::EAMK_EmbPot(double rho){
	Tools tls;
	double rhoIntervals[8] = {1.00000,  0.86000,  0.760000,  0.560000,  0.440000,  0.280000,  1.150000,  0.000000};
	double aValues[8]      = {0.00000, -0.24030, -0.237254, -0.241903, -0.250201, -0.238019, -0.206929, -0.236804};
	double bValues[8]      = {0.00000,  0.00000, -0.043512,  0.136488, -0.053512, -0.149512, -0.239112,  0.046620};
	double cValues[8]      = {0.00000,  0.15540, -0.900000,  0.475000,  0.400000,  0.280000, -1.200000,  0.072000};
	double phi = 0;
	int mValue = 2;

	if (rhoIntervals[1] <= rho && rho < rhoIntervals[6]) phi = aValues[1] + cValues[1]*tls.Pow(rho-rhoIntervals[0],2);
	for (int i=2; i<6; i++){
		if (rhoIntervals[i] <= rho && rho < rhoIntervals[i-1]){
			phi = aValues[i] + bValues[i]*(rho-rhoIntervals[i-1]) + cValues[i]*tls.Pow(rho-rhoIntervals[i-1],2);
			break;
		}
	}
	if (rho < rhoIntervals[5]){
		phi = aValues[6] + bValues[6]*(rho-rhoIntervals[5]) + cValues[6]*tls.Pow(rho-rhoIntervals[5],2);
	}
	if (rho >= rhoIntervals[6]){
		phi = aValues[7] + bValues[7]*(rho-rhoIntervals[6]) + cValues[7]*tls.Pow(rho-rhoIntervals[6], mValue);
	}
	return phi;
}
double MC::EAMK_eDens(double radius){
	double pValues[3] = {0, 3.7461, 0.8400};
	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::EAMK_PairPot(double radius){
	Tools tls;
	double phi=0.;

	if (radius > 3.6){
		phi  = -0.84432914577684e2                   + 0.16516874958198e4/radius            - 0.18050026911001e5/tls.Pow(radius,2);
		phi +=  0.12158970338571e6/tls.Pow(radius,3) - 0.52145119211160e6/tls.Pow(radius,4) + 0.13954378228164e7/tls.Pow(radius, 5);
		phi += -0.21278483968480e7/tls.Pow(radius,6) + 0.14119174994173e7/tls.Pow(radius,7) + 0.18351869172267e1*radius;
	}else phi = 0.183253 + 0.201781*(3.60-radius) + 0.14*(exp(1.96*(3.60-radius)) - 1);
	return phi;
}
// EAM Rb potential //
// Paper:
// Belashchenko, D.K., 2012.
// Electron Contribution to Energy of Alkali Metals in the Scheme of an Embedded Atom Model.
// High Temp., 50(3), pp.331-339.
// Cut-off radius: 14.35 Angstrom.
double* MC::EAMRb_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Particle ithPart, jthPart;
	static double energy[2];
	double rij;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= 14.35){ //Cut-off radius is 14.35 AA.
				rho += EAMRb_eDens(rij);
				phiLow += EAMRb_PairPot(rij);
			}
		}
	}
	energy[0] = EAMRb_EmbPot(rho)*eVToK;
	energy[1] = phiLow*eVToK;
	return energy;
}
double MC::EAMRb_EmbPot(double rho){
	Tools tls;
	double rhoIntervals[8] = {1.00000,  0.90000,  0.760000,  0.650000,  0.500000,  0.280000,  1.450000,  0.000000};
	double aValues[8]      = {0.00000, -0.39630, -0.395189, -0.394587, -0.390822, -0.367575, -0.192569, -0.373802};
	double bValues[8]      = {0.00000,  0.00000, -0.022220,  0.013620, -0.082080, -0.227880, -1.363080,  0.099990};
	double cValues[8]      = {0.00000,  0.11110, -0.128000,  0.435000,  0.486000,  2.580000,  1.200000,  0.045000};
	double mValue=1.5, phi=0;

	if (rhoIntervals[1] <= rho && rho < rhoIntervals[6]) phi = aValues[1] + cValues[1]*tls.Pow(rho-rhoIntervals[0],2);
	for (int i=2; i<6; i++){
		if (rhoIntervals[i] <= rho && rho < rhoIntervals[i-1]){
			phi = aValues[i] + bValues[i]*(rho-rhoIntervals[i-1]) + cValues[i]*tls.Pow(rho-rhoIntervals[i-1],2);
			break;
		}
	}
	if (rho < rhoIntervals[5]){
		phi = aValues[6] + bValues[6]*(rho-rhoIntervals[5]) + cValues[6]*tls.Pow(rho-rhoIntervals[5],2);
	}
	if (rho >= rhoIntervals[6]){
		phi = aValues[7] + bValues[7]*(rho-rhoIntervals[6]) + cValues[7]*pow(rho-rhoIntervals[6], mValue);
	}
	return phi;
}
double MC::EAMRb_eDens(double radius){
	double pValues[3] = {0, 4.4417, 0.8192};
	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::EAMRb_PairPot(double radius){
	Tools tls;
	double phi=0.;

	if (radius > 3.7){
		phi  =  0.31906474916390e2           - 0.81743392122167e3/radius    + 0.11276093229851e5/tls.Pow(radius,2);
		phi += -0.91581958116768e5/tls.Pow(radius,3) + 0.45029175720846e6/tls.Pow(radius,4) - 0.13200477580999e7/tls.Pow(radius,5);
		phi +=  0.21301998697094e7/tls.Pow(radius,6) - 0.14604995446672e7/tls.Pow(radius,7) - 0.51504928921871*radius;
	}else phi = 0.132908 + 0.040299*(3.70-radius) + 0.15*(exp(1.96*(3.70-radius)) - 1);
	return phi;
}
// EAM Ga potential //
// Paper:
// Belashchenko, D.K., 2012.
// Computer Simulation of the Properties of Liquid Metals: Gallium, Lead, and Bismuth.
// Russ. J. Phys. Chem. A, 86, pp.779-790.
// Cut-off radius: 8.3 Angstrom.
double* MC::EAMGa_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Particle ithPart, jthPart;
	static double energy[2];
	double rij;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	for (int j=1; j<=box[ithBox].fluid[jthSpecies].nParts; j++){
		jthPart = box[ithBox].fluid[jthSpecies].particle[j];
		if (ithPart != jthPart){
			rij = NeighDistance(ithBox, ithPart, jthPart);
			if (rij <= 8.3){ //Cut-off radius is 8.3 AA.
				rho += EAMGa_eDens(rij);
				phiLow += EAMGa_PairPot(rij);
			}
		}
	}
	energy[0] = EAMGa_EmbPot(rho)*eVToK;
	energy[1] = phiLow*eVToK;
	return energy;
}
double MC::EAMGa_EmbPot(double rho){
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
double MC::EAMGa_eDens(double radius){
	double pValues[3] = {0, 2.24450, 1.2};
	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::EAMGa_PairPot(double radius){
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
// ---------- Fluid-Fluid potentials ---------- //

// ---------- Solid-Fluid potentials ---------- //
// Paper:
// Steele, W.A., 1973.
// The Physical Interaction of Gases with Crystalline Solids: I. Gas-Solid Energies and Properties of Isolated Adsorbed Atoms.
// Surf. Sci., 36(1), pp.317-352.
// Steele 10-4-3 potential extended by Jason for multiple layers.
// Domain of potential: z in [0,H], where H is the pore size.
double MC::SlitLJ_Pot(int ithBox, int ithSpecies, int index){
	Tools tls;
	double z, t1, t2, t3, t4;
	int nLayers = box[ithBox].nLayersPerWall;
	double sizeR = box[ithBox].width[2]*0.5; // AA
	double dens = box[ithBox].solidDens; // AA^-2
	double eps = box[ithBox].fluid[ithSpecies].epsilon[0]; // K
	double sig = box[ithBox].fluid[ithSpecies].sigma[0]; // AA
	double delta = box[ithBox].deltaLayers; // AA
	double usf = 0;

	z = box[ithBox].fluid[ithSpecies].particle[index].z-sizeR;
	for (int i=0; i<nLayers; i++){
		t1 = tls.Pow(sig/(sizeR+i*delta+z),10);
		t2 = tls.Pow(sig/(sizeR+i*delta-z),10);
		t3 = tls.Pow(sig/(sizeR+i*delta+z),4);
		t4 = tls.Pow(sig/(sizeR+i*delta-z),4);
		usf += 0.2*(t1+t2)-0.5*(t3+t4);
	}
	usf *= 4.*pi*eps*dens*sig*sig; // K
	return usf;
}
// Paper:
// Tjatjopoulos et al., 1988.
// Molecule-Micropore Interaction Potentials.
// J. Phys. Chem., 92(13), pp.4006-4007.
// Source for hypergeometric function:
// Press, W.H., 2007.
// Numerical recipes 3rd edition: The art of scientific computing.
// Cambridge university press.
double MC::CylindricalLJ_Pot(int ithBox, int ithSpecies, int index){
	Tools tls;
	double yPos, zPos;
	double x, u1, u2, usf, rho;
	double sizeR = box[ithBox].width[2]*0.5; // AA
	double dens = box[ithBox].solidDens; // AA^-2
	double eps = box[ithBox].fluid[ithSpecies].epsilon[0]; // K
	double sig = box[ithBox].fluid[ithSpecies].sigma[0]; // AA

	yPos = box[ithBox].fluid[ithSpecies].particle[index].y - 0.5*box[ithBox].width[1];
	zPos = box[ithBox].fluid[ithSpecies].particle[index].z - 0.5*box[ithBox].width[2];
	rho = sqrt(tls.Pow(yPos,2) + tls.Pow(zPos,2));
	// Reducing units/
	dens *= sig*sig;
	sizeR /= sig;
	rho /= sig;
	x = sizeR-rho;
	u1 = 63. /32. * 1./(tls.Pow(x,10)*tls.Pow(2.0-x/sizeR,10)) * hypgeo(-4.5,-4.5,1.,tls.Pow(1-x/sizeR,2));
	u2 = 3. / (tls.Pow(x,4)*tls.Pow(2.0-x/sizeR,4)) * hypgeo(-1.5,-1.5,1.,tls.Pow(1.0-x/sizeR,2));
	usf = (u1-u2)*pi*pi*dens*eps; // K
	return usf;
}
// Paper:
// Baksh, M.S.A. and Yang, R.T., 1991.
// Model for Spherical Cavity Radii and Potential Functions of Sorbates in Zeolites.
// AIChE J., 37(6), pp.923-930.
double MC::SphericalLJ_Pot(int ithBox, int ithSpecies, int index){
	Tools tls;
	double xPos, yPos, zPos;
	double r, x, usf;
	double sizeR = box[ithBox].width[2]*0.5;
	double dens = box[ithBox].solidDens;
	double eps = box[ithBox].fluid[ithSpecies].epsilon[0];
	double sig = box[ithBox].fluid[ithSpecies].sigma[0];
	double u1=0, u2=0;

	xPos = box[ithBox].fluid[ithSpecies].particle[index].x - 0.5*box[ithBox].width[0];
	yPos = box[ithBox].fluid[ithSpecies].particle[index].y - 0.5*box[ithBox].width[1];
	zPos = box[ithBox].fluid[ithSpecies].particle[index].z - 0.5*box[ithBox].width[2];
	r = sqrt(tls.Pow(xPos,2) + tls.Pow(yPos,2) + tls.Pow(zPos,2));
	// Reducing units/
	dens *= sig*sig;
	sizeR /= sig;
	r /= sig;
	x = sizeR-r;
	for (int i=0; i<10; i++){
		u1 += 1. / tls.Pow(sizeR,i) / tls.Pow(x,10-i);
		u1 += tls.Pow(-1.0,i) / tls.Pow(sizeR,i) / tls.Pow(x-2.0*sizeR,10-i);
		if (i < 4){
			u2 += 1. / tls.Pow(sizeR,i) / tls.Pow(x,4-i);
			u2 += tls.Pow(-1.0,i) / tls.Pow(sizeR,i) / tls.Pow(x-2.0*sizeR,4-i);
		}
	}
	usf = (2./5.) * u1 - u2;
	usf *= 2.*pi*eps*dens; //K
	return usf;
}
// ---------- Solid-Fluid potentials ---------- //

