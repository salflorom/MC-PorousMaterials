//Author: Santiago A. Flores Roman

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
// EAM Ga potential //
// Paper:
// Belashchenko, D.K., 2012.
// Computer Simulation of the Properties of Liquid Metals: Gallium, Lead, and Bismuth.
// Russ. J. Phys. Chem. A, 86, pp.779-790.
// Cut=off radius: 8.3 Angstrom.
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
double MC::EAMGa_StepUnit(double radius, double leftLim, double rightLim){
	if (leftLim <= radius && radius < rightLim) return 1.;
	else return 0.;
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
				phi += aValues[m][i] * pow((radius-rIntervals[i+1]),m) * EAMGa_StepUnit(radius, rIntervals[i], rIntervals[i+1]);
			}
		}
	}else if (rIntervals[0] < radius && radius <= rIntervals[1]){
		phi = 0.619588 - 51.86268*(2.15-radius) + 27.8*(exp(1.96*(2.15-radius))-1);
	}
	return phi;
}
// EAM Rb potential //
// Paper:
// Belashchenko, D.K., 2006.
// Embedded Atom Model Application to Liquid Metals: Liquid rubidium.
// Russ. J. Phys. Chem., 80(10), pp.1567-1577.
// Available interval (K): (313, 2000).
// Cut-off radius: 14.35 Angstrom.
double* MC::EAMRb_Pot(int ithBox, int ithSpecies, int jthSpecies, int index){
	Particle ithPart, jthPart;
	static double energy[2];
	double rij, rcut;
	double rho=0, phiLow=0;
	double eVToK = (801088317./5.0e27)/kb;

	ithPart = box[ithBox].fluid[ithSpecies].particle[index];
	rcut = fluid[ithSpecies].rcut[jthSpecies]; // AA
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
	double rhoIntervals[2] = { 1.000000,  0.900000};
	double aValues[4]      = { 0.000000, -0.406100,  0.111100, -0.008000};
	double bValues[4]      = {-0.404981, -0.022460, -0.091400, -0.860667};
	double temp = thermoSys.temp;
	double phi = 0;

	if (rho >= rhoIntervals[1]) phi = aValues[1] + aValues[2]*tls.Pow(rho-rhoIntervals[0],2) + aValues[3]*tls.Pow(rho-rhoIntervals[0],3);
	else phi = bValues[0] + bValues[1]*(rho-rhoIntervals[1]) + bValues[2]*tls.Pow(rho-rhoIntervals[1],2) + bValues[3]*(rho-rhoIntervals[1],3);
	return phi;
}
double MC::EAMRb_eDens(double radius){
	double pValues[3] = {0, 4.4400, 0.8192};
	double temp = thermoSys.temp;

	return pValues[1]*exp(-pValues[2]*radius);
}
double MC::EAMRb_PairPot(double radius){
	int ithRadius;
	double schommers[227][2] = {{3.100,  3.4941300e+00},
								{3.150,  3.1969900e+00},
								{3.200,  2.9268800e+00},
								{3.250,  1.0476900e+00},
								{3.300,  4.3920000e-01},
								{3.350,  4.0885000e-01},
								{3.400,  2.5805000e-01},
								{3.450,  2.4818000e-01},
								{3.500,  2.3567000e-01},
								{3.550,  2.0578000e-01},
								{3.600,  1.7769100e-01},
								{3.650,  1.8910700e-01},
								{3.700,  1.3265600e-01},
								{3.750,  1.1538400e-01},
								{3.800,  1.1373200e-01},
								{3.850,  7.2463100e-02},
								{3.900,  7.1904900e-02},
								{3.950,  5.0585800e-02},
								{4.000,  4.3854200e-02},
								{4.050,  3.6086500e-02},
								{4.100,  2.4898800e-02},
								{4.150,  1.3997300e-02},
								{4.200,  3.8783100e-02},
								{4.250, -5.3519100e-02},
								{4.300, -1.3902700e-02},
								{4.350, -2.2013100e-02},
								{4.400, -2.9561700e-02},
								{4.450, -3.7592700e-02},
								{4.500, -4.4626000e-02},
								{4.550, -5.0976600e-02},
								{4.600, -5.6596800e-02},
								{4.650, -6.0931100e-02},
								{4.700, -6.5009900e-02},
								{4.750, -6.8162600e-02},
								{4.800, -7.1236600e-02},
								{4.850, -7.3332400e-02},
								{4.900, -7.5456500e-02},
								{4.950, -7.6865200e-02},
								{5.000, -7.7807600e-02},
								{5.050, -7.8660000e-02},
								{5.100, -7.9080200e-02},
								{5.150, -7.9317100e-02},
								{5.200, -7.9670000e-02},
								{5.250, -7.9706300e-02},
								{5.300, -7.9241500e-02},
								{5.350, -7.8614900e-02},
								{5.400, -7.8224900e-02},
								{5.450, -7.7142900e-02},
								{5.500, -7.6468700e-02},
								{5.550, -7.5162300e-02},
								{5.600, -7.3925900e-02},
								{5.650, -7.2673800e-02},
								{5.700, -7.1144000e-02},
								{5.750, -6.9688000e-02},
								{5.800, -6.8014100e-02},
								{5.850, -6.6102400e-02},
								{5.900, -6.4178600e-02},
								{5.950, -6.2284600e-02},
								{6.000, -5.9969600e-02},
								{6.050, -5.7450800e-02},
								{6.100, -5.5560500e-02},
								{6.150, -5.3395500e-02},
								{6.200, -5.1058600e-02},
								{6.250, -4.8689900e-02},
								{6.300, -4.6535000e-02},
								{6.350, -4.4248000e-02},
								{6.400, -4.2408300e-02},
								{6.450, -4.0071100e-02},
								{6.500, -3.8135600e-02},
								{6.550, -3.6491700e-02},
								{6.600, -3.3512900e-02},
								{6.650, -3.2244700e-02},
								{6.700, -3.0038100e-02},
								{6.750, -2.7806600e-02},
								{6.800, -2.5773500e-02},
								{6.850, -2.4274600e-02},
								{6.900, -2.2814200e-02},
								{6.950, -2.1176500e-02},
								{7.000, -1.9699200e-02},
								{7.050, -1.7650100e-02},
								{7.100, -1.6637200e-02},
								{7.150, -1.5166500e-02},
								{7.200, -1.3816800e-02},
								{7.250, -1.2044100e-02},
								{7.300, -1.0851200e-02},
								{7.350, -9.4040900e-03},
								{7.400, -8.2554000e-03},
								{7.450, -6.9271000e-03},
								{7.500, -5.9576500e-03},
								{7.550, -4.4748000e-03},
								{7.600, -3.4041900e-03},
								{7.650, -2.4339400e-03},
								{7.700, -1.1721200e-03},
								{7.750, -1.8191700e-04},
								{7.800,  5.7357900e-04},
								{7.850,  1.7539600e-03},
								{7.900,  2.5720000e-03},
								{7.950,  3.6192400e-03},
								{8.000,  4.3219700e-03},
								{8.050,  5.0186600e-03},
								{8.100,  5.8141600e-03},
								{8.150,  6.2828800e-03},
								{8.200,  6.9115300e-03},
								{8.250,  7.7310300e-03},
								{8.300,  8.1589900e-03},
								{8.350,  8.5557300e-03},
								{8.400,  8.8646100e-03},
								{8.450,  9.0379500e-03},
								{8.500,  9.7870000e-03},
								{8.550,  9.7653700e-03},
								{8.600,  9.7845800e-03},
								{8.650,  9.9648200e-03},
								{8.700,  9.9365600e-03},
								{8.750,  1.0061700e-02},
								{8.800,  1.0227500e-02},
								{8.850,  9.7450200e-03},
								{8.900,  1.0036300e-02},
								{8.950,  9.7751900e-03},
								{9.000,  9.6744300e-03},
								{9.050,  9.4183200e-03},
								{9.100,  9.5039600e-03},
								{9.150,  9.0345800e-03},
								{9.200,  8.8524100e-03},
								{9.250,  8.7414900e-03},
								{9.300,  8.3813500e-03},
								{9.350,  8.1424700e-03},
								{9.400,  7.7455200e-03},
								{9.450,  7.3538000e-03},
								{9.500,  7.0507100e-03},
								{9.550,  6.4152900e-03},
								{9.600,  6.5860400e-03},
								{9.650,  6.2360200e-03},
								{9.700,  6.0886700e-03},
								{9.750,  5.7449800e-03},
								{9.800,  5.6021200e-03},
								{9.850,  5.1657900e-03},
								{9.900,  5.0184900e-03},
								{9.950,  4.6897100e-03},
								{10.00,  4.4517000e-03},
								{10.05,  4.0207800e-03},
								{10.10,  3.5952500e-03},
								{10.15,  3.4573800e-03},
								{10.20,  3.2455500e-03},
								{10.25,  2.8338800e-03},
								{10.30,  2.6146600e-03},
								{10.35,  2.3418200e-03},
								{10.40,  1.7358700e-03},
								{10.45,  1.6351600e-03},
								{10.50,  1.1644800e-03},
								{10.55,  1.0370600e-03},
								{10.60,  6.9196700e-04},
								{10.65,  3.9450400e-04},
								{10.70,  2.9390300e-05},
								{10.75, -2.3650500e-04},
								{10.80, -7.2125300e-04},
								{10.85, -8.0213800e-04},
								{10.90, -1.0311700e-03},
								{10.95, -1.2667200e-03},
								{11.00, -1.6033000e-03},
								{11.05, -1.6252600e-03},
								{11.10, -2.0686100e-03},
								{11.15, -2.1527400e-03},
								{11.20,  2.3679000e-03},
								{11.25, -2.5607600e-03},
								{11.30, -2.9458900e-03},
								{11.35, -2.6662500e-03},
								{11.40, -2.8480300e-03},
								{11.45, -2.8798400e-03},
								{11.50, -2.9133600e-03},
								{11.55, -3.0798100e-03},
								{11.60, -3.1140600e-03},
								{11.65, -3.1287200e-03},
								{11.70, -3.0777100e-03},
								{11.75, -3.3027200e-03},
								{11.80, -3.0491700e-03},
								{11.85, -3.4808900e-03},
								{11.90, -3.0896200e-03},
								{11.95, -3.1461600e-03},
								{12.00, -3.1636700e-03},
								{12.05, -3.1591900e-03},
								{12.10, -3.1163800e-03},
								{12.15, -3.1804400e-03},
								{12.20, -2.9218900e-03},
								{12.25, -3.0356500e-03},
								{12.30, -3.0122500e-03},
								{12.35, -3.0819500e-03},
								{12.40, -2.9736600e-03},
								{12.45, -2.8802200e-03},
								{12.50, -2.9095500e-03},
								{12.55, -2.6343300e-03},
								{12.60, -2.7071100e-03},
								{12.65, -2.5578300e-03},
								{12.70, -2.4176500e-03},
								{12.75, -2.1604100e-03},
								{12.80, -2.3342100e-03},
								{12.85, -1.9505500e-03},
								{12.90, -2.1415200e-03},
								{12.95, -1.9726900e-03},
								{13.00, -1.7187300e-03},
								{13.05, -1.8171600e-03},
								{13.10, -1.7118400e-03},
								{13.15, -1.6313200e-03},
								{13.20, -1.4445300e-03},
								{13.25, -1.4670200e-03},
								{13.30, -1.2981000e-03},
								{13.35, -1.4856900e-03},
								{13.40, -1.1247700e-03},
								{13.45, -1.1761600e-03},
								{13.50, -1.0486000e-03},
								{13.55, -8.7051200e-04},
								{13.60, -7.7697900e-04},
								{13.65, -7.5154700e-04},
								{13.70, -8.1809500e-04},
								{13.75, -5.9900800e-04},
								{13.80, -6.8876100e-04},
								{13.85, -4.6956400e-04},
								{13.90, -4.3730600e-04},
								{13.95, -3.8378500e-04},
								{14.00, -3.0979800e-04},
								{14.05, -3.0162700e-04},
								{14.10, -2.6248900e-04},
								{14.15, -2.3979500e-04},
								{14.20, -1.7069100e-04},
								{14.25, -1.3420600e-04},
								{14.30,  1.0163200e-04},
								{14.35,  1.0163200e-04},
								{14.40,  1.0163200e-04}};
	double phi=0.;

	// y = y1 + (x-x1)(y2-y1)/(x2-x1)
	if (radius < 3.1) phi = 1e10;
	else{
		ithRadius = int((radius-3.10)/0.05);
		phi = schommers[ithRadius][1] + (radius-schommers[ithRadius][0])*(schommers[ithRadius+1][1]-schommers[ithRadius][1])/(schommers[ithRadius+1][0]-schommers[ithRadius][0]);
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

