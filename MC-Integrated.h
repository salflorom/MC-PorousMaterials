#ifndef _MC_INTEGRATED_H_
#define _MC_INTEGRATED_H_

#include <climits>

#define NBINS 100
#define MAXPART INT_MAX //Max. num. of particles.
#define pi 3.14159265358979323846
#define kb 1.38062e-23
#define na 6.02214e23

using namespace std;

void ReadInputFile(void);
void PrintInputFile(void);
void ComputeBoxSize(void);
void AssignPossitions(void);
int SelectParticle(void);
void MoveParticle(int index);
void PBC(int index);
void RDF(void);
void PrintRDF(void);
void Metropolis(int index);
double EnergyChange(int index);
double Distance(int i, int j);
void ResetParticle(int index);
void CreateEXYZ(int step);
void PrintStats(int step);
double Random(void){return (double)rand()/RAND_MAX;} //random in the interval [0,1].

struct Stats{
	double binWidth;
	int acceptance=0, rejection=0, obstruction=0;
	long int rdfx[NBINS+1]={0}, rdfy[NBINS+1]={0}, rdfz[NBINS+1]={0};
} stats;
struct Simulation{
	int step, nSteps, nEquilSteps, printEvery;
	int obstructionLimit=10;
	double dx=1, dy=1, dz=1;
	bool PBC[3]={true, true, true};
	FILE *file;
	char command[35];
} sim;
struct Fluid{
	string symbol;
	double temp, dens, molarMass, eTotal, boxWidth, rcut, sigma, epsilon;
	int nParts;
} fluid;
struct Particle{
	double x, y, z;
	int obstructed=0;
} part[INT_MAX];

#endif /* _MC_INTEGRATED_H_ */
