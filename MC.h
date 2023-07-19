#ifndef _MC_H_
#define _MC_H_

#include <string>
#include <cstdlib>
#include <unordered_map>
#include <vector>

#define NBINS 100
#define MAXPART 99999 //Max. num. of particles.
#define pi 3.14159265358979323846
#define kb 1.38062e-23
#define na 6.02214e23

using namespace std;

/* **************************************************************************** */
class MC {
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		MC(void);
		int SelectParticle(void);
		void ReadInputFile(string inFileName);
		void PrintParams(void);
		void ComputeBoxSize(void);
		void AssignPossitions(void);
		void MoveParticle(int index);
		void AdjustStepSize(void);
		void PBC(int index);
		void RDF(void);
		void PrintRDF(void);
		void Metropolis(int index);
		void ResetParticle(int index);
		void CreateEXYZ(int step);
		void PrintStats(int step);
		double EnergyOfParticle(int index);
		double Distance(int i, int j);
		unordered_map<string,double> GetFluidParams(void);
		vector<vector<double>> GetPositions(void);
		int GetNSteps(void){return int(sim.nSteps);}
		int GetNEquilSteps(void){return int(sim.nEquilSteps);}
		int GetPrintEvery(void){return int(sim.printEvery);}
		int GetObstructionLimit(void){return sim.obstructionLimit;}
		int GetObstruction(void){return stats.obstruction;}
		int GetObstructed(int index){return part[index].obstructed;}
		void IncrementObstruction(void){stats.obstruction++;}
		double Random(void){return (double)rand()/RAND_MAX;} //random in the interval [0,1].

	protected:
/* **************************************************************************** */
		struct Stats{
			int acceptance, rejection, obstruction, nDisplacements;
			double binWidth;
			long int rdfx[NBINS+1], rdfy[NBINS+1], rdfz[NBINS+1];
		} stats;
		struct Simulation{
			double dx, dy, dz;
			double nSteps, nEquilSteps, printEvery;
			int step;
			int obstructionLimit;
			bool PBC[3];
		} sim;
		struct Fluid{
			string symbol;
			double temp, dens, molarMass, boxWidth, rcut, sigma, epsilon;
			double initialE, finalE;
			int nParts;
		} fluid;
		struct Particle{
			double x, y, z;
			int obstructed;
		} part[MAXPART];
/* **************************************************************************** */
};
/* **************************************************************************** */

#endif /* _MC_H_ */

