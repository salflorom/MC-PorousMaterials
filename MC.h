//Author: Santiago A. Flores Roman

#ifndef _MC_H_
#define _MC_H_

#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

#define NBINS 100
#define MAXPART 99999 //Max. num. of particles.
#define pi M_PI
#define kb 1.38062e-23 // J/K
#define na 6.02214e23 // mol^-1
#define plank 1.054571817e-34 // J*s

using namespace std;

/* **************************************************************************** */
class MC {
/* **************************************************************************** */
	public:
/* **************************************************************************** */
		MC(void);
		void ReadInputFile(string inFileName);
		void PrintParams(void);
		void ComputeBoxSize(void);
		void AssignPossitions(void);
		void MoveParticle(int index);
		void PBC(int index);
		void RDF(void);
		void PrintRDF(int step);
		void Metropolis(int index);
		void ResetParticle(int index);
		void CreateEXYZ(int step);
		void CreateLogFile(int step);
		void PrintStats(int step);
		void IncrementObstruction(void){stats.obstruction++;}
		double SystemEnergy(void);
		void ComputeChemicalPotential(void);
		double EnergyOfParticle(int index);
		double LJ_Energy(int i, int j);
		// For EAM Ga potential //
		double EAMGA_Energy(int index);
		double StepUnit(double radius, double leftLim, double rightLim);
		double EmbPot(double rho);
		double eDens(double radius);
		double PairPot(double radius);
		// For EAM Ga potential //
		double NeighDistance(int i, int j);
		double Random(void){return (double)rand()/RAND_MAX;} //random in the interval [0,1].
		bool Overlap(int index);
		int SelectParticle(void);
		int GetNSteps(void){return int(sim.nSteps);}
		int GetNInitSteps(void){return int(sim.nInitSteps);}
		int GetNEquilSteps(void){return int(sim.nEquilSteps);}
		int GetPrintEvery(void){return int(sim.printEvery);}
		int GetObstruction(void){return stats.obstruction;}
		int GetObstructed(int index){return part[index].obstructed;}
		string* GetCompute(void){return sim.compute;}
/* **************************************************************************** */
	protected:
/* **************************************************************************** */
		struct Stats{
			int acceptance, rejection, obstruction, nDisplacements;
			int widomInsertions;
			double binWidth, insertionParam;
			long int rdfx[NBINS+1], rdfy[NBINS+1], rdfz[NBINS+1];
		} stats;
		struct Simulation{
			string compute[3];
			double dx, dy, dz;
			double nSteps, nInitSteps, nEquilSteps, printEvery;
			int step;
			bool PBC[3];
		} sim;
		struct Fluid{
			string name, vdwPot;
			double temp, dens, molarMass, boxWidth, rcut, sigma, epsilon;
			double initialE, finalE, mu;
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

