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
#define kb 1.380649e-23 // J/K
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
		void ResetStats(void);
		void PrintParams(void);
		void ComputeBoxSize(void);
		void InitialConfig(void);
		void MoveParticle(void);
		void ChangeVolume(void);
		void ExchangeParticle(void);
		void PBC(int index);
		void RDF(void);
		void PrintRDF(int step);
		void Metropolis(int index);
		void ResetParticle(int index);
		void CreateEXYZ(int step);
		void CreateLogFile(int step);
		void PrintStats(int step);
		void ComputeWidom(void);
		void ComputeChemicalPotential(void);
		void IncrementObstruction(void){stats.obstruction++;}
		double SystemEnergy(void);
		double EnergyOfParticle(int index);
		double LJ_Energy(int i, int j);
		// For EAM Ga potential vvvvv //
		double EAMGA_Energy(int index);
		double StepUnit(double radius, double leftLim, double rightLim);
		double EmbPot(double rho);
		double eDens(double radius);
		double PairPot(double radius);
		// For EAM Ga potential ^^^^^ //
		double NeighDistance(int i, int j);
		double Random(void){return (double)rand()/RAND_MAX;} //random in the interval [0,1].
		double* GetMCMoveProbabilities(void);
		int GetNSets(void){return int(sim.nSets);}
		int GetNEquilSets(void){return int(sim.nEquilSets);}
		int GetNStepsPerSet(void){return int(sim.nStepsPerSet);}
		int GetPrintEvery(void){return int(sim.printEvery);}
		int GetObstruction(void){return stats.obstruction;}
		int GetObstructed(int index){return part[index].obstructed;}
		string* GetCompute(void){return sim.compute;}
/* **************************************************************************** */
	protected:
/* **************************************************************************** */
		struct Stats{
			int acceptance, rejection, obstruction, nDisplacements;
			int acceptanceVol, rejectionVol, nVolChanges;
			int insertion, deletion, nExchanges;
			int widomInsertions;
			double binWidth, widom;
			long int rdfx[NBINS+1], rdfy[NBINS+1], rdfz[NBINS+1];
		} stats;
		struct Simulation{
			string compute[3];
			double dx, dy, dz, dv;
			double nSets, nEquilSets, nStepsPerSet, printEvery;
			double displaceProb, exchangeProb, volumeProb;
			int step;
			bool PBC[3];
		} sim;
		struct Fluid{
			string name, vdwPot;
			double temp, dens, molarMass, extPress, boxWidth, mu;
			double rcut, sigma, epsilon;
			double systemEnergy;
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

