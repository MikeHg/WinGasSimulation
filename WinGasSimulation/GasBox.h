#pragma once
// the number of particles in the box
#define NUM_PARTS 32 * 32

// Time step for verlet algorithm
#define DELTA_T 0.0005L/3.0L

// Standard Mass of Atoms
#define MASS 0.1L
// how much heavier is atom 1
#define MASS1_FACTOR 20.0L

// the size of the box where the particles are confined to
#define BOX_X 640.0L
#define BOX_Y 480.0L
#define BOX_Z 500.0L

// number of intervals to partition the Maxwell velocity distribution (see V_MAX)
#define NUM_BOXES 20
// maximum velocity we want to consider for Maxwell velocity distribution
#define V_MAX 500.0L

// initial minimal distance between particles
#define INI_MIN_DIST 5.0L

// medium start velocity for initial Gauss distribution
#define V_START 205.0L

#define NUM_THREADS 6

//#define SINGLE

#include "Particle.h"
#include "BS_thread_pool.hpp"
#include <cstddef>
#include <vector>

class GasBox {
private:
	long double initialEnergy_;
	long double totalEnergy_;
	long double relativeEnergyError_;
	long double fx_pt[NUM_PARTS], fy_pt[NUM_PARTS], fz_pt[NUM_PARTS];
	Particle* myParticles = new Particle[NUM_PARTS];

	long double tempPotEnergy[NUM_PARTS]; // can be lowered to NUM_THREADS
	long double fx_tmp[NUM_THREADS][NUM_PARTS], fy_tmp[NUM_THREADS][NUM_PARTS], fz_tmp[NUM_THREADS][NUM_PARTS];

	void Gauss(long double* x1, long double* x2);
	long double Forces();
	long double ForcesSingleThread();
	void periodic_bounds(int i);
	unsigned int stepNo;

	// for maxwell drawing
	double MaxwPoints[NUM_BOXES][2];
	double firstControlPoints[NUM_BOXES][2];
	double secondControlPoints[NUM_BOXES][2];
	void initMaxwell();
	int setDistribV();

	int calcDistribution[NUM_BOXES];

	BS::thread_pool pool = BS::thread_pool (NUM_THREADS);

public:
	GasBox();

	void verlet();
	void switchVelocities();
	long double getTotalEnergy();
	long double getEnergyError();
	int getNumberOfParticles();
	unsigned int getStep();
	double getMaxwell(int point, int coordinate);
	double getFirstCP(int point, int coordinate);
	double getSecondCP(int point, int coordinate);
	int getDistribV(int n);

	std::vector<double> getParticleLocation(int i);
	std::vector<double> getParticleVelocity(int i);
	std::vector<double> getParticleForces(int i);
};
