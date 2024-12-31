#include "pch.h"
#include "GasBox.h"
#include <math.h>
#include <stdlib.h>
#include "Bezier.h"

//#include <chrono>
//#include <iostream>
//#include <thread>

/*
* initializes all coordinates and velocities of the NUM_PARTS particles
*
* The first particle is the special particle ("Aufteilchen") with a different mass
* then all other particles,
*
* For simplicity the z-location is the same for all particles BOX_Z / 2.0L
* and the z-velocity is zero, i.e. no particle will move out of the plane,
*
* All calculations are done 3-dimensional. I.e. initializing "z" also would be ok.
*
*/
GasBox::GasBox() {
	int i, j, distflag;
	long double totalKineticEnergy, totalPotentialEnergy, tmp_rnd[3], r;
	long double pi = 3.1415926535897932384626433832795028841971L;

	totalKineticEnergy = 0.0L;
	// starting without seed, i.e. srand(1)

	// handle the first atom differently, i.e. place it on a fixed location with a fixed velocity
	myParticles[0].setLocation(BOX_X / 2.0L - 30.0L, BOX_Y / 2.0L + 25.0L, BOX_Z / 2.0L);
	myParticles[0].setVelocity(150.0L, 0.0L, 0.0L);
	myParticles[0].setMass(MASS * MASS1_FACTOR);
	

	totalKineticEnergy += myParticles[0].getKineticEnergy();

	for (i = 1; i < NUM_PARTS; i++)
	{
		do {
			distflag = 0;

			tmp_rnd[0] = (long double)rand() / (long double)RAND_MAX;
			tmp_rnd[0] = tmp_rnd[0] * (BOX_X - 2.0L) + 1.0L;

			tmp_rnd[1] = (long double)rand() / (long double)RAND_MAX;
			tmp_rnd[1] = tmp_rnd[1] * (BOX_Y - 2.0L) + 1.0L;

			//tmp_rnd[2] = (long double)rand() / (long double)RAND_MAX;
			//z_pt[i] = tmp_rnd * (BOX_Z - 2.0L) + 1.0L;
			tmp_rnd[2] = BOX_Z / 2.0L;

			// ensure that distance between particles is at least INI_MIN_DIST
			for (j = 0; j < i; j++) {
				r = myParticles[j].getDistance(tmp_rnd[0], tmp_rnd[1], tmp_rnd[2]);
				if (r < INI_MIN_DIST) {
					distflag = 1;
					break;
				}
			}
		} while (distflag == 1);

		myParticles[i].setLocation(tmp_rnd[0], tmp_rnd[1], tmp_rnd[2]);
		myParticles[i].setMass(MASS);

		// One could either start with a random velocity
		//tmp_rnd = 2.0 * pi * (long double)rand()/(double)RAND_MAX;
		//vx_pt[i] = V_START * sin(tmp_rnd);
		//vy_pt[i] = V_START * cos(tmp_rnd);

		// or use a normal distribution
		Gauss(&tmp_rnd[0], &tmp_rnd[1]);


		tmp_rnd[0] *= sqrt(V_START) * pi * pi;
		tmp_rnd[1] *= sqrt(V_START) * pi * pi;
		tmp_rnd[2] = 0.0L;

		myParticles[i].setVelocity(tmp_rnd[0], tmp_rnd[1], tmp_rnd[2]);

		totalKineticEnergy += myParticles[i].getKineticEnergy();
	}

#ifdef SINGLE
	totalPotentialEnergy = ForcesSingleThread();
#else
	totalPotentialEnergy = Forces();
#endif // SINGLE

	initialEnergy_ = totalKineticEnergy + totalPotentialEnergy;
	
	// initialize GasBox Energies
	totalEnergy_ = initialEnergy_;
	relativeEnergyError_ = 0.0L;

	stepNo = 0;

	initMaxwell();
	setDistribV();
	switchVelocities();
}


/*
* Implementation of the Marsaglia Polar Method
* see https://en.wikipedia.org/wiki/Marsaglia_polar_method
* Derives two independent random values for a normal distribution
* Mean Value = 0.0
* Std Deviation = 1.0
*
* for one value it should be similar to:
* #include <random>
*
*    std::default_random_engine generator;
*    std::normal_distribution<long double> distribution(0.0, 1.0);
*    long double rValue = distribution(generator);
*/
void GasBox::Gauss(long double* x1, long double* x2)
{
	long double u1, u2, v1, v2, s;

	s = 2.0L;

	while (s >= 1.0L)
	{
		u1 = (long double)rand() / (long double)RAND_MAX;
		u2 = (long double)rand() / (long double)RAND_MAX;

		v1 = 2.0L * u1 - 1.0L;
		v2 = 2.0L * u2 - 1.0L;
		s = v1 * v1 + v2 * v2;
	}

	*x1 = v1 * sqrtl(-2.0L * logl(s) / s);
	*x2 = v2 * sqrtl(-2.0L * logl(s) / s);
}

// calculate forces and return potential energy
// use a single thread approach which is clearly structured but a little slow
long double GasBox::ForcesSingleThread() {
	int i, j;
	long double fx, fy, fz;
	long double pot;

	
	for (i = 0; i < NUM_PARTS; i++)
	{
		fx_pt[i] = 0.0L;
		fy_pt[i] = 0.0L;
		fz_pt[i] = 0.0L;
	}

	pot = 0.0L;

	for (i = 0; i < NUM_PARTS; i++)
	{
		for (j = i + 1; j < NUM_PARTS; j++)
		{
			pot+= myParticles[i].getForceAndPotEnergy(myParticles[j].getX(), myParticles[j].getY(), myParticles[j].getZ(),
				&fx, &fy, &fz);
			fx_pt[i] += fx;
			fy_pt[i] += fy;
			fz_pt[i] += fz;

			fx_pt[j] -= fx;
			fy_pt[j] -= fy;
			fz_pt[j] -= fz;
		}
	}
	
	return pot;
}

// calculate forces and return potential energy
// Use the Multi-Thread Approach to speed-up calculation
long double GasBox::Forces() {
	int i, j;
	long double pot;

	for (i = 0; i < NUM_PARTS; i++)
	{
		fx_pt[i] = 0.0L;
		fy_pt[i] = 0.0L;
		fz_pt[i] = 0.0L;
		tempPotEnergy[i] = 0.0L;
		for (j = 0; j < NUM_THREADS; j++) {
			fx_tmp[j][i] = 0.0L; fy_tmp[j][i] = 0.0L; fz_tmp[j][i] = 0.0L;
		}
	}

	pot = 0.0L;

	// Thread Pool implementation starts here
	const BS::multi_future<void> loop_future = pool.submit_loop<unsigned int>(0, NUM_THREADS,
		[myParticles = myParticles, &tempPotEnergy = tempPotEnergy, &fx_pt = fx_pt, &fy_pt = fy_pt, &fz_pt = fz_pt,
		&fx_tmp = fx_tmp, &fy_tmp = fy_tmp, &fz_tmp = fz_tmp](const unsigned int k) {

			long double fx, fy, fz, r2, r2t, r6, _tmp, x1, y1, z1, x2, y2, z2;

			for (unsigned int m = k; m < NUM_PARTS; m += NUM_THREADS) {
				x1 = myParticles[m].getX(); y1 = myParticles[m].getY(); z1 = myParticles[m].getZ();
				for (unsigned int n = m + 1; n < NUM_PARTS; n++)
				{
					x2 = myParticles[n].getX(); y2 = myParticles[n].getY(); z2 = myParticles[n].getZ();
					r2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
					r2t = 1.0 / r2;
					r2 = (SIGMA * SIGMA) * r2t;

					r6 = r2 * r2 * r2;
					_tmp = 6.0L * E4 * r6 * (2.0L * r6 - 1.0L) * r2t;

					fx = _tmp * (x1 - x2); fy = _tmp * (y1 - y2); fz = _tmp * (z1 - z2);
				
					tempPotEnergy[m] += E4 * r6 * (r6 - 1.0L);

					fx_pt[m] += fx; fy_pt[m] += fy; fz_pt[m] += fz;

					// This part in particular needs to be thread safe
					fx_tmp[k][n] -= fx; fy_tmp[k][n] -= fy; fz_tmp[k][n] -= fz;
				}
			}
		});

    loop_future.wait();

	for (i = 0; i < NUM_PARTS; i++) {
		pot += tempPotEnergy[i];
		for (j = 0; j < NUM_THREADS; j++) {
			fx_pt[i] += fx_tmp[j][i]; fy_pt[i] += fy_tmp[j][i]; fz_pt[i] += fz_tmp[j][i];
		}
	}
	
	return pot;
}

/*
* check if a particle moves out of the box limits and if so,
* change velocities and position. I.e. reflect the particle on the wall
*/
void GasBox::periodic_bounds(int i)
{
	if (myParticles[i].getX() < 0.0L)
	{
		myParticles[i].setLocation(-myParticles[i].getX(), myParticles[i].getY(), myParticles[i].getZ());
		myParticles[i].setVelocity(-myParticles[i].getVx(), myParticles[i].getVy(), myParticles[i].getVz());
	}
	if (myParticles[i].getX() > BOX_X)
	{
		myParticles[i].setLocation(2.0L * BOX_X - myParticles[i].getX(), myParticles[i].getY(), myParticles[i].getZ());
		myParticles[i].setVelocity(-myParticles[i].getVx(), myParticles[i].getVy(), myParticles[i].getVz());
	}

	if (myParticles[i].getY() < 0.0L)
	{
		myParticles[i].setLocation(myParticles[i].getX(), -myParticles[i].getY(), myParticles[i].getZ());
		myParticles[i].setVelocity(myParticles[i].getVx(), -myParticles[i].getVy(), myParticles[i].getVz());
	}
	if (myParticles[i].getY() > BOX_Y)
	{
		myParticles[i].setLocation(myParticles[i].getX(), 2.0L * BOX_Y - myParticles[i].getY(), myParticles[i].getZ());
		myParticles[i].setVelocity(myParticles[i].getVx(), -myParticles[i].getVy(), myParticles[i].getVz());
	}

	if (myParticles[i].getZ() < 0.0L)
	{
		myParticles[i].setLocation(myParticles[i].getX(), myParticles[i].getY(), -myParticles[i].getZ());
		myParticles[i].setVelocity(myParticles[i].getVx(), myParticles[i].getVy(), -myParticles[i].getVz());
	}
	if (myParticles[i].getZ() > BOX_Z)
	{
		myParticles[i].setLocation(myParticles[i].getX(), myParticles[i].getY(), 2.0L * BOX_Z - myParticles[i].getZ());
		myParticles[i].setVelocity(myParticles[i].getVx(), myParticles[i].getVy(), -myParticles[i].getVz());
	}
}

/*
*
* Verlet-Algorithm to calculate new positions and velocities of the particles
* An improved version is Beeman's Algorithm
*
* Verlet:
* r_{n+1} = r_{n} + v_{n} dt + 1/2 a_{n} dt²
*
* v_{n+1} = v_{n} + 1/2 (a_{n} + a_{n+1}) dt
*
* returns
* 1 if the particles have not been initialized yet
* 2 if velocities will be switched in a step
* 0 if the algorithm was normally executed
*/
void GasBox::verlet() 
{
	unsigned int i;
	long double totalPotentialEnergy, totalKineticEnergy, tmp_en;
	long double _dt_2m, _dt2_2m;

	// Could also be parallelized but does not bring any measurable performance increase

	for (i = 0; i < NUM_PARTS; i++) {
		_dt2_2m = DELTA_T * DELTA_T / (2.0L * myParticles[i].getMass());
		_dt_2m = DELTA_T / (2.0L * myParticles[i].getMass());

		myParticles[i].setLocation(
			myParticles[i].getX() + myParticles[i].getVx() * DELTA_T + fx_pt[i] * _dt2_2m,
			myParticles[i].getY() + myParticles[i].getVy() * DELTA_T + fy_pt[i] * _dt2_2m,
			myParticles[i].getZ() + myParticles[i].getVz() * DELTA_T + fz_pt[i] * _dt2_2m
		);

		myParticles[i].setVelocity(
			myParticles[i].getVx() + fx_pt[i] * _dt_2m,
			myParticles[i].getVy() + fy_pt[i] * _dt_2m,
			myParticles[i].getVz() + fz_pt[i] * _dt_2m
		);

		periodic_bounds(i);
	}

#ifdef SINGLE
	totalPotentialEnergy = ForcesSingleThread();
#else
	totalPotentialEnergy = Forces();
#endif // SINGLE

	totalKineticEnergy = 0.0L;

	// Second part of Verlet
	for (i = 0; i < NUM_PARTS; i++) {

		_dt_2m = DELTA_T / (2.0L * myParticles[i].getMass());

		myParticles[i].setVelocity(
			myParticles[i].getVx() + fx_pt[i] * _dt_2m,
			myParticles[i].getVy() + fy_pt[i] * _dt_2m,
			myParticles[i].getVz() + fz_pt[i] * _dt_2m
		);

		totalKineticEnergy += myParticles[i].getMass() * ( myParticles[i].getVx() * myParticles[i].getVx() 
			+ myParticles[i].getVy() * myParticles[i].getVy()
			+ myParticles[i].getVz() * myParticles[i].getVz() ) / 2.0L;
	}

	totalEnergy_ = totalKineticEnergy + totalPotentialEnergy;

	tmp_en = fabsl((totalEnergy_ / initialEnergy_ - 1.0L) * 100.0L);
	if (tmp_en > relativeEnergyError_) relativeEnergyError_ = tmp_en;

	stepNo++;
	setDistribV();
}

void GasBox::switchVelocities() {
	for (int i = 0; i < NUM_PARTS; i++)
	{
		myParticles[i].setVelocity(-myParticles[i].getVx(), -myParticles[i].getVy(), -myParticles[i].getVz());
	}
}

unsigned int GasBox::getStep() {
	return stepNo;
}

long double GasBox::getTotalEnergy() {
	return totalEnergy_;
}

long double GasBox::getEnergyError() {
	return relativeEnergyError_;
}

int GasBox::getNumberOfParticles() {
	return NUM_PARTS;
}

std::vector<double> GasBox::getParticleLocation(int i) {
	std::vector<double> Point = { 0.0L, 0.0L, 0.0L };
	if (i < 0 || i >= NUM_PARTS) return Point;
	Point[0] = (double)myParticles[i].getX();
	Point[1] = (double)myParticles[i].getY();
	Point[2] = (double)myParticles[i].getZ();

	return Point;
}

std::vector<double> GasBox::getParticleVelocity(int i) {
	std::vector<double> Velocity = { 0.0L, 0.0L, 0.0L };
	if (i < 0 || i >= NUM_PARTS) return Velocity;
	Velocity[0] = (double)myParticles[i].getVx();
	Velocity[1] = (double)myParticles[i].getVy();
	Velocity[2] = (double)myParticles[i].getVz();

	return Velocity;
}

std::vector<double> GasBox::getParticleForces(int i) {
	std::vector<double> Forces = { 0.0L, 0.0L, 0.0L };
	if (i < 0 || i >= NUM_PARTS) return Forces;
	Forces[0] = (double)fx_pt[i];
	Forces[1] = (double)fy_pt[i];
	Forces[2] = (double)fz_pt[i];

	return Forces;
}

void GasBox::initMaxwell()
{
	double a, b, yval, vdiff, v, width;

	b = 3.0f / (2.0f * V_START * V_START);
	a = 100.0f * b * exp(1.0f);

	vdiff = V_MAX / (double)NUM_BOXES;
	// just use BOX_X as a measure of width for the distribution drawing area
	width = BOX_X / (double)NUM_BOXES;

	for (int i = 0; i < NUM_BOXES; i++)
	{
		v = ((double)i + 0.5L) * vdiff;

		yval = a * v * v * exp(-b * v * v);

		// use BOX_Y as a measure of height for the distribution drawing area
		MaxwPoints[i][0] = (width * (double)i);
		MaxwPoints[i][1] = (int)(BOX_Y + 150.0 - yval);
	}

	BezierSpline::GetCurveControlPoints(MaxwPoints, NUM_BOXES, firstControlPoints, secondControlPoints);
}

double GasBox::getMaxwell(int point, int coordinate) {
	return MaxwPoints[point][coordinate];
}

double GasBox::getFirstCP(int point, int coordinate) {
	return firstControlPoints[point][coordinate];
}

double GasBox::getSecondCP(int point, int coordinate) {
	return secondControlPoints[point][coordinate];
}

/*
* Calculate the velocity distribution of all particles into NUM_BOXES velocity sections
* and ignoring all particles faster then V_MAX
* The number of particles faster then V_MAX will be returned.
*/
int GasBox::setDistribV()
{
	int i, nmax, pos, nparts;
	double v, vdiff;

	vdiff = V_MAX / (double)NUM_BOXES;

	for (i = 0; i < NUM_BOXES; i++)
		calcDistribution[i] = 0.0;

	for (i = 0; i < NUM_PARTS; i++)
	{
		v = (double)(myParticles[i].getVx() * myParticles[i].getVx() 
			+ myParticles[i].getVy() * myParticles[i].getVy() 
			+ myParticles[i].getVz() * myParticles[i].getVz());
		v = sqrt(v);
		pos = (int)floor(v / vdiff);
		if (pos < NUM_BOXES)
			calcDistribution[pos]+=1.0;
	}

	nmax = 0;
	nparts = 0;
	for (i = 0; i < NUM_BOXES; i++)
	{
		if (calcDistribution[i] > nmax) nmax = calcDistribution[i];
		nparts += calcDistribution[i];
	}

	nparts = NUM_PARTS - nparts;

	// this should not happen, that all parts are out of the distribution
	if (nmax == 0) return nparts;

	for (i = 0; i < NUM_BOXES; i++)
	{
		calcDistribution[i] = (int)(calcDistribution[i] * 100 / nmax);
	}

	return nparts;
}

int GasBox::getDistribV(int n)
{
	if (n < 0 || n >= NUM_BOXES) return 0;
	return calcDistribution[n];
}


