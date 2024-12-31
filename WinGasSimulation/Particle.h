#pragma once
#define E4	4.0L
#define SIGMA 5.0L
class Particle
{
private:
	long double x_, y_, z_;
	long double vx_, vy_, vz_;
	long double mass_;

public:
	Particle();

	void setLocation(long double x, long double y, long double z);
	void setVelocity(long double vx, long double vy, long double vz);
	long double getX();
	long double getY();
	long double getZ();
	long double getVx();
	long double getVy();
	long double getVz();
	void setMass(long double mass);
	long double getMass();

	long double getKineticEnergy();
	long double getDistance(long double xref, long double yref, long double zref);

	void getForce(long double xref, long double yref, long double zref,
		long double* fx, long double* fy, long double* fz);

	long double getPotentialEnergy(long double xref, long double yref, long double zref);

	long double getForceAndPotEnergy(long double xref, long double yref, long double zref,
		long double* fx, long double* fy, long double* fz);

};
