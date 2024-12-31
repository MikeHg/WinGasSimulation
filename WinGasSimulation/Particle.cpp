#include "pch.h"
#include "Particle.h"
#include <math.h>

Particle::Particle() {
	x_ = 0.0L; y_ = 0.0L;	z_ = 0.0L;
	vx_ = 0.0L; vy_ = 0.0L;	vz_ = 0.0L;
	mass_ = 0.1L;
}

void Particle::setLocation(long double x, long double y, long double z) {
	x_ = x; y_ = y;	z_ = z;
}

void Particle::setVelocity(long double vx, long double vy, long double vz) {
	vx_ = vx; vy_ = vy;	vz_ = vz;
}

long double Particle::getX() { return x_; }
long double Particle::getY() { return y_; }
long double Particle::getZ() { return z_; }
long double Particle::getVx() { return vx_; }
long double Particle::getVy() { return vy_; }
long double Particle::getVz() { return vz_; }

void Particle::setMass(long double mass) { mass_ = mass; }
long double Particle::getMass() { return mass_; }

long double Particle::getKineticEnergy() { 
	return mass_ * (vx_ * vx_ + vy_ * vy_ + vz_ * vz_) / 2.0L; 
}

long double Particle::getDistance(long double xref, long double yref, long double zref) {
	long double r2;
	r2 = (xref - x_) * (xref - x_) +
		(yref - y_) * (yref - y_) +
		(zref - z_) * (zref - z_);
	return sqrtl(r2);
}


// force = - dV / dr
// force[r] = - 4*epsilon*((-12*sigma^12)/r^13 + (6*sigma^6)/r^7)
// keep in Mind:
// dV / dx = dV / dr * dr / dx
// for r = sqrt( Delta_x^2 + Delta_y^2 +Delta_z^2) -> dr/dx = Delta_x / r, i.e.
// force_x = +6 (4 * epsilon)* (sigma^6)/(r^6)*( 2 * sigma^6 / r^6 - 1 ) * Delta_x/r^2 

// calculate force on particle ref in potential of this particle
void Particle::getForce(long double xref, long double yref, long double zref,
	long double* fx, long double* fy, long double* fz) 
{
	long double r2, r6, r2t, _tmp;

	r2 = (x_ - xref) * (x_ - xref)
		+ (y_ - yref) * (y_ - yref)
		+ (z_ - zref) * (z_ - zref);

	r2t = 1.0 / r2;
	r2 = (SIGMA * SIGMA) * r2t;

	r6 = r2 * r2 * r2;

	_tmp = 6.0L * E4 * r6 * (2.0L * r6 - 1.0L) * r2t;

	*fx = _tmp * (x_ - xref);
	*fy = _tmp * (y_ - yref);
	*fz = _tmp * (z_ - zref);
}

// potential[r] = 4 epsilon ((sigma/r)^12 - (sigma/r)^6)
// minimum at r = 2^(1/6) * sigma = 1.122462 * sigma

// returns potential energy (lennard jones potential)
// of particle ref in potential of this particle
long double Particle::getPotentialEnergy(long double xref, long double yref, long double zref)
{
	long double r2, r6, pot;

	r2 = (x_ - xref) * (x_ - xref)
		+ (y_ - yref) * (y_ - yref)
		+ (z_ - zref) * (z_ - zref);

	r2 = (SIGMA * SIGMA) / r2;

	r6 = r2 * r2 * r2;

	pot = E4 * r6 * (r6 - 1.0L);

	return pot;
}

long double Particle::getForceAndPotEnergy(long double xref, long double yref, long double zref,
	long double* fx, long double* fy, long double* fz)
{
	long double r2, r6, r2t, _tmp, pot;

	r2 = (x_ - xref) * (x_ - xref)
		+ (y_ - yref) * (y_ - yref)
		+ (z_ - zref) * (z_ - zref);

	r2t = 1.0 / r2;
	r2 = (SIGMA * SIGMA) * r2t;

	r6 = r2 * r2 * r2;

	_tmp = 6.0L * E4 * r6 * (2.0L * r6 - 1.0L) * r2t;

	*fx = _tmp * (x_ - xref);
	*fy = _tmp * (y_ - yref);
	*fz = _tmp * (z_ - zref);

	pot = E4 * r6 * (r6 - 1.0L);
	return pot;
}

