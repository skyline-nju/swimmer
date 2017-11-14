#ifndef LBM_2D_H
#define LBM_2D_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

const double PI = 3.14159265358979323846;

	// domain size
	const unsigned int scale = 1;
	const unsigned int NX = 1024 * scale;
	const unsigned int NY = NX;

	// the number of directions in the lattice
	const unsigned int ndir = 9;

	// the memory size for the populations and scalar values
	const size_t size_ndir = NX * NY * ndir;
	const size_t size_scalar = NX * NY;

	// The lattice weights
	const double w0 = 4.0 / 9.0;   // zero weight
	const double ws = 1.0 / 9.0;   // adjacent weight
	const double wd = 1.0 / 36.0;  // diagonal weight

																 // arrays of the lattice weights and direction components
	const double wi[] = { w0, ws, ws, ws, ws, wd, wd, wd, wd };
	const int dirx[] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
	const int diry[] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

	// The kinematic viscosity nu and the corresponding relaxation parameter tau
	const double nu = 1.0 / 6.0;
	const double tau = 3.0 * nu + 0.5;
	const double tauinv = 2.0 / (6.0 * nu + 1);
	const double omtauinv = 1.0 - tauinv;

	// The maximum flow speed
	const double u_max = 0.04 / scale;

	// The fluid density
	const double rho0 = 1.0;

	// The number of time steps in the simulation
	const unsigned int NSTEPS = 2000 * scale * scale;

	// functions for computing linear array indexes from two-dimensional coordinates
	inline size_t scalar_index(unsigned int x, unsigned int y) {
		return NX * y + x;
	}
	inline size_t field_index(unsigned int x, unsigned int y, unsigned int d) {
		return NX * (NY * d + y) + x;
	}
	inline size_t field0_index(unsigned int x, unsigned int y) {
		return scalar_index(x, y);
	}
	inline size_t fieldn_index(unsigned int x, unsigned int y, unsigned int d) {
		return field_index(x, y, d - 1);
	}
	inline size_t fieldn_index2(unsigned int x, unsigned int y, unsigned int d) {
		return (ndir - 1) * (NX * y + x) + (d - 1);
	}

	void init_equilibrium(double *f, double *r, double *u, double *v);
	void init_equilibrium(double *f0, double *f1, double *r, double *u, double *v);
	void compute_rho_u(double *f, double *r, double *u, double *v);
	void stream(double *f, unsigned int x, unsigned int y);
	void collide(double *f0, double *f1, unsigned int x, unsigned int y, double *r, double *u, double *v, bool save);
	void stream_collide(double *f0, double *f1, double *r, double *u, double *v, bool save);
	void stream_collide(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save);
	void compute_flow_properties(unsigned int t, double *r, double *u, double *v, double *prop);

// functions used to compute the exact solution for Taylor-Green vortex decay
void taylor_green(unsigned int t, unsigned int x, unsigned int y,
	double *r, double *u, double *v);
void taylor_green(unsigned int t, double *r, double *u, double *v);
#endif
