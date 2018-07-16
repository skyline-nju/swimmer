#pragma once

#include "D2Q9.h"

namespace lbm {

namespace D2Q9 {

void taylor_green(unsigned int t, unsigned int x, unsigned int y,
                  int NX, int NY, double rho0, double u_max, double nu,
                  double * r, double * u, double * v);

void taylor_green(unsigned int t, double rho0, double u_max, double nu,
                  Lattice &lat2d);

void compute_flow_properties(unsigned int t, Lattice &lat, double * prop,
                             double rho0, double u_max, double nu);

void benchmark_taylor_green(int nx, int ny);
} // namespace D2Q9

} // namespace lbm
