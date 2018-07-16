#include "benchmark.h"
#include <cmath>
#include <ctime>
#include <iostream>

void lbm::D2Q9::taylor_green(unsigned int t, unsigned int x, unsigned int y,
                             int NX, int NY, double rho0, double u_max, double nu,
                             double * r, double * u, double * v) {
  const double PI = 3.14159265358979323846;
  double kx = 2.0 * PI / NX;
  double ky = 2.0 * PI / NY;
  double td = 1.0 / (nu * (kx * kx + ky * ky));
  double X = x + 0.5;
  double Y = y + 0.5;
  double ux = -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y) * exp(-1.0 * t / td);
  double uy = u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y) * exp(-1.0 * t / td);
  double P = -0.25 * rho0 * u_max * u_max * (
    (ky / kx) * cos(2.0 * kx * X) + (kx / ky) * cos(2.0 * ky * Y)
    ) * exp(-2.0 * t / td);
  double rho = rho0 + 3.0 * P;

  *r = rho;
  *u = ux;
  *v = uy;
}

void lbm::D2Q9::taylor_green(unsigned int t, double rho0, double u_max, double nu,
                             Lattice &lat2d) {
  for (int y = 0; y < lat2d.ny; y++) {
    for (int x = 0; x < lat2d.nx; x++) {
      size_t sidx = lat2d.scalar_idx(x, y);
      taylor_green(t, x, y, lat2d.nx, lat2d.ny, rho0, u_max, nu,
                   &lat2d.rho[sidx], &lat2d.ux[sidx], &lat2d.uy[sidx]);
    }
  }
  lat2d.init_pdf();
}

void lbm::D2Q9::compute_flow_properties(unsigned int t, Lattice & lat, double * prop,
                                        double rho0, double u_max, double nu) {
  // prop must point to space for 4 doubles:
  // 0: energy
  // 1: L2 error in rho
  // 2: L2 error in ux
  // 3: L2 error in uy
  double E = 0.0;
  double sumrhoe2 = 0.0;
  double sumuxe2 = 0.0;
  double sumuye2 = 0.0;
  double sumrhoa2 = 0.0;
  double sumuxa2 = 0.0;
  double sumuya2 = 0.0;
  for (int y = 0; y < lat.ny; ++y) {
    for (int x = 0; x < lat.nx; ++x) {
      double rho = lat.get_rho(x, y);
      double ux = lat.get_ux(x, y);
      double uy = lat.get_uy(x, y);
      E += rho * (ux*ux + uy * uy);
      double rhoa, uxa, uya;
      taylor_green(t, x, y, lat.nx, lat.ny, rho0, u_max, nu, &rhoa, &uxa, &uya);
      sumrhoe2 += (rho - rhoa)*(rho - rhoa);
      sumuxe2 += (ux - uxa)*(ux - uxa);
      sumuye2 += (uy - uya)*(uy - uya);
      sumrhoa2 += (rhoa - rho0)*(rhoa - rho0);
      sumuxa2 += uxa * uxa;
      sumuya2 += uya * uya;
    }
  }
  prop[0] = E;
  prop[1] = sqrt(sumrhoe2 / sumrhoa2);
  prop[2] = sqrt(sumuxe2 / sumuxa2);
  prop[3] = sqrt(sumuye2 / sumuya2);
}

void lbm::D2Q9::benchmark_taylor_green(int nx, int ny) {
  const int scale = 1;
  const double nu = 1 / 6.0;
  const double tau = 3.0 * nu + 0.5;
  const double u_max = 0.04 / scale;
  const double rho0 = 1.0;
  const int n_step = 2000 * scale * scale;

  BGK bgk(tau);
  Lattice lattice(nx, ny);
  taylor_green(0, rho0, u_max, nu, lattice);

  clock_t start = clock();

  for (int i = 0; i < n_step; i++) {
    bool save = i % 500 == 0;
    lattice.stream_collide(bgk, save);
    if (save) {
      double prop[4];
      compute_flow_properties(i, lattice, prop, rho0, u_max, nu);
      std::cout << "t = " << i << "\n";
      std::cout << "energy = " << prop[0] << "\n";
      std::cout << "L2 error in rho = " << prop[1] << "\n";
      std::cout << "L2 error in ux = " << prop[2] << "\n";
      std::cout << "L2 error in uy = " << prop[3] << std::endl;
    }
  }
  double runtime = double(clock() - start) / CLOCKS_PER_SEC;
  size_t nodes_updated = n_step * size_t(nx * ny);
  double speed = nodes_updated / (1e6 * runtime);
  double bytesPerGiB = 1024. * 1024. * 1024.;
  double bandwidth = nodes_updated * 2 * ndir * sizeof(double) / (runtime * bytesPerGiB);
  std::cout << "speed = " << speed << std::endl;
  std::cout << "bandwidth = " << bandwidth << std::endl;
}
