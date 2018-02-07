#include "neighborList.h"

NeighborList_2::NeighborList_2(double Lx, double Ly, double r_cut,
  double r_buf_ratio, int nPar): pbc_2(Lx, Ly) {
  cl = new CellList_list_2(Lx, Ly, r_cut, 0);

  double r_buf = r_cut * r_buf_ratio;
  half_r_buf_square = 0.25 * r_buf * r_buf;
  double r_verlet = r_cut + r_buf;
  r_verlet_square = r_verlet * r_verlet;

  double rmin = 0.5;
  int n1 = int(r_verlet / rmin) + 2;
  mat_bins.x = n1 * n1;        // max neighbors of one particle
  std::cout << "Max length of one neighbor list is "
            << mat_bins.x << std::endl;
  mat_bins.y = nPar;
  int ntot = mat_bins.x * mat_bins.y;
  mat.reserve(ntot);
  for (int i = 0; i < ntot; i++)
    mat.push_back(0);
  last_pos.reserve(nPar);

  refresh_interval = 0;
  refresh_count = 0;
  dt = 0;
  dt_square = 0;
  dt_min = 0;
  dt_max = 0;
  max_neighbor = 0;

  std::cout << "r_cut = " << r_cut << "\n";
  std::cout << "r_verlet = " << r_verlet << "\n";
  std::cout << "r_buf = " << r_buf << "\n";
  std::cout << "half of r_buf square = " << half_r_buf_square << "\n";
  std::cout << "matrix: " << mat_bins.x << " x " << mat_bins.y << "\n";
}

NeighborList_2::~NeighborList_2() {
  delete cl;
  double mean = dt / refresh_count;
  double std = sqrt(dt_square / refresh_count - mean * mean);
  std::cout << "neighbor list had been refreshed for " << refresh_count << " times\n";
  std::cout << "mean interval: " << mean << "\n";
  std::cout << "std: " << std << "\n";
  std::cout << "min interval: " << dt_min << "\n";
  std::cout << "max interval: " << dt_max << "\n";
  std::cout << "max neighbor: " << max_neighbor << "\n";
}

void NeighborList_2::count_time() {
  refresh_count++;
  dt += refresh_interval;
  dt_square += refresh_interval * refresh_interval;
  if (refresh_interval > dt_max)
    dt_max = refresh_interval;
  if (dt_min == 0 || dt_min > refresh_interval)
    dt_min = refresh_interval;

  double mean = dt / refresh_count;
  double err = std::sqrt(dt_square / refresh_count - mean * mean);
  safe_inteval = mean - 4 * err;
  refresh_interval = 0;
}

