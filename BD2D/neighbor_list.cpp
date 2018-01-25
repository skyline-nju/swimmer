#include "neighbor_list.h"


NeighborList_2::NeighborList_2(double Lx0, double Ly0, double r_cut0):
                               Lx(Lx0), Ly(Ly0), r_cut(r_cut0), mat(nullptr){
  half_Lx = 0.5 * Lx;
  half_Ly = 0.5 * Ly;
  cl = new CellList_w_list(Lx, Ly, r_cut0);
}

NeighborList_2::NeighborList_2(double Lx0, double Ly0, double r_cut0,
                               double r_buf_ratio, int nPar): 
                               Lx(Lx0), Ly(Ly), r_cut(r_cut0), nrows(nPar) {
  half_Lx = 0.5 * Lx;
  half_Ly = 0.5 * Ly;
  cl = new CellList_w_list(Lx, Ly, r_cut0);
  mat = nullptr;
  set_r_buf(r_buf_ratio, nPar);
}

NeighborList_2::~NeighborList_2() {
  delete[] mat;
  delete cl;
  double mean = sum_interval / count_refresh;
  double std = sqrt(sum_interval_square / count_refresh - mean * mean);
  std::cout << "neighbor list had been refreshed for " << count_refresh << " times\n";
  std::cout << "mean interval: " << mean << "\n";
  std::cout << "std: " << std << "\n";
  std::cout << "min interval: " << min_interval << "\n";
  std::cout << "max interval: " << max_interval << "\n";
  std::cout << "max neighbor: " << max_neighbor << "\n";
}

void NeighborList_2::set_r_buf(double r_buf_ratio, int nPar) {
  r_buf = r_buf_ratio * r_cut;
  half_r_buf_square = 0.25 * r_buf * r_buf;
  r_verlet = r_cut + r_buf;
  r_verlet_square = r_verlet * r_verlet;

  double rmin = 0.5;
  int n1 = int(r_verlet / rmin) + 2;
  ncols = n1 * n1;        // max neighbors of one particle
  std::cout << "Max length of one neighbor list is " << ncols << std::endl;

  nrows = nPar;
  ntot = ncols * nrows;
  if (mat)
    delete[] mat;
  mat = new int[ntot];

  std::cout << "nrows = " << nrows << "\n";
  std::cout << "ncols = " << ncols << "\n";
  std::cout << "r_cut = " << r_cut << "\n";
  std::cout << "r_verlet = " << r_verlet << "\n";
  std::cout << "r_buf = " << r_buf << std::endl;

  last_update_time = 0;
  count_refresh = 0;
  sum_interval = 0;
  sum_interval_square = 0;
  min_interval = 0;
  max_interval = 0;

  max_neighbor = 0;
}

