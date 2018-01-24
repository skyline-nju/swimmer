#include "force2d.h"
#include <iostream>
#include <cmath>

F_base_2::F_base_2(double Lx0, double Ly0) : Lx(Lx0), Ly(Ly0) {
  half_Lx = 0.5 * Lx;
  half_Ly = 0.5 * Ly;
}

F_WCA_2::F_WCA_2(double Lx0, double Ly0, double eps) : F_base_2(Lx0, Ly0) {
  eps24 = eps * 24;
  r_cut = std::pow(2.0, 1. / 6.);
  r_cut_square = r_cut * r_cut;
  set_sigma(1);
  std::cout << "the cutoff distance for WCA force is " << r_cut << std::endl;
}

void::F_WCA_2::set_sigma(double sigma0) {
  sigma = sigma0;
  sigma_inverse_square = 1 / (sigma * sigma);
}