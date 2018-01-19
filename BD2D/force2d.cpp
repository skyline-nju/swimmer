#include "force2d.h"
#include <iostream>
#include <cmath>

F_base_2::F_base_2(double Lx0, double Ly0) : Lx(Lx0), Ly(Ly0) {
  half_Lx = 0.5 * Lx;
  half_Ly = 0.5 * Ly;
}

F_WCA_2::F_WCA_2(double Lx0, double Ly0, double eps) : F_base_2(Lx0, Ly0) {
  eps24 = eps * 24;
  r_cut = std::pow(2.0, 1.0 / 6.0);
  r_cut_square = r_cut * r_cut;
  std::cout << "the cutoff distance for WCA force is " << r_cut << std::endl;
}
