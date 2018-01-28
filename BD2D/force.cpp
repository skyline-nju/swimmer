#include "force.h"

WCAForce::WCAForce(double eps, double sigma) {
  eps24 = eps * 24;
  double r_cut = std::pow(2.0, 1. / 6.);
  r_cut_square = r_cut * r_cut;
  set_sigma(sigma);
}

void WCAForce::set_sigma(double sigma) {
  sigma_inverse_square = 1 / (sigma * sigma);
}

