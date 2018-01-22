#include "particle2d.h"
#include "comn.h"

bool overlap_2(double dx, double dy, double sigma, double Lx, double Ly) {
  double half_Lx = 0.5 * Lx;
  double half_Ly = 0.5 * Ly;

  if (dx < -half_Lx) {
    dx += Lx;
  } else if (dx > half_Lx) {
    dx -= Lx;
  }
  if (dy < -half_Ly) {
    dy += Ly;
  } else if (dy > half_Ly) {
    dy -= Ly;
  }

  if (dx * dx + dy * dy < sigma * sigma) {
    return true;
  } else {
    return false;
  }
}

void BP_2::set_data(double x0, double y0, Ran * myran) {
  x = x0;
  y = y0;
}

void BP_with_ori_2::set_data(double x0, double y0, Ran * myran) {
  x = x0;
  y = y0;
  theta = myran->doub() * 2 * PI;
  std::cout << "theta = " << theta << std::endl;
}
