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
