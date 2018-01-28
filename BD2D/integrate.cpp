#include "integrate.h"
#include <cmath>
#include <iostream>

EulerMethod::EulerMethod(double h0): h(h0) {
  sqrt_24h = std::sqrt(24 * h);
  sqrt_72h = std::sqrt(72 * h);
  trible_h = h * 3;
  std::cout << "h = " << h << std::endl;
}
