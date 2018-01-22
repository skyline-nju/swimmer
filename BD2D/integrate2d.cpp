#include "integrate2d.h"

Int_EM_2::Int_EM_2(double H, double LX, double LY) : Int_base_2(H, LX, LY) {
  sqrt_24h = std::sqrt(24 * h);
  sqrt_72h = std::sqrt(72 * h);
  trible_h = h * 3;
  std::cout << "h = " << h << "\tLx = " << Lx << "\tLy = " 
            << Ly << std::endl;
}