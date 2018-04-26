#include "uniDomain.h"

UniDomain_2::UniDomain_2(const cmdline::parser & cmd)
  : flag_bc_(cmd.get<int>("bc_x"), cmd.get<int>("bc_y")),
    sigma_(cmd.get<double>("sigma")), myran_(cmd.get<ull>("seed")),
    log_(nullptr), nc_(nullptr), profile_(nullptr) {
  l_.x = cmd.get<double>("Lx");
  l_.y = cmd.exist("Ly") ? cmd.get<double>("Ly") : l_.x;
  half_l_ = 0.5 * l_;
  const auto phi = cmd.get<double>("phi");
  n_par_ = cal_particle_number_2(phi, l_.x, l_.y, sigma_);
  std::cout << "particle number = " << n_par_ << std::endl;
  // initial output
  set_output_2(cmd, &log_, &nc_, &profile_);
}

UniDomain_2::~UniDomain_2() {
  delete log_;
  delete nc_;
  delete profile_;
}
