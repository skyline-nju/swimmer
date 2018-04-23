#include "uniDomain.h"

UniDomain_2::UniDomain_2(const cmdline::parser & cmd)
  : flag_bc_(cmd.get<int>("bc_x"), cmd.get<int>("bc_y")),
    sigma_(cmd.get<double>("sigma")), myran_(cmd.get<ull>("seed")),
    log_(nullptr), xy_(nullptr), nc_(nullptr), profile_(nullptr) {
  l_.x = cmd.get<double>("Lx");
  l_.y = cmd.exist("Ly") ? cmd.get<double>("Ly") : l_.x;
  half_l_ = 0.5 * l_;
  const auto phi = cmd.get<double>("phi");
  n_par_ = cal_particle_number_2(phi, l_.x, l_.y, sigma_);

  // initial output
  if (cmd.exist("output")) {
    log_ = new LogExporter_2(cmd);
    const auto snap_format = cmd.get<std::string>("snap_fmt");
    if (snap_format == "xy") {
      xy_ = new XyExporter(cmd);
    } else if (snap_format == "nc") {
      nc_ = new NcParExporter_2(cmd);
    } else if (snap_format == "both") {
      xy_ = new XyExporter(cmd);
      nc_ = new NcParExporter_2(cmd);
    }
    if (cmd.exist("profile")) {
      profile_ = new ProfileExporter(cmd);
    }
  }
}


UniDomain_2::~UniDomain_2() {
  delete log_;
  delete xy_;
  delete nc_;
  delete profile_;
}

