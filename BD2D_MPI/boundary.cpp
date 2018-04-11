#include "boundary.h"

BoundaryBase_2::BoundaryBase_2(const cmdline::parser& cmd) {
  lx_ = cmd.get<double>("Lx");
  ly_ = cmd.exist("Ly") ? cmd.get<double>("Ly") : lx_;
  x_min_ = 0;
  y_min_ = 0;
}

void BoundaryBase_2::set_para0(double lx, double ly, double x0, double y0) {
  lx_ = lx;
  ly_ = ly;
  x_min_ = x0;
  y_min_ = y0;
}

PBC_x_2::PBC_x_2(const cmdline::parser & cmd): BoundaryBase_2(cmd) {
  half_lx_ = 0.5 * lx_;
  x_max_ = x_min_ + lx_;
}

PBC_x_2::PBC_x_2(double lx, double ly, double x0, double y0) :
  BoundaryBase_2(lx, ly, x0, y0) {
  half_lx_ = 0.5 * lx_;
  x_max_ = x_min_ + lx_;
}

void PBC_x_2::set_para(double lx, double ly, double x0, double y0) {
  set_para0(lx, ly, x0, y0);
  half_lx_ = 0.5 * lx_;
  x_max_ = x_min_ + lx_;
}

void PBC_x_2::new_rand_pos(double & x, double & y, Ran &myran) const {
  x = x_min_ + myran.doub() * lx_;
  y = y_min_ + myran.doub() * ly_;
}

PBC_xy_2::PBC_xy_2(const cmdline::parser & cmd): BoundaryBase_2(cmd) {
  half_lx_ = 0.5 * lx_;
  half_ly_ = 0.5 * ly_;
  x_max_ = x_min_ + lx_;
  y_max_ = y_min_ + ly_;
}

PBC_xy_2::PBC_xy_2(double lx, double ly, double x0, double y0) :
  BoundaryBase_2(lx, ly, x0, y0) {
  half_lx_ = 0.5 * lx_;
  half_ly_ = 0.5 * ly_;
  x_max_ = x_min_ + lx_;
  y_max_ = y_min_ + ly_;
}

void PBC_xy_2::set_para(double lx, double ly, double x0, double y0) {
  set_para0(lx, ly, x0, y0);
  half_lx_ = 0.5 * lx_;
  half_ly_ = 0.5 * ly_;
  x_max_ = x_min_ + lx_;
  y_max_ = y_min_ + ly_;
}

void PBC_xy_2::new_rand_pos(double & x, double & y, Ran &myran) const {
  x = x_min_ + myran.doub() * lx_;
  y = y_min_ + myran.doub() * ly_;
}

Wall_x_2::Wall_x_2(const cmdline::parser & cmd): PBC_x_2(cmd) {
  k_ = cmd.get<double>("k_wall");
  const auto a = cmd.get<double>("sigma") * 0.5;
  xl_ = x_min_ + a;
  xr_ = x_max_ - a;
}

void Wall_x_2::new_rand_pos(double & x, double & y, Ran & myran) const {
  const auto lx = xr_ - xl_;
  x = xl_ + myran.doub() * lx;
  y = y_min_ + myran.doub() * ly_;
}

Wall_x_PBC_y_2::Wall_x_PBC_y_2(const cmdline::parser & cmd)
  : Wall_x_2(cmd) {
  half_ly_ = 0.5 * ly_;
  y_max_ = y_min_ + ly_;
}
