#ifndef MIXTURE_H
#define MIXTURE_H
#include "brownian.h"

template <class Par, class Tlist>
class CW_CCW_AB_2 : public DynamicBase_2 {
public:
  CW_CCW_AB_2(const cmdline::parser &cmd);
  ~CW_CCW_AB_2() { delete tlist; }
  void run(int nsteps);
protected:
  std::vector<Par> p_arr;
  std::vector<Vec_2<double>>f_arr;
  Tlist *tlist;
  double Pe;
  double tau;
};

template<class Par, class Tlist>
CW_CCW_AB_2<Par, Tlist>::CW_CCW_AB_2(const cmdline::parser & cmd):
  DynamicBase_2(cmd) {
  Pe = cmd.get<double>("Pe");
  tau = cmd.get<double>("tau");
  p_arr.reserve(nPar);
  ini_pos_rand(p_arr, nPar);
  f_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    p_arr[i].theta = myran->doub() * PI * 2.;
    f_arr.emplace_back();
  }
  WCAForce fwca(1, 1);
  tlist = new Tlist(Lx, Ly, fwca.get_r_cut(), 0.35, nPar);
  tlist->create(p_arr);
  if (output_on) {
    xy_out->write_mix(0, p_arr, nPar / 2);
  }
}

template<class Par, class Tlist>
void CW_CCW_AB_2<Par, Tlist>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  MySpatialSortingTraits<Par> sst;
  int half_n = nPar / 2;
  auto lambda = [&](int i, int j) {
    Vec_2<double> dis(p_arr[i].x - p_arr[j].x, p_arr[i].y - p_arr[j].y);
    pbc2.nearest_dis(dis);
    fwca(f_arr[i], f_arr[j], dis);
  };
  for (int i = 1; i <= nsteps; i++) {
    if (spatial_sort_on)
      tlist->cal_force(p_arr, lambda, sst);
    else
      tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < half_n; ip++) {
      euler.update_xy_theta(p_arr[ip], f_arr[ip], Pe, tau, pbc2, myran);
    }
    for (int ip = half_n; ip < nPar; ip++) {
      euler.update_xy_theta(p_arr[ip], f_arr[ip], Pe, -tau, pbc2, myran);
      //euler.update_xy(p_arr[ip], f_arr[ip], pbc2, myran);
    }
    if (output_on) {
      xy_out->write_mix(i, p_arr, half_n);
      (*log_out)(i);
    }
  }
}

template <class Par, class Tlist>
class CW_CCW_AB_Dipole_2 : public DynamicBase_2 {
public:
  CW_CCW_AB_Dipole_2(const cmdline::parser &cmd);
  ~CW_CCW_AB_Dipole_2() { delete tlist; }
  void run(int nsteps);
protected:
  std::vector<Par> p_arr;
  std::vector<Vec_3<double>> f_arr;
  Tlist *tlist;
  double r_cut;
  double epsilon;
  double ratio;
  double tau;
  double Pe;
};

template<class Par, class Tlist>
CW_CCW_AB_Dipole_2<Par, Tlist>::CW_CCW_AB_Dipole_2(const cmdline::parser & cmd):
DynamicBase_2(cmd) {
  p_arr.reserve(nPar);
  ini_pos_rand(p_arr, nPar);
  f_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    f_arr.emplace_back();
    double theta = myran->doub() * 2 * PI;
    p_arr[i].u.x = std::cos(theta);
    p_arr[i].u.y = std::sin(theta);
  }
  r_cut = 3;
  epsilon = cmd.get<double>("dipole_strength");
  ratio = cmd.get<double>("dipole_ratio");
  tau = cmd.get<double>("tau");
  Pe = cmd.get<double>("Pe");

  tlist = new Tlist(Lx, Ly, r_cut, 0.3, nPar);
  tlist->create(p_arr);
  if (output_on) {
    xy_out->write_theta(0, p_arr, nPar / 2);
  }
}

template<class Par, class Tlist>
void CW_CCW_AB_Dipole_2<Par, Tlist>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  double qh, qt;
  if (ratio < 0) {
    qh = std::sqrt(-ratio);
    qt = -1 / qh;
  } else {
    qh = std::sqrt(ratio);
    qt = 1 / qh;
  }
  ExtDipoleForce fed(epsilon, qh, qt, r_cut, 4, 3.0 / 16);
  auto lambda = [&](int i, int j) {
    fed(f_arr[i], f_arr[j], p_arr[i], p_arr[j], pbc2, fwca);
  };
  int half_nPar = nPar / 2;
  for (int i = 1; i <= nsteps; i++) {
    tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < half_nPar; ip++) {
      euler.update_xy_uxuy(p_arr[ip], f_arr[ip], Pe, tau, pbc2, myran);
    }
    for (int ip = half_nPar; ip < nPar; ip++) {
      euler.update_xy_uxuy(p_arr[ip], f_arr[ip], Pe, -tau, pbc2, myran);
    }
    if (output_on) {
      xy_out->write_theta(i, p_arr, half_nPar);
      (*log_out)(i);
    }
  }
}

#endif