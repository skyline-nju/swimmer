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



#endif
