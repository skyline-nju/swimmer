#ifndef BROWNIAN_H
#define BROWNIAN_H
#include "cmdline.h"
#include "force.h"
#include "vect.h"
#include "boundary.h"
#include "rand.h"
#include "cellList.h"
#include "ioData.h"
#include "integrate.h"
#include "particle2d.h"
#include "neighborList.h"

class DynamicBase_2 {
public:
  DynamicBase_2(const cmdline::parser &cmd);
  virtual ~DynamicBase_2();

  virtual void run(int steps)=0;

  void run();

  template <class Par>
  void ini_pos_rand(std::vector<Par> &p_arr, int np, double sigma = 1) const;

  template <class Par>
  void output(int i_step, std::vector<Par> &p_arr);

  void set_spatial_sort(bool flag) { USE_SPATIAL_SORT = flag; }


protected:
  double Lx;
  double Ly;
  int nPar;
  int nstep;
  double packing_fraction;
  double h;
  bool OUTPUT_ON;
  bool USE_SPATIAL_SORT;

  Ran *myran;

  XY_Writer *xy_out;
  LogWriter *log_out;

  PBC_2 pbc2;

  std::vector<std::ofstream> fout;
};

template<class Par>
void DynamicBase_2::ini_pos_rand(std::vector<Par>& p_arr, int np, double sigma) const {
  if (packing_fraction < 0.5) {
    create_rand_2(p_arr, np, sigma, Lx, Ly, myran);
  } else {
    double sigma0 = sigma;
    sigma = 0.5;
    create_rand_2(p_arr, np, sigma, Lx, Ly, myran);
    
    std::vector< Vec_2<double> > f_arr;
    f_arr.reserve(np);
    for (int i = 0; i < np; i++) {
      f_arr.emplace_back();
    }

    EulerMethod euler(h);
    WCAForce f_WCA(1, 1);
    auto f_ij = [&p_arr, &f_arr, this, &f_WCA](auto i, auto j) {
      Vec_2<double> distance(p_arr[i].x - p_arr[j].x, p_arr[i].y - p_arr[j].y);
      pbc2.nearest_dis(distance);
      f_WCA(f_arr[i], f_arr[j], distance);
    };
    double d = 0.01;
    int n = int((sigma0 - sigma) / d);
    int delta_t = int(d / h);
    std::cout << "The diameter is increased by " << d << " every "
              << delta_t << " time steps" << std::endl;

    CellList_list_2 cl(Lx, Ly, f_WCA.get_r_cut(), 0);
    cl.create(p_arr);
    for (int i = 0; i < n; i++) {
      sigma += d;
      f_WCA.set_sigma(sigma);
      for (int t = 0; t < delta_t; t++) {
        cl.for_each_pair(f_ij);
        for (int ip = 0; ip < np; ip++) {
          euler.update_xy(p_arr[ip], f_arr[ip], pbc2, myran);
        }
        cl.update(p_arr);
      }
    }
    std::cout << "The diameter is increased to " << sigma << std::endl;
  }
}

template<class Par>
inline void DynamicBase_2::output(int i_step, std::vector<Par>& p_arr) {
  if (OUTPUT_ON) {
    (*xy_out)(i_step, p_arr);
    (*log_out)(i_step);
  }
}
template <class Par, class TList>
class BD_2 : public DynamicBase_2 {
public:
  BD_2(const cmdline::parser &cmd);
  ~BD_2() { delete tlist; }
  void run(int nsteps);
protected:
  std::vector<Par> p_arr;
  std::vector<Vec_2<double>> f_arr;
  TList *tlist;
};

template<class Par, class TList>
BD_2<Par, TList>::BD_2(const cmdline::parser & cmd):DynamicBase_2(cmd) {
  p_arr.reserve(nPar);
  ini_pos_rand(p_arr, nPar);
  f_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    f_arr.emplace_back();
  }
  WCAForce fwca(1, 1);
  tlist = new TList(Lx, Ly, fwca.get_r_cut(), 0.4, nPar);
  tlist->create(p_arr);
  if (OUTPUT_ON) {
    (*xy_out)(0, p_arr);
  }
}

template<class Par, class TList>
void BD_2<Par, TList>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  MySpatialSortingTraits<Par> sst;
  auto lambda = [&](int i, int j) {
    Vec_2<double> dis(p_arr[i].x - p_arr[j].x, p_arr[i].y - p_arr[j].y);
    pbc2.nearest_dis(dis);
    fwca(f_arr[i], f_arr[j], dis);
  };
  for (int i = 1; i <= nsteps; i++) {
    if (USE_SPATIAL_SORT)
      tlist->cal_force(p_arr, lambda, sst);
    else
      tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < nPar; ip++) {
      euler.update_xy(p_arr[ip], f_arr[ip], pbc2, myran);
    }
    output(i, p_arr);
  }
}

template <class Par, class TList>
class ABD_2 : public DynamicBase_2 {
public:
  ABD_2(const cmdline::parser &cmd);
  ~ABD_2() { delete tlist; }
  void run(int nsteps);

protected:
  std::vector<Par> p_arr;
  std::vector<Vec_2<double>> f_arr;
  TList *tlist;
  double Pe;
};

template<class Par, class TList>
ABD_2<Par, TList>::ABD_2(const cmdline::parser & cmd): DynamicBase_2(cmd) {
  Pe = cmd.get<double>("Pe");
  p_arr.reserve(nPar);
  ini_pos_rand(p_arr, nPar);
  f_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    p_arr[i].theta = myran->doub() * PI * 2.;
    f_arr.emplace_back();
  }
  WCAForce fwca(1, 1);
  tlist = new TList(Lx, Ly, fwca.get_r_cut(), 0.35, nPar);
  tlist->create(p_arr);
  if (OUTPUT_ON) {
    (*xy_out)(0, p_arr);
  }

}

template<class Par, class TList>
void ABD_2<Par, TList>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  MySpatialSortingTraits<Par> sst;
  auto lambda = [&](int i, int j) {
    Vec_2<double> dis(p_arr[i].x - p_arr[j].x, p_arr[i].y - p_arr[j].y);
    pbc2.nearest_dis(dis);
    fwca(f_arr[i], f_arr[j], dis);
  };
  for (int i = 1; i <= nsteps; i++) {
    if (USE_SPATIAL_SORT)
      tlist->cal_force(p_arr, lambda, sst);
    else
      tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < nPar; ip++) {
      euler.update_xy_theta(p_arr[ip], f_arr[ip], Pe, pbc2, myran);
    }
    output(i, p_arr);
  }

}

#endif
