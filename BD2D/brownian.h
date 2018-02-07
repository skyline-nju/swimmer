#ifndef BROWNIAN_H
#define BROWNIAN_H
#include <iomanip>
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
  void run_simple_bd(int n_steps, std::vector<Par> &p_arr,
                     std::vector< Vec_2<double> > &f_arr, EulerMethod &euler,
                     WCAForce &f_wca, CellList_list_2 &cl);

  template <class Par>
  void ini_pos_rand(std::vector<Par> &p_arr, int np, double sigma = 1.0);

  template <class Par>
  void output(int i_step, std::vector<Par> &p_arr);

  void set_spatial_sort(bool flag) { spatial_sort_on = flag; }

  template <class BiFunc>
  void for_each_pair(BiFunc f2) {
    for (int i = 0; i < nPar - 1; i++) {
      for (int j = i + 1; j < nPar; j++) {
        f2(i, j);
      }
    }
  }


protected:
  double Lx;
  double Ly;
  int nPar;
  int nstep;
  double packing_fraction;
  double h;
  bool output_on;
  bool spatial_sort_on;

  Ran *myran;

  XY_Writer *xy_out;
  LogWriter *log_out;

  PBC_2 pbc2;

  std::vector<std::ofstream> fout;
};

template <class Par>
void DynamicBase_2::run_simple_bd(int n_steps, std::vector<Par> &p_arr, 
                                  std::vector< Vec_2<double> > &f_arr,
                                  EulerMethod &euler, WCAForce &f_wca,
                                  CellList_list_2 &cl) {
  auto f_ij = [&p_arr, &f_arr, &f_wca, this](int i, int j) {
    Vec_2<double> dR(p_arr[i].x - p_arr[j].x, p_arr[i].y - p_arr[j].y);
    pbc2.nearest_dis(dR);
    f_wca(f_arr[i], f_arr[j], dR);
  };
  for (int i = 0; i < n_steps; i++) {
    cl.for_each_pair(f_ij);
    cl.for_each_pair(f_ij);
    for (int ip = 0; ip < nPar; ip++) {
      euler.update_xy(p_arr[ip], f_arr[ip], pbc2, myran);
    }
    cl.update(p_arr);
  }
}
template<class Par>
void DynamicBase_2::ini_pos_rand(std::vector<Par>& p_arr, int np, double sigma) {
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
    double d = 0.01;
    int n = int((sigma0 - sigma) / d);
    int delta_t = int(d / h);

    std::cout << "The diameter is increased by "  << d << " every "
              << delta_t << " time steps" << std::endl;

    CellList_list_2 cl(Lx, Ly, f_WCA.get_r_cut(), 0);
    cl.create(p_arr);
    for (int i = 0; i < n; i++) {
      sigma += d;
      f_WCA.set_sigma(sigma);
      run_simple_bd(delta_t, p_arr, f_arr, euler, f_WCA, cl);
    }
    f_WCA.set_sigma(sigma0);
    std::cout << "The diameter is increased to " << sigma
              << " after " << n << " steps\n";
  }
}

template<class Par>
inline void DynamicBase_2::output(int i_step, std::vector<Par>& p_arr) {
  if (output_on) {
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
  if (output_on) {
    (*xy_out)(0, p_arr);
  }
}

template<class Par, class TList>
void BD_2<Par, TList>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  MySpatialSortingTraits<Par> sst;
  auto lambda = [&](int i, int j) {
    fwca(f_arr[i], f_arr[j], p_arr[i], p_arr[j], pbc2);
  };
  for (int i = 1; i <= nsteps; i++) {
    if (spatial_sort_on)
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
  if (output_on) {
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
    if (spatial_sort_on)
      tlist->cal_force(p_arr, lambda, sst);
    else
      tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < nPar; ip++) {
      euler.update_xy_theta(p_arr[ip], f_arr[ip], Pe, pbc2, myran);
    }
    output(i, p_arr);
  }

}

template <class Par, class TList>
class BD_dipole_2 : public DynamicBase_2 {
public:
  BD_dipole_2(const cmdline::parser &cmd);
  ~BD_dipole_2() { delete tlist; }
  void run(int nsteps);
private:
  std::vector<Par> p_arr;
  std::vector<Vec_3<double>> f_arr;
  TList *tlist;
};

template<class Par, class TList>
BD_dipole_2<Par, TList>::BD_dipole_2(const cmdline::parser & cmd):
  DynamicBase_2(cmd) {
  p_arr.reserve(nPar);
  ini_pos_rand(p_arr, nPar);
  f_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    f_arr.emplace_back();
    p_arr[i].theta = myran->doub() * 2 * PI;
  }
  double r_cut = 4;
  tlist = new TList(Lx, Ly, r_cut, 0.4, nPar);
  tlist->create(p_arr);
  if (output_on) {
    xy_out->write_theta(0, p_arr);
  }
}

template<class Par, class TList>
void BD_dipole_2<Par, TList>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  //DipoleForce fd(4, 2.5);
  ExtDipoleForce fed(80, 1, -1, 4, 4, 3.0 / 16);
  MySpatialSortingTraits<Par> sst;
  auto lambda = [&](int i, int j) {
    fed(f_arr[i], f_arr[j], p_arr[i], p_arr[j], pbc2, fwca);
  };
  for (int i = 1; i <= nsteps; i++) {
    if (spatial_sort_on)
      tlist->cal_force(p_arr, lambda, sst);
    else
      tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < nPar; ip++) {
      euler.update_xy_theta(p_arr[ip], f_arr[ip], pbc2, myran);
    }
    if (output_on) {
      xy_out->write_theta(i, p_arr);
      (*log_out)(i);
    }
  }
}

template <class Par, class TList>
class CABD_dipole_2 : public DynamicBase_2 {
public:
  CABD_dipole_2(const cmdline::parser &cmd);
  ~CABD_dipole_2() { delete tlist; }
  void run(int nsteps);
private:
  std::vector<Par> p_arr;
  std::vector<Vec_3<double>> f_arr;
  TList *tlist;
  double r_cut;
  double epsilon;
  double ratio;
  double tau;
  double Pe;
};

template<class Par, class TList>
CABD_dipole_2<Par, TList>::CABD_dipole_2(const cmdline::parser & cmd):
  DynamicBase_2(cmd) {
  p_arr.reserve(nPar);
  ini_pos_rand(p_arr, nPar);
  f_arr.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    f_arr.emplace_back();
    //p_arr[i].theta = myran->doub() * 2 * PI;
    double theta = myran->doub() * 2 * PI;
    p_arr[i].u.x = std::cos(theta);
    p_arr[i].u.y = std::sin(theta);
  }
  r_cut = 3;
  epsilon = cmd.get<double>("dipole_strength");
  ratio = cmd.get<double>("dipole_ratio");
  tau = cmd.get<double>("tau");
  Pe = cmd.get<double>("Pe");

  tlist = new TList(Lx, Ly, r_cut, 0.3, nPar);
  tlist->create(p_arr);
  if (output_on) {
    xy_out->write_theta(0, p_arr);
  }
}

template<class Par, class TList>
void CABD_dipole_2<Par, TList>::run(int nsteps) {
  EulerMethod euler(h);
  WCAForce fwca(1, 1);
  double qh, qt;
  if (ratio < 0) {
    qh = std::sqrt(-ratio);
    qt = - 1 / qh;
  } else {
    qh = std::sqrt(ratio);
    qt = 1 / qh;
  }
  ExtDipoleForce fed(epsilon, qh, qt, r_cut, 4, 3.0 / 16);
  MySpatialSortingTraits<Par> sst;
  auto lambda = [&](int i, int j) {
    fed(f_arr[i], f_arr[j], p_arr[i], p_arr[j], pbc2, fwca);
  };
  for (int i = 1; i <= nsteps; i++) {
    if (spatial_sort_on)
      tlist->cal_force(p_arr, lambda, sst);
    else
      tlist->cal_force(p_arr, lambda);
    for (int ip = 0; ip < nPar; ip++) {
      euler.update_xy_uxuy(p_arr[ip], f_arr[ip], Pe, tau, pbc2, myran);
    }
    if (output_on) {
      //xy_out->write_double(i, p_arr);
      xy_out->write_theta(i, p_arr);
      (*log_out)(i);
    }
  }
  
}

#endif

