#ifndef DYNAMIC_H
#define DYNAMIC_H
#include <fstream>
#include <vector>
#include <iostream>
#include <functional>
#include "cmdline.h"
#include "integrate2d.h"
#include "force2d.h"
#include "particle2d.h"
#include "io_data.h"
#include "rand.h"
#include "cell_list.h"
#include "neighbor_list.h"

// visting all particles
template<class Par, class UFunc>
void for_each_particle(Par *p, int n, UFunc f1) {
  for (int i = 0; i < n; i++) {
    f1(p[i]);
  }
}

// visting all pairs by two loops over all particles
template<class Par, class BiFunc>
void for_each_pair(Par *p, int n, BiFunc f_pair) {
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      f_pair(p[i], p[j]);
    }
  }
}

class BaseDynamic_2 {
public:
  BaseDynamic_2(const cmdline::parser &cmd);
  virtual ~BaseDynamic_2();

  virtual void run(int steps) = 0;
  void run();

  template <class Par>
  void ini_rand(Par *p, double sigma);

protected:
  double Lx;
  double Ly;
  int nPar;
  int nstep;
  double packing_fraction;
  bool out_on;
  Ran *myran;
  F_WCA_2 * fwca;
  Int_EM_2 *euler;

  XY_Writer *xy_out;
  LogWriter *log_out;
  std::vector<std::ofstream> fout;
};

template <class Par>
void BaseDynamic_2::ini_rand(Par *p, double sigma) {
  if (packing_fraction < 0.5) {
    create_rand_2(p, nPar, sigma, Lx, Ly, myran);
  } else {
    double sigma0 = sigma;
    sigma = 0.5;
    create_rand_2(p, nPar, sigma, Lx, Ly, myran);
    auto force = [=](auto i, auto j) {
      fwca->pair(p, i, j);
    };
    auto integ = [=](auto &p1) {
      euler->int_T(p1, myran);
    };
    double d = 0.01;
    int n = int((sigma0 - sigma) / d);
    int delta_t = int(d / (euler->get_h()));
    std::cout << "The diameter is increased by " << d << " every "
              << delta_t << " time steps" << std::endl;

    CellList_w_list clist(Lx, Ly, fwca->get_r_cut());
    clist.create(p, nPar);
    for (int i = 0; i < n; i++) {
      sigma += d;
      fwca->set_sigma(sigma);
      for (int t = 0; t < delta_t; t++) {
        clist.for_each_pair(force);
        for_each_particle(p, nPar, integ);
        clist.update(p, nPar);
      }
    }
    std::cout << "The diameter is increased to " << sigma << std::endl;
  }
}

template<class Par, class List>
class BrownianDynamic : public BaseDynamic_2 {
public:
  BrownianDynamic(const cmdline::parser &cmd);
  ~BrownianDynamic();
  void run(int n_steps);

protected:
  Par * par;
  List *nlist;
};

template<class Par, class List>
BrownianDynamic<Par,List>::BrownianDynamic(const cmdline::parser & cmd):
                                      BaseDynamic_2(cmd) {
  par = new Par[nPar];
  double sigma = cmd.get<double>("sigma");
  ini_rand(par, 1);
  std::cout << typeid(Par).name() << std::endl;
  if (out_on)
    (*xy_out)(0, par);
  double l0 = fwca->get_r_cut();
  nlist = new List(Lx, Ly, l0);
  if (typeid(List) == typeid(CellList_w_list)) {

  } else {
    nlist->set_r_buf(0.4, nPar);
  }
  nlist->create(par, nPar);
  std::cout << "Finish initilization." << std::endl;
}

template<class Par, class List>
BrownianDynamic<Par, List>::~BrownianDynamic() {
  delete[] par;
  delete nlist;
}

template<class Par, class List>
void BrownianDynamic<Par, List>::run(int n_steps) {
  auto force = [=](auto i, auto j) {
    fwca->pair(par, i, j);
  };
  auto integ = [=](auto &p1) {
    euler->int_T(p1, myran);
  };

  for (int i = 1; i <= n_steps; i++) {
    nlist->for_each_pair(force);
    for_each_particle(par, nPar, integ);
    nlist->update(par, nPar);
    if (out_on) {
      (*xy_out)(i, par);
      (*log_out)(i);
    }
  }
}

#endif