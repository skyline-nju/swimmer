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
    auto force = [=](auto &p1, auto &p2) {
      fwca->pair(p1, p2);
    };
    auto integ = [=](auto &p1) {
      euler->int_T(p1, myran);
    };
    double d = 0.01;
    int n = int((sigma0 - sigma) / d);
    int delta_t = int(d / (euler->get_h()));
    std::cout << "The diameter is increased by " << d << " every "
              << delta_t << " time steps" << std::endl;
    for (int i = 0; i < n; i++) {
      sigma += d;
      fwca->set_sigma(sigma);
      for (int t = 0; t < delta_t; t++) {
        for_each_pair(p, nPar, force);
        for_each_particle(p, nPar, integ);
      }
    }
    std::cout << "The diameter is increased to " << sigma << std::endl;
  }
}

template<class Par>
class BD_2: public BaseDynamic_2 {
public:
  BD_2(const cmdline::parser &cmd);
  ~BD_2() { delete[] par; }
  void run(int steps);

protected:
  Par *par;
  bool cell_list_on;
  Cell_list_2<Par> *clist;
};

template<class Par>
BD_2<Par>::BD_2(const cmdline::parser & cmd):
    BaseDynamic_2(cmd) {
  par = new Par[nPar];
  double sigma = cmd.get<double>("sigma");
  ini_rand(par, 1);
  std::cout << typeid(Par).name() << std::endl;
  if (typeid(Par) == typeid(BP_2)) {
    cell_list_on = false;
  } else {
    cell_list_on = true;
    double r_cut = fwca->get_r_cut();
    clist = new Cell_list_2<Par>(Lx, Ly, r_cut, r_cut);
  }
  if (out_on) 
    (*xy_out)(0, par);
  std::cout << "Finish initilization." << std::endl;
}

template<class Par>
void BD_2<Par>::run(int steps) {
  auto force = [=](auto &p1, auto &p2) {
    fwca->pair(p1, p2);
  };
  auto integ = [=](auto &p1) {
    euler->int_T(p1, myran);
  };
  for (int i = 1; i <= steps; i++) {
    if (cell_list_on) {
      clist->refresh(par, nPar);
      clist->for_each_pair(force);
    } else {
      for_each_pair(par, nPar, force);
    }
    for_each_particle(par, nPar, integ);
    if (out_on) {
      (*xy_out)(i, par);
      (*log_out)(i);
    }
  }
}

#endif
