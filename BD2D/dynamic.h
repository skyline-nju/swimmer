#ifndef DYNAMIC_H
#define DYNAMIC_H
#include <fstream>
#include <vector>
#include "cmdline.h"
#include "integrate2d.h"
#include "force2d.h"
#include "particle2d.h"
#include "io_data.h"
#include "rand.h"

// visting all particles
template<class Par, class UFunc>
void for_each_particle(Par *p, int n, UFunc f1) {
  for (int i = 0; i < n; i++) {
    f1(&p[i]);
  }
}

// visting all pairs by two loop over all particles
template<class Par, class BiFunc>
void for_each_pair(Par *p, int n, BiFunc f_pair) {
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      f_pair(p[i], p[j]);
    }
  }
}


class BD_2 {
public:
  BD_2(const cmdline::parser &cmd);
  ~BD_2();
  void set_output(const cmdline::parser &cmd);
  virtual void create_particles(const cmdline::parser &cmd);
  virtual void run();

protected:
  double Lx;
  double Ly;
  int nPar;
  int nstep;
  double packing_fraction;
  Ran *myran;

  BP_2 *par;
  F_WCA_2 * fwca;
  Int_EM_2 *euler;

  XY_Writer *xy_out;
  LogWriter *log_out;
  std::vector<std::ofstream> fout;
};
#endif
