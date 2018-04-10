#ifndef SINGLE_DOMAIN_2_H
#define SIGNLE_DOMAIN_2_H
#include "cellList2D.h"
#include "rand.h"
#include "boundary.h"
#include "particle.h"
#include "integrate.h"
#include "force.h"

template <typename _TNode, typename _TBC>
class Single_domain_2 {
public:
  typedef unsigned long long int Ull;
  Single_domain_2(double Lx, double Ly, double phi, Ull seed0, double r_cut = 1);

  void ini_rand(double sigma = 1);

  template <typename _TForce>
  void cal_force(const _TForce &force);

  template <typename _TIntegrate>
  void integrate(const _TIntegrate &integ);

  template <typename _TIntegrate>
  void integrate2(const _TIntegrate &integ);

  template <typename _TIntegrate>
  void integrate3(const _TIntegrate &integ);

  void eval_elapsed_time(int n_steps, int t_start, int t_sep,
    int mode = 1, double h = 0.01);
private:
  Ran myran;
  Cell_list_2<_TNode> cell_list;
  std::vector<_TNode> p_arr;
  _TBC bc;

  int nPar;
  double packing_frac;
};

template<typename _TNode, typename _TBC>
Single_domain_2<_TNode, _TBC>::Single_domain_2(double Lx, double Ly,
    double phi, Ull seed0, double r_cut): myran(seed0),
    cell_list(Lx, Ly, 0, 0, r_cut), bc(Lx, Ly), packing_frac(phi) {
  nPar = cal_particle_number_2(packing_frac, Lx, Ly, r_cut);
  p_arr.reserve(nPar);
}

template<typename _TNode, typename _TBC>
void Single_domain_2<_TNode, _TBC>::ini_rand(double sigma) {
  create_rand_2(p_arr, nPar, sigma, myran, bc);
  cell_list.create(p_arr);
}

template<typename _TNode, typename _TBC>
template<typename _TForce>
inline void Single_domain_2<_TNode, _TBC>::cal_force(const _TForce &force) {
  auto f_ij = [this, &force](_TNode *pi, _TNode *pj) {
    force(*pi, *pj, bc);
  };
  cell_list.for_each_pair(f_ij);
}

template<typename _TNode, typename _TBC>
template<typename _TIntegrate>
void Single_domain_2<_TNode, _TBC>::integrate(const _TIntegrate &integ) {
  for (int i = 0; i < nPar; i++) {
    integ(p_arr[i], bc, myran);
  }
  cell_list.recreate(p_arr);
}

template<typename _TNode, typename _TBC>
template<typename _TIntegrate>
void Single_domain_2<_TNode, _TBC>::integrate2(const _TIntegrate &integ) {
  auto move = [this, &integ](_TNode *p) {
    integ(*p, bc, myran);
  };
  cell_list.update(move);
}

template<typename _TNode, typename _TBC>
template<typename _TIntegrate>
void Single_domain_2<_TNode, _TBC>::integrate3(const _TIntegrate &integ) {
  auto move = [this, &integ](_TNode *p) {
    integ(*p, bc, myran);
  };
  cell_list.update_by_row(move);
}

template<typename _TNode, typename _TBC>
inline void Single_domain_2<_TNode, _TBC>::eval_elapsed_time(
    int n_steps, int t_start, int t_sep, int mode, double h) {
  Run_and_tumble move(h, 0.0001);
  SpringForce f_spring(100);

  auto t1 = std::chrono::system_clock::now();
  int count = 0;
  for (int i = 0; i < n_steps; i++) {
    cal_force(f_spring);
    if (mode == 1)
      integrate(move);
    else if (mode == 3)
      integrate2(move);
    else if (mode == 2)
      integrate3(move);
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
  std::cout << std::endl;
}

#endif
