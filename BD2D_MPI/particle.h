#ifndef PARTICLE_H
#define PARTICLE_H
#include "rand.h"
#include "comn.h"
#include "vect.h"
/********************** Create particles with random position ***************/
// check whether the new particle is overlapped with existed particle.
template <typename _TPar, typename _TBC>
bool check_overlap(const _TPar &p_new, const std::vector<_TPar> &p_arr,
    double sigma_square, const _TBC &bc) {
  bool flag_overlap = false;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    if (get_dis_square(p_new, *it, bc) < sigma_square) {
      flag_overlap = true;
      break;
    }
  }
  return flag_overlap;
}

// Add new particle with random position, avoiding overlapped with existed
// particle.
template <typename _TPar, typename _TBC>
bool add_new_rand_2(std::vector<_TPar> &p_arr, int n, double sigma_square,
                    Ran &myran, const _TBC &bc) {
  bool flag = false;
  int count = 0;
  while (count < 10 * n) {
    _TPar p_new(myran, bc);
    if (!check_overlap(p_new, p_arr, sigma_square, bc)) {
      p_arr.push_back(p_new);
      flag = true;
      break;
    } else {
      count++;
    }
  }
  return flag;
}

// Create particle randomly without overlapping
template <class _TPar, typename _TBC>
void create_rand_2(std::vector<_TPar> &p_arr, int n, double sigma,
                   Ran &myran, const _TBC &bc) {
  double sigma_square = sigma * sigma;
  for (int i = 0; i < n; i++) {
    if (!add_new_rand_2(p_arr, n, sigma, myran, bc)) {
      std::cout << "Failed to add the " << i << "-th particles" << std::endl;
      exit(1);
    }
  }
  std::cout << "Create " << n << " particles with diameter " << sigma
    << " and random positions" << std::endl;
}

/* Class for 2D Brownian particle */
class BrownPar_2 {
public:
  BrownPar_2() {}
  BrownPar_2(double x0, double y0) : x(x0), y(y0) { fx = fy = 0; }
  template <typename _TBC>
  BrownPar_2(Ran &myran, const _TBC &bc);

  double x;
  double y;
  double fx;
  double fy;
};

template<typename _TBC>
inline BrownPar_2::BrownPar_2(Ran &myran, const _TBC & bc) {
  bc.new_rand_pos(x, y, myran);
  fx = fy = 0;
}

/* Class for 2D Active Brownian particle */
class ActiveBrownPar_2{
public:
  ActiveBrownPar_2() {}
  ActiveBrownPar_2(double x0, double y0, double ux0, double uy0) :
    x(x0), y(y0), ux(ux0), uy(uy0), fx(0), fy(0) {}
  template <typename _TBC>
  ActiveBrownPar_2(Ran &myran, const _TBC &bc);

  double x;
  double y;
  double ux;
  double uy;
  double fx;
  double fy;
};


template<typename _TBC>
inline ActiveBrownPar_2::ActiveBrownPar_2(Ran &myran, const _TBC & bc) {
  bc.new_rand_pos(x, y, myran);
  double theta = 2 * PI * myran.doub();
  ux = std::cos(theta);
  uy = std::sin(theta);
  fx = fy = 0;
}
#endif
