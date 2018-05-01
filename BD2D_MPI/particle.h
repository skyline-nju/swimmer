#ifndef PARTICLE_H
#define PARTICLE_H
#include "rand.h"
#include "comn.h"
#include "vect.h"

/* Check whether the new particle is overlapped with existed particle. */
template <typename TPar>
bool check_overlap_2(const TPar &p_new, const std::vector<TPar> &p_arr,
                     double sigma_square, const Vec_2<double> &l,
                     const Vec_2<int> &flag_bc) {
  bool flag_overlap = false;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    double dx = p_new.x - (*it).x;
    double dy = p_new.y - (*it).y;
    if (flag_bc.x == 0)
      untangle_1(dx, l.x, 0.5 * l.x);
    if (flag_bc.y == 0)
      untangle_1(dy, l.y, 0.5 * l.y);
    if ( dx * dx + dy * dy < sigma_square) {
      flag_overlap = true;
      break;
    }
  }
  return flag_overlap;
}

/* Add new particle with random position, avoiding overlapped with existed. */
template <typename TPar>
bool add_new_rand_2(std::vector<TPar> &p_arr, int n, double sigma_square,
                    Ranq2 &myran, const Vec_2<double> &l,
                    const Vec_2<double> &origin,
                    const Vec_2<int> &flag_bc) {
  bool flag = false;
  int count = 0;
  while (count < 10 * n) {
    TPar p_new(myran, l, origin);
    if (!check_overlap_2(p_new, p_arr, sigma_square, l, flag_bc)) {
      p_arr.push_back(p_new);
      flag = true;
      break;
    } else {
      count++;
    }
  }
  return flag;
}

 /* Create particle randomly without overlapping. */
template <class TPar>
void create_rand_2(std::vector<TPar> &p_arr, int n, double sigma,
                   Ranq2 &myran, const Vec_2<double> &l,
                   const Vec_2<double> &origin,
                   const Vec_2<int> &flag_bc) {
  auto sigma_square = sigma * sigma;
  for (int i = 0; i < n; i++) {
    if (!add_new_rand_2(p_arr, n, sigma_square, myran, l, origin, flag_bc)) {
      std::cout << "Failed to add the " << i << "-th particles" << std::endl;
      exit(1);
    }
  }
  std::cout << "Create " << n << " particles with diameter " << sigma
    << " and random positions" << std::endl;
}

/**
 * \brief Class for 2D Brownain Particles
 */
class BrownPar_2 {
public:
  BrownPar_2():x(0), y(0), fx(0), fy(0) {}

  BrownPar_2(double x0, double y0) : x(x0), y(y0) { fx = fy = 0; }

  BrownPar_2(double x0, double y0, double fx0, double fy0)
    : x(x0), y(y0), fx(fx0), fy(fy0) {}

  BrownPar_2(Ranq2 &myran, const Vec_2<double> &l, const Vec_2<double> &origin)
    : x(origin.x + myran.doub() * l.x), y(origin.y + myran.doub() * l.y),
      fx(0), fy(0) {}

  double x;
  double y;
  double fx;
  double fy;
};

/**
 * \brief Class for Active Brownian particles in 2D
 */
class ActiveBrownPar_2: public BrownPar_2{
public:
  ActiveBrownPar_2() { ux = uy = 0; }

  ActiveBrownPar_2(double x0, double y0, double ux0, double uy0)
    : BrownPar_2(x0, y0), ux(ux0), uy(uy0) {}

  ActiveBrownPar_2(Ranq2 &myran, const Vec_2<double> &l,
                   const Vec_2<double> &origin);

  double ux;
  double uy;
};

inline ActiveBrownPar_2::ActiveBrownPar_2(Ranq2& myran, const Vec_2<double>& l,
                                          const Vec_2<double>& origin)
    : BrownPar_2(myran, l, origin) {
  const auto theta = 2 * PI * myran.doub();
  ux = std::cos(theta);
  uy = std::sin(theta);
}
#endif
