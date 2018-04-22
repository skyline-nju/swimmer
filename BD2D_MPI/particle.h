#ifndef PARTICLE_H
#define PARTICLE_H
#include "rand.h"
#include "comn.h"
#include "vect.h"

inline void wrap_1(double &x, double x_min, double x_max, double l) {
  if (x < x_min) {
    x += l;
  } else if (x >= x_max) {
    x -= l;
  }
}

inline void unwrap_1(double &dx, double l, double half_l) {
  if (dx < -half_l) {
    dx += l;
  } else if (dx > half_l) {
    dx -= l;
  }
}

/* Check whether the new particle is overlapped with existed particle. */
template <typename TPar, typename TBc>
bool check_overlap(const TPar &p_new, const std::vector<TPar> &p_arr,
                   double sigma_square, const TBc &bc) {
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
      unwrap_1(dx, l.x, 0.5 * l.x);
    if (flag_bc.y == 0)
      unwrap_1(dy, l.y, 0.5 * l.y);
    if ( dx * dx + dy * dy < sigma_square) {
      flag_overlap = true;
      break;
    }
  }
  return flag_overlap;
}

/* Add new particle with random position, avoiding overlapped with existed. */
template <typename TPar, typename TBc>
bool add_new_rand_2(std::vector<TPar> &p_arr, int n, double sigma_square,
                    Ran &myran, const TBc &bc) {
  bool flag = false;
  int count = 0;
  while (count < 10 * n) {
    TPar p_new(myran, bc);
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

template <typename TPar>
bool add_new_rand_2(std::vector<TPar> &p_arr, int n, double sigma_square,
                    Ran &myran, const Vec_2<double> &l,
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
template <class TPar, typename TBc>
void create_rand_2(std::vector<TPar> &p_arr, int n, double sigma,
                   Ran &myran, const TBc &bc) {
  auto sigma_square = sigma * sigma;
  for (int i = 0; i < n; i++) {
    if (!add_new_rand_2(p_arr, n, sigma_square, myran, bc)) {
      std::cout << "Failed to add the " << i << "-th particles" << std::endl;
      exit(1);
    }
  }
  std::cout << "Create " << n << " particles with diameter " << sigma
    << " and random positions" << std::endl;
}

template <class TPar>
void create_rand_2(std::vector<TPar> &p_arr, int n, double sigma,
                   Ran &myran, const Vec_2<double> &l,
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

  template <typename TBc>
  BrownPar_2(Ran &myran, const TBc &bc);

  BrownPar_2(Ran &myran, const Vec_2<double> &l, const Vec_2<double> &origin)
    : x(origin.x + myran.doub() * l.x), y(origin.y + myran.doub() * l.y),
      fx(0), fy(0) {}

  double x;
  double y;
  double fx;
  double fy;
};

template<typename TBc>
BrownPar_2::BrownPar_2(Ran &myran, const TBc & bc) {
  x = y = 0;
  bc.new_rand_pos(x, y, myran);
  fx = fy = 0;
}

/**
 * \brief Class for Active Brownian particles in 2D
 */
class ActiveBrownPar_2: public BrownPar_2{
public:
  ActiveBrownPar_2() { ux = uy = 0; }

  ActiveBrownPar_2(double x0, double y0, double ux0, double uy0) :
    BrownPar_2(x0, y0), ux(ux0), uy(uy0) {}

  template <typename TBc>
  ActiveBrownPar_2(Ran &myran, const TBc &bc);

  ActiveBrownPar_2(Ran &myran, const Vec_2<double> &l,
                   const Vec_2<double> &origin);

  double ux;
  double uy;
};



template<typename TBc>
ActiveBrownPar_2::ActiveBrownPar_2(Ran &myran, const TBc & bc):
                                   BrownPar_2(myran, bc) {
  const auto theta = 2 * PI * myran.doub();
  ux = std::cos(theta);
  uy = std::sin(theta);
}

inline ActiveBrownPar_2::ActiveBrownPar_2(Ran& myran, const Vec_2<double>& l,
                                          const Vec_2<double>& origin)
  : BrownPar_2(myran, l, origin) {
  const auto theta = 2 * PI * myran.doub();
  ux = std::cos(theta);
  uy = std::sin(theta);
}
#endif
