#ifndef FORCE_H
#define FORCE_H
#include "vect.h"
#include <cmath>

class SpringForce {
public:
  SpringForce() {}
  SpringForce(double k0, double sigma0 = 1):
    spring_const(k0), sigma(sigma0) {sigma_square = sigma * sigma;}
  void ini(double k0, double sigma0 = 1);
  void set_sigma(double sigma0);
  bool within_range(double dd) const { return dd < sigma_square; }

  template<typename _TPar1, typename _TPar2, typename _TBC>
  void eval(_TPar1 &p1, _TPar2 &p2, const _TBC &bc) const;

  template<typename _TPar1, typename _TPar2, typename _TBC>
  void operator()(_TPar1 &p1, _TPar2 &p2, const _TBC &bc) const {
    eval(p1, p2, bc);
  }

private:
  double spring_const;
  double sigma;
  double sigma_square;
};

inline void SpringForce::ini(double k0, double sigma0) {
  spring_const = k0;
  set_sigma(sigma0);
}

inline void SpringForce::set_sigma(double sigma0) {
  sigma = sigma0;
  sigma_square = sigma * sigma;
}

template<typename _TPar1, typename _TPar2, typename _TBC>
void SpringForce::eval(_TPar1 & p1, _TPar2 & p2, const _TBC &bc) const {
  Vec_2<double> R12_vec(p2.x - p1.x, p2.y - p1.y);
  bc.nearest_dis(R12_vec);
  double r12_square = R12_vec.square();
  if (r12_square < sigma_square) {
    double r12 = std::sqrt(r12_square);
    double f = (sigma - r12) * spring_const;
    double fx = f * R12_vec.x / r12;
    double fy = f * R12_vec.y / r12;
    p1.fx -= fx;
    p1.fy -= fy;
    p2.fx += fx;
    p2.fy += fy;
  }
}


#endif
