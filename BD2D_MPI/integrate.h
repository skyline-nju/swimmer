/**
 * @brief Integration methods for swimmers.
 * 
 * @file integrate.h
 * @author your name
 * @date 2018-04-26
 */
#ifndef INTEGRATE_H
#define INTEGRATE_H
#include "rand.h"

class Run_and_tumble {
public:
  Run_and_tumble(double h, double alpha, double v0 = 1., double mu = 1.)
    : h_(h), tumbling_rate_(alpha * h), v0_(v0), mu_(mu) {}

  template <typename TPar, typename UniFunc>
  void move (TPar &p, Ranq2 &myran, UniFunc tangle) const;

  template<typename TPar, typename UniFunc>
  void operator() (TPar &p, Ranq2 &myran, UniFunc tangle) const {
    move(p, myran, tangle);
  }
  template <typename TPar, typename UniF1, typename UniF2>
  void operator() (TPar &p, Ranq2 &myran, UniF1 tangle, UniF2 wall_force) const {
    wall_force(p);
    move(p, myran, tangle);
  }
private:
  double h_;
  double tumbling_rate_;
  double v0_;
  double mu_;
};

template <typename TPar, typename UniFunc>
void Run_and_tumble::move (TPar& p, Ranq2& myran, UniFunc tangle) const {
  if (myran.doub() < tumbling_rate_) {
    circle_point_picking(p.ux, p.uy, myran);
  }
  p.x += (p.ux * v0_ + mu_ * p.fx) * h_;
  p.y += (p.uy * v0_ + mu_ * p.fy) * h_;
  tangle(p);
  p.fx = p.fy = 0;
}

#endif
