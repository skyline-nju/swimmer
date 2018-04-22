#ifndef INTEGRATE_H
#define INTEGRATE_H
#include <rand.h>
#include "boundary.h"

class Run_and_tumble {
public:
  Run_and_tumble(double h, double alpha, double v0 = 1);

  template<typename TPar, typename TBc>
  void update(TPar &p, const TBc &bc, Ran &myran) const;

  template<typename TPar, typename TBc>
  void operator()(TPar &p, const TBc &bc, Ran &myran) const {
    update(p, bc, myran);
  }

  template<typename TPar>
  void operator()(TPar &p, const Wall_x_PBC_y_2 &bc, Ran &myran) const {
    bc.wall_force(p);
    update(p, bc, myran);
  }

  template <typename TPar, typename UniFunc1, typename UniFunc2>
  void operator() (TPar &p, Ran &myran, UniFunc1 wrap,
                   UniFunc2 wall_force) const;

private:
  double h_;
  double tumbling_rate_;
  double v0_;
  double reduced_tumbling_rate_;
  double mu_;
};

inline Run_and_tumble::Run_and_tumble(double h, double alpha, double v0) :
    h_(h), tumbling_rate_(alpha), v0_(v0) {
  reduced_tumbling_rate_ = h_ * tumbling_rate_;
  mu_ = 1;
}

template<typename TPar, typename TBc>
void Run_and_tumble::update(TPar & p, const TBc & bc, Ran & myran) const {
  if (myran.doub() < reduced_tumbling_rate_) {
    circle_point_picking(p.ux, p.uy, myran);
  }
  p.x += (p.ux * v0_ + mu_ * p.fx) * h_;
  p.y += (p.uy * v0_ + mu_ * p.fy) * h_;
  bc.wrap(p);
  p.fx = p.fy = 0;
}

template <typename TPar, typename UniFunc1, typename UniFunc2>
void Run_and_tumble::operator()(TPar& p, Ran& myran,
                                UniFunc1 wrap, UniFunc2 wall_force) const {
  wall_force(p);
  if (myran.doub() < reduced_tumbling_rate_) {
    circle_point_picking(p.ux, p.uy, myran);
  }
  p.x += (p.ux * v0_ + mu_ * p.fx) * h_;
  p.y += (p.uy * v0_ + mu_ * p.fy) * h_;
  wrap(p);
  p.fx = p.fy = 0;
}

#endif
