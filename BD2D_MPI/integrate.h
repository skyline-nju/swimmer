#ifndef INTEGRATE_H
#define INTEGRATE_H
#include "boundary.h"
#include "rand.h"

class Run_and_tumble {
public:
  Run_and_tumble(double h0, double alpha, double v = 1);

  template<typename _TPar, typename _TBC>
  void update(_TPar &p, const _TBC &bc, Ran &myran) const;

  template<typename _TPar, typename _TBC>
  void operator()(_TPar &p, const _TBC &bc, Ran &myran) const {
    update(p, bc, myran);
  }

private:
  double h;
  double tumbling_rate;
  double v0;
  double reduced_tumbling_rate;
  double mu;
};

inline Run_and_tumble::Run_and_tumble(double h0, double alpha, double v) :
    h(h0), tumbling_rate(alpha), v0(v) {
  reduced_tumbling_rate = h * tumbling_rate;
  mu = 1;
}

template<typename _TPar, typename _TBC>
inline void Run_and_tumble::update(_TPar & p, const _TBC & bc, Ran & myran) const {
  if (myran.doub() < reduced_tumbling_rate) {
    circle_point_picking(p.ux, p.uy, myran);
  }
  p.x += (p.ux + mu * p.fx) * h;
  p.y += (p.uy + mu * p.fy) * h;
  bc.wrap(p);
  p.fx = p.fy = 0;
}

#endif
