#ifndef INTEGRATE_H
#define INTEGRATE_H
#include "boundary.h"
#include "rand.h"
class EulerMethod {
public:
  EulerMethod(double h0);

  template <class Par>
  void update_xy(Par &p, Vec_2<double> &f,
                 const PBC_2 &pbc2, Ran *myran) const;

protected:
  double h;
  double sqrt_24h;
  double sqrt_72h;
  double trible_h;
};



template<class Par>
inline void EulerMethod::update_xy(Par & p, Vec_2<double>& f,
                                   const PBC_2 &pbc2, Ran * myran) const {
  p.x += f.x * h + (myran->doub() - 0.5) * sqrt_24h;
  p.y += f.y * h + (myran->doub() - 0.5) * sqrt_24h;
  pbc2.wrap(p);
  f.x = 0;
  f.y = 0;
}

#endif
