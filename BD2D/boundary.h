#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "vect.h"

struct PBC_2 {
  Vec_2<double> L;
  Vec_2<double> half_L;

  PBC_2() : L(), half_L() {}
  PBC_2(double Lx, double Ly) :
    L(Lx, Ly), half_L(0.5 * Lx, 0.5 * Ly) {}
  PBC_2(const Vec_2<double> &L0) :
    L(L0), half_L(0.5 * L0.x, 0.5 * L0.y) {}
  void ini(double Lx, double Ly);

  void nearest_dis(Vec_2<double> &dis) const;

  template<class Par1, class Par2>
  void nearest_dis(Vec_2<double> &dis, const Par1 &p1, const Par2 &p2) const;

  template<class Par1, class Par2>
  double nearest_dis_square(const Par1 &p1, const Par2 &p2) const;

  template <class Par>
  void wrap(Par &R) const;
};

inline void PBC_2::ini(double Lx, double Ly) {
  L.x = Lx;
  L.y = Ly;
  half_L.x = 0.5 * Lx;
  half_L.y = 0.5 * Ly;
}

inline void PBC_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_L.x) {
    dis.x += L.x;
  } else if (dis.x > half_L.x) {
    dis.x -= L.x;
  }
  if (dis.y < -half_L.y) {
    dis.y += L.y;
  } else if (dis.y > half_L.y) {
    dis.y -= L.y;
  }
}

template<class Par1, class Par2>
inline void PBC_2::nearest_dis(Vec_2<double>& dis,
  const Par1 & p1, const Par2 & p2) const {
  dis.x = p1.x - p2.x;
  dis.y = p1.y - p2.y;
  return nearest_dis(dis);
}

template<class Par1, class Par2>
inline double PBC_2::nearest_dis_square(const Par1 &p1, const Par2 &p2) const {
  Vec_2<double> dR(p1.x - p2.x, p1.y - p2.y);
  nearest_dis(dR);
  return dR.square();
}

template <class Par>
inline void PBC_2::wrap(Par &R) const {
  if (R.x < 0) {
    R.x += L.x;
  } else if (R.x >= L.x) {
    R.x -= L.x;
  }
  if (R.y < 0) {
    R.y += L.y;
  } else if (R.y >= L.y) {
    R.y -= L.y;
  }
}
#endif
