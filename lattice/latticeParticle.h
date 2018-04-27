#ifndef LATTICEPARTICLE_H
#define LATTICEPARTICLE_H
#include "vect.h"

namespace lattice {
template <typename T1, typename T2>
class Par_2 {
public:
  Par_2() = default;

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l,
        const Vec_2<int> &origin = Vec_2<int>());

  template <typename TRan>
  void tumble(TRan &myran);

  T1 x;
  T1 y;
  T2 ux;
  T2 uy;
};

template <typename T1, typename T2>
template <typename TRan>
Par_2<T1, T2>::Par_2(TRan& myran, const Vec_2<int>& l, const Vec_2<int>& origin)
  : ux(0), uy(0) {
  x = origin.x + int(myran.doub() * l.x);
  y = origin.y + int(myran.doub() * l.y);
  tumble(myran);
}

template <typename T1, typename T2>
template<typename TRan>
void Par_2<T1, T2>::tumble(TRan & myran) {
  const auto state = int(myran.doub() * 4);
  switch (state) {
  case 0: ux = -1; uy = 0;  break;
  case 1: ux = 1;  uy = 0;  break;
  case 2: ux = 0;  uy = -1; break;
  case 3: ux = 0;  uy = 1;  break;
  default: break;
  }
}

} // end of namespace
#endif
