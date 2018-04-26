#ifndef LATTICEPARTICLE_H
#define LATTICEPARTICLE_H
#include <cstdint>
#include "vect.h"

class LatticePar_2 {
public:
  LatticePar_2() = default;

  template <typename TRan>
  LatticePar_2(TRan &myran, const Vec_2<int> &l,
               const Vec_2<int> &origin = Vec_2<int>());

  template <typename TRan>
  void tumble(TRan &myran);

  int x;
  int y;
  int8_t ux;
  int8_t uy;
};

template <typename TRan>
LatticePar_2::LatticePar_2(TRan& myran, const Vec_2<int>& l,
                           const Vec_2<int>& origin): ux(0), uy(0) {
  x = origin.x + int(myran.doub() * l.x);
  y = origin.y + int(myran.doub() * l.y);
  tumble(myran);
}

template<typename TRan>
void LatticePar_2::tumble(TRan & myran) {
  const auto state = int(myran.doub() * 4);
  switch (state) {
  case 0: ux = -1; uy = 0;  break;
  case 1: ux = 1;  uy = 0;  break;
  case 2: ux = 0;  uy = -1; break;
  case 3: ux = 0;  uy = 1;  break;
  default: break;
  }
}
#endif
