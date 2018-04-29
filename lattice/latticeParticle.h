/**
 * @brief lattice particles
 * 
 * @file latticeParticle.h
 * @author skyline-nju
 * @date 2018-04-27
 */
#ifndef LATTICEPARTICLE_H
#define LATTICEPARTICLE_H
#include "vect.h"

namespace lattice {
/**
 * @brief Simple lattice particle
 * 
 * For saving memories, use uint16_t and int8_t for T1, T2 respectively.
 * 
 * @tparam T1  Integer type of coordinates
 * @tparam T2  Integer type of orientations
 */
template <typename T1, typename T2>
class Par_2 {
public:
  Par_2() = default;

  template <typename TRan>
  Par_2(TRan &myran, const Vec_2<int> &l,
        const Vec_2<int> &origin = Vec_2<int>());

  template <typename TRan>
  void tumble(TRan &myran);

  T2 get_ux() const { return ux; }
  T2 get_uy() const { return uy; }

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

const int state_2[4][2] = {
  {1, 0},       //! right 
  {0, 1},       //! up
  {-1, 0},      //! left
  {0, -1}       //! down
};

template <typename T>
class Par_s_2{
public:
  Par_s_2()=default;
  template <typename TRan>
  Par_s_2(TRan &myran, const Vec_2<int> &l,
        const Vec_2<int> &origin = Vec_2<int>());

  template <typename TRan>
  void tumble(TRan &myran);

  int get_ux()const { return state_2[s][0]; }
  int get_uy()const { return state_2[s][1]; }

  T x;
  T y;
  int8_t s;
};

template <typename T>
template <typename TRan>
Par_s_2<T>::Par_s_2(TRan& myran, const Vec_2<int>& l,
                    const Vec_2<int>& origin): s(0) {
  x = origin.x + int(myran.doub() * l.x);
  y = origin.y + int(myran.doub() * l.y);
  tumble(myran);
}

template <typename T>
template <typename TRan>
void Par_s_2<T>::tumble(TRan& myran) {
  s = myran.doub() * 4;
}

} // end of namespace

template <typename T1, typename T2>
std::ostream &operator << (std::ostream &os, const lattice::Par_2<T1, T2> &p) {
  os << "x=" << int(p.x) << "\ty=" << int(p.y)
     << "\tux=" << int(p.ux) << "\tuy=" << int(p.uy);
  return os;
}
#endif
