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

  void rot90(double r_value, double prob_sep);

  void rot90(double r_value, const double *prob_arr);

  void jump_rand(int &x_new, int &y_new, double r_val,
                 const double *prob_arr) const;

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

template <typename T1, typename T2>
void Par_2<T1, T2>::rot90(double r_value, double prob_sep) {
  auto tmp = ux;
  if (r_value < prob_sep) {
    ux = uy;
    uy = -tmp;
  } else {
    ux = -uy;
    uy = tmp;
  }
}

  template <typename T1, typename T2>
void Par_2<T1, T2>::rot90(double r_value, const double* prob_arr) {
  if (r_value < prob_arr[0]) {
    auto tmp = ux;
    if (r_value < prob_arr[1]) {
      ux = uy;
      uy = -tmp;
    } else {
      ux = -uy;
      uy = tmp;
    }
  }
}

template<typename T1, typename T2>
void Par_2<T1, T2>::jump_rand(int & x_new, int & y_new, double r_val,
                              const double * prob_arr) const {
  if (r_val < prob_arr[0]) {
    x_new = x+ux;
    y_new = y+uy;
  } else if (r_val < prob_arr[1]) {
    x_new = x-ux;
    y_new = y-uy;
  } else if (r_val < prob_arr[2]) {
    x_new = x-uy;
    y_new = y+ux;
  } else {
    x_new = x+uy;
    y_new = y-ux;
  }
}

const int state_2[4][2] = {{1, 0},       //! right 
                           {0, 1},       //! up
                           {-1, 0},      //! left
                           {0, -1}};     //! down

template <typename T>
class Par_s_2{
public:
  Par_s_2(): s(0) {}
  template <typename TRan>
  Par_s_2(TRan &myran, const Vec_2<int> &l,
        const Vec_2<int> &origin = Vec_2<int>());

  template <typename TRan>
  void tumble(TRan &myran) { s = myran.doub() * 4; }

  void rot90(double r_value, double prob_sep);

  void rot90(double r_value, const double* prob_arr);

  void jump_rand(int & x_new, int & y_new, double r_value,
                 const double * prob_arr) const;

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
void Par_s_2<T>::rot90(double r_value, double prob_sep) {
  if (r_value < prob_sep) {
    s -= 1;
    if (s < 0)
      s = 3;
  } else {
    s += 1;;
    if (s > 3)
      s = 0;
  }
}

template <typename T>
void Par_s_2<T>::rot90(double r_value, const double* prob_arr) {
  if (r_value < prob_arr[0]) {
    if (r_value < prob_arr[1]) {
      s -= 1;
      if (s < 0)
        s = 3;
    } else {
      s += 1;
      if (s > 3)
        s = 0;
    }
  }
}

template <typename T>
void Par_s_2<T>::jump_rand(int& x_new, int& y_new, double r_value,
                           const double* prob_arr) const {
  if (r_value < prob_arr[0]) {
    x_new = x + get_ux();
    y_new = y + get_uy();
  } else if (r_value < prob_arr[1]) {
    x_new = x - get_ux();
    y_new = y - get_uy();
  } else if (r_value < prob_arr[2]) {
    x_new = x - get_uy();
    y_new = y + get_ux();
  } else {
    x_new = x + get_uy();
    y_new = y - get_ux();
  }
}

} // end of namespace

template <typename T1, typename T2>
std::ostream &operator << (std::ostream &os, const lattice::Par_2<T1, T2> &p) {
  os << "x=" << int(p.x) << "\ty=" << int(p.y)
     << "\tux=" << int(p.ux) << "\tuy=" << int(p.uy);
  return os;
}
#endif
