#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "vect.h"
#include "rand.h"
#include "cmdline.h"

template <typename T1, typename T2, typename TBc>
double get_dis_square(const T1 &p1, const T2 &p2, const TBc &bc) {
  Vec_2<double> dis(p1.x - p2.x, p1.y - p2.y);
  bc.nearest_dis(dis);
  return dis.square();
}

class BoundaryBase_2 {
public:
  BoundaryBase_2() { lx_ = ly_ = x_min_ = y_min_ = 0; }
  explicit BoundaryBase_2(const cmdline::parser &cmd);
  BoundaryBase_2(double lx, double ly, double x0 = 0, double y0 = 0):
    lx_(lx), ly_(ly), x_min_(x0), y_min_(y0) { }
  double get_Lx() const { return lx_; }
  double get_Ly() const { return ly_; }
protected:
  void set_para0(double lx, double ly, double x0 = 0, double y0 = 0);
  double lx_;
  double ly_;
  double x_min_;
  double y_min_;
};

class PBC_x_2 : public BoundaryBase_2 {
public:
  PBC_x_2() { half_lx_ = x_max_ = 0; }
  PBC_x_2(double lx, double ly, double x0 = 0, double y0 = 0);
  explicit PBC_x_2(const cmdline::parser &cmd);
  void set_para(double lx, double ly, double x0 = 0, double y0 = 0);
  void new_rand_pos(double &x, double &y, Ran &myran) const;
  void nearest_dis(Vec_2<double> &dis) const;
  template <typename T>
  void wrap(T &p) const;

protected:
  double half_lx_;
  double x_max_;

};

inline void PBC_x_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_lx_) {
    dis.x += lx_;
  } else if (dis.x > half_lx_) {
    dis.x -= lx_;
  }
}

template <typename T>
inline void PBC_x_2::wrap(T &p) const {
  if (p.x < x_min_) {
    p.x += lx_;
  } else if (p.x >= x_max_) {
    p.x -= lx_;
  }
}

class PBC_xy_2 : public BoundaryBase_2 {
public:
  PBC_xy_2() { half_lx_ = half_ly_ = x_max_ = y_max_ = 0; }
  PBC_xy_2(double lx, double ly, double x0 = 0, double y0 = 0);
  explicit PBC_xy_2(const cmdline::parser &cmd);
  void set_para(double lx, double ly, double x0 = 0, double y0 = 0);
  void new_rand_pos(double &x, double &y, Ran &myran) const;
  void nearest_dis(Vec_2<double> &dis) const;
  template <typename T>
  void wrap(T &p) const;

protected:
  double half_lx_;
  double half_ly_;
  double x_max_;
  double y_max_;
};

inline void PBC_xy_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_lx_) {
    dis.x += lx_;
  } else if (dis.x > half_lx_) {
    dis.x -= lx_;
  }
  if (dis.y < -half_ly_) {
    dis.y += ly_;
  } else if (dis.y > half_ly_) {
    dis.y -= ly_;
  }
}

template<typename T>
void PBC_xy_2::wrap(T & p) const {
  if (p.x < x_min_) {
    p.x += lx_;
  } else if (p.x >= x_max_) {
    p.x -= lx_;
  }
  if (p.y < x_min_) {
    p.y += ly_;
  } else if (p.y >= x_max_) {
    p.y -= ly_;
  }
}

#endif