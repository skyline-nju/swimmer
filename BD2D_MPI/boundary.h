#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "vect.h"
#include "rand.h"
#include "cmdline.h"

/**
 * \brief Base class for boundary condition.
 */
class BoundaryBase_2 {
public:
  BoundaryBase_2() { lx_ = ly_ = x_min_ = y_min_ = 0; }
  explicit BoundaryBase_2(const cmdline::parser &cmd);
  BoundaryBase_2(double lx, double ly, double x0 = 0, double y0 = 0):
    lx_(lx), ly_(ly), x_min_(x0), y_min_(y0) { }
  double get_lx() const { return lx_; }
  double get_ly() const { return ly_; }
  double get_xmin() const {return x_min_;}
  double get_ymin() const {return y_min_;}
protected:
  void set_para0(double lx, double ly, double x0 = 0, double y0 = 0);
  double lx_;
  double ly_;
  double x_min_;
  double y_min_;
};

/**
 * \brief Periodic boundary condition in the x direction and free in the y direction. 
 */
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

template <typename T>
void PBC_x_2::wrap(T &p) const {
  if (p.x < x_min_) {
    p.x += lx_;
  } else if (p.x >= x_max_) {
    p.x -= lx_;
  }
}

/**
 * \brief Place two walls at x = xmin, xmax, and free boundary condition in y.
 */
class Wall_x_2 : public PBC_x_2 {
public:
  Wall_x_2() { k_ = xl_ = xr_= 0; }
  Wall_x_2(double k, double a, double lx, double ly, double x0=0, double y0=0)
    : PBC_x_2(lx, ly, x0, y0), k_(k), xl_(x0 + a), xr_(x0 + lx - a) {}
  explicit Wall_x_2(const cmdline::parser &cmd);
  void new_rand_pos(double &x, double &y, Ran &myran) const;
  template <typename T> void wrap(const T &p) const {}
  
  template<typename TPar>
  void force(TPar &p) const;
  
protected:
  double k_;
  double xl_;
  double xr_;
};

template<typename TPar>
void Wall_x_2::force(TPar &p) const {
  if (p.x < xl_) {
    p.fx += (xl_ - p.x) * k_;
  } else if (p.x > xr_) {
    p.fx += (xr_ - p.x) * k_;
  }
}

class Wall_x_PBC_y_2 : public Wall_x_2 {
public:
  Wall_x_PBC_y_2() { half_ly_ = y_max_ = 0; }
  Wall_x_PBC_y_2(double k, double a, double lx, double ly, double x0=0, double y0=0)
    : Wall_x_2(k, a, lx, ly, x0, y0), half_ly_(ly * 0.5), y_max_(y0+ly) { }
  explicit Wall_x_PBC_y_2(const cmdline::parser &cmd);

  void nearest_dis(Vec_2<double>& dis) const;

  template <typename T> void wrap(T &p) const;

protected:
  double half_ly_;
  double y_max_;
};


inline void PBC_x_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_lx_) {
    dis.x += lx_;
  } else if (dis.x > half_lx_) {
    dis.x -= lx_;
  }
}

inline void Wall_x_PBC_y_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.y < -half_ly_) {
    dis.y += ly_;
  } else if (dis.y > half_ly_) {
    dis.y -= ly_;
  }
}

template<typename T>
void Wall_x_PBC_y_2::wrap(T & p) const {
  if (p.y < y_min_) {
    p.y += ly_;
  } else if (p.y >= y_max_) {
    p.y -= ly_;
  }
}

/**
 * \brief Periodic boundary condition in both x, y direction.
 */
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


/**
* \brief Get the nearest distance between two particle
* \tparam T1 Template class for Partcile 1
* \tparam T2 Template class for Particle 2
* \tparam TBc Template class for the used boundary condition
* \param p1
* \param p2
* \param bc
* \return distance
*/
template <typename T1, typename T2, typename TBc>
double get_dis_square(const T1 &p1, const T2 &p2, const TBc &bc) {
  Vec_2<double> dis(p1.x - p2.x, p1.y - p2.y);
  bc.nearest_dis(dis);
  return dis.square();
}

template <typename T1, typename T2>
double get_dis_square(const T1 &p1, const T2 &p2, const Wall_x_2 &bc) {
  Vec_2<double> dis(p1.x - p2.x, p1.y - p2.y);
  return dis.square();
}
#endif


