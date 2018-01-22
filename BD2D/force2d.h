#ifndef FORCE2D_H
#define FORCE2D_H

class F_base_2 {
public:
  F_base_2(double Lx0, double Ly0);
  
  // periodic boundary condition
  void PBC(double &dx, double &dy) const;

  // cal distance with periodic boundary condition
  template <class ParA, class ParB>
  void cal_dis_w_PBC(const ParA &a, const ParB &b, double &dx, double &dy) const;

protected:
  double Lx;
  double Ly;
  double half_Lx;
  double half_Ly;
};

// periodic boundary condition
inline void F_base_2::PBC(double &dx, double &dy) const {
  if (dx < -half_Lx) {
    dx += Lx;
  } else if (dx > half_Lx) {
    dx -= Lx;
  }
  if (dy < -half_Ly) {
    dy += Ly;
  } else if (dy > half_Ly) {
    dy -= Ly;
  }
}

// cal distance with periodic boundary condition
template <class ParA, class ParB>
inline void F_base_2::cal_dis_w_PBC(const ParA &a, const ParB &b,
                                    double &dx, double &dy) const {
  dx = a.x - b.x;
  dy = a.y - b.y;
  PBC(dx, dy);
}

/*****************************************************************************
*   The WCA force on particle 1 exerted by particle 2 with a distance r is   *
*   F_21 = 24 * epsilon * [2 (1/r)^12 - (1/r)^6] * r_21 / r^2,               *
*   where r_21 is a vector pointing from particle 2 to particle 1.           *
*****************************************************************************/

class F_WCA_2 : public F_base_2 {
public:
  F_WCA_2(double Lx0, double Ly0, double eps);

  template <class ParA, class ParB>
  int pair(ParA &a, ParB &b, double dx, double dy) const;

  template <class ParA, class ParB>
  int pair(ParA &a, ParB &b) const;

  template <class ParA, class ParB>
  int operator() (ParA &a, ParB &b) const { return pair(a, b); }

protected:
  double eps24;
  double r_cut;
  double r_cut_square;
};

template <class ParA, class ParB>
int F_WCA_2::pair(ParA &a, ParB &b, double dx, double dy) const {
  double dr_square = dx * dx + dy * dy;
  if (dr_square > r_cut_square) {
    return 0;
  } else {
    double r_2 = 1 / dr_square;       // 1 / dr ^ 2
    double r_6 = r_2 * r_2 * r_2;
    double tmp = eps24 * (2 * r_6 * r_6 - r_6) * r_2;
    double fx = tmp * dx;
    double fy = tmp * dy;
    a.fx += fx;
    a.fy += fy;
    b.fx -= fx;
    b.fy -= fy;
    return 1;
  }
}

template <class ParA, class ParB>
int F_WCA_2::pair(ParA &a, ParB &b) const {
  double dx, dy;
  cal_dis_w_PBC(a, b, dx, dy);
  return pair(a, b, dx, dy);
}

/*****************************************************************************
*                 Calculate all forces by visiting all pairs.                *
*****************************************************************************/
template<class Par, class Force>
void cal_pair_force_simply(Par *p, int n, const Force *force) {
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      force->pair(p[i], p[j]);
    }
  }
}

#endif
