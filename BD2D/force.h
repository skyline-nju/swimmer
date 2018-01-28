#ifndef FORCE_H
#define FORCE_H
#include <vector>
#include <cmath>
#include "vect.h"
#include "boundary.h"

class WCAForce {
public:
  WCAForce(double eps, double sigma = 1);

  void set_sigma(double sigma);

  double get_r_cut() { return std::sqrt(r_cut_square); }

  bool operator()(Vec_2<double> &fi, Vec_2<double> &fj,
            const Vec_2<double> &dR) const;

protected:
  double eps24;
  double r_cut_square;
  double sigma_inverse_square;
};

inline bool WCAForce::operator() (Vec_2<double>& fi, Vec_2<double>& fj,
                                  const Vec_2<double>& dR) const {
  double dR_square = dR.square() * sigma_inverse_square;
  if (dR_square < r_cut_square) {
    double r_2 = 1 / dR_square;    //  dr ^ -2
    double r_6 = r_2 * r_2 * r_2;  //  dr ^ -6
    double tmp = eps24 * (2 * r_6 * r_6 - r_6) * r_2;
    double fx = tmp * dR.x;
    double fy = tmp * dR.y;
    fi.x += fx;
    fi.y += fy;
    fj.x -= fx;
    fj.y -= fy;
    return true;
  } else {
    return false;
  }
}

#endif