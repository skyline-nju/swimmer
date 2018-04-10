#ifndef FORCE_H
#define FORCE_H
#include <vector>
#include <cmath>
#include "vect.h"
#include "boundary.h"

class WCAForce {
public:
  WCAForce() {}
  WCAForce(double eps, double sigma = 1) { ini(eps, sigma); }

  void ini(double eps, double sigma = 1);

  void set_sigma(double sigma);

  double get_r_cut() { return std::sqrt(r_cut_square); }

  bool within_range(double dd) const { return dd < r_cut_square; }

  template <class T>
  bool operator()(T &fi, T &fj, const Vec_2<double> &dR) const;

  template <class T1, class T2>
  bool operator()(T1 &fi, T2 &fj, const Vec_2<double> &dR) const;

  template <class T>
  void eval(T &fi, T &fj, const Vec_2<double> &dR, double dR_square) const;

  template <class Force, class Par>
  void operator()(Force &fi, Force &fj, Par &pi, Par &pj, const PBC_2 &bc) const;

private:
  double eps24;
  double r_cut_square;
  double sigma_inverse_square;
};

template <class T>
inline bool WCAForce::operator() (T& fi, T& fj, const Vec_2<double>& dR) const {
  double dR_square = dR.square() * sigma_inverse_square;
  if (dR_square < r_cut_square) {
    double r_2 = 1 / dR_square;    //  dr ^ -2
    double r_6 = r_2 * r_2 * r_2;  //  dr ^ -6
    Vec_2<double> f(eps24 * (2 * r_6 * r_6 - r_6) * r_2 * dR);
    fi += f;
    fj -= f;
    return true;
  } else {
    return false;
  }
}

template <class T1, class T2>
inline bool WCAForce::operator() (T1& fi, T2& fj, const Vec_2<double>& dR) const {
  double dR_square = dR.square() * sigma_inverse_square;
  if (dR_square < r_cut_square) {
    double r_2 = 1 / dR_square;    //  dr ^ -2
    double r_6 = r_2 * r_2 * r_2;  //  dr ^ -6
    Vec_2<double> f(eps24 * (2 * r_6 * r_6 - r_6) * r_2 * dR);
    fi += f;
    fj -= f;
    return true;
  } else {
    return false;
  }
}

template<class T>
inline void WCAForce::eval(T & fi, T & fj, const Vec_2<double>& dR,
                           double dR_square) const {
  double r_2 = 1 / dR_square;    //  dr ^ -2
  double r_6 = r_2 * r_2 * r_2;  //  dr ^ -6
  Vec_2<double> f(eps24 * (2 * r_6 * r_6 - r_6) * r_2 * dR);
  fi += f;
  fj -= f;
}

template<class Force, class Par>
inline void WCAForce::operator()(Force & fi, Force & fj,
                                 Par & pi, Par & pj, const PBC_2 &bc) const {
  Vec_2<double> dR(pi.x - pj.x, pi.y - pj.y);
  bc.nearest_dis(dR);
  double dR_square = dR.square();
  if (dR_square < r_cut_square) {
    eval(fi, fj, dR, dR_square);
  }
}

class DipoleForce {
public:
  DipoleForce(double MM, double r_cut);

  void eval(Vec_3<double> &fi, Vec_3<double> &fj,
            const Vec_2<double> &dR, double dR_square,
            double theta_i, double theta_j) const;

  template <class Force, class Par>
  void operator() (Force &fi, Force &fj, Par &pi, Par &pj,
    const PBC_2 &bc, const WCAForce &fwca) const;
  

private:
  double mm;
  double r_cut_square;
};

template<class Force, class Par>
inline void DipoleForce::operator()(Force & fi, Force & fj, Par & pi, Par & pj,
                                    const PBC_2 & bc, const WCAForce & fwca) const {
  Vec_2<double> dR(pi.x - pj.x, pi.y - pj.y);
  bc.nearest_dis(dR);
  double dR_square = dR.square();
  if (dR_square < r_cut_square) {
    eval(fi, fj, dR, dR_square, pi.theta, pj.theta);
    if (fwca.within_range(dR_square)) {
      fwca.eval(fi, fj, dR, dR_square);
    }
  }
}

class ExtDipoleForce {
public:
  ExtDipoleForce(double eps, double qh, double qt,
                 double r_cut, double k, double d);

  template <class Force, class Par>
  void operator()(Force &fi, Force &fj, Par &pi, Par &pj,
    const PBC_2 &bc, const WCAForce &fwca) const;

  double cal_potential(const Vec_2<double> &dR,
    const Vec_2<double> &u1, const Vec_2<double> &u2) const;

  void eval(Vec_3<double> &f1, Vec_3<double> &f2, const Vec_2<double> &dR,
    const Vec_2<double> &u1, const Vec_2<double> &u2) const;

  void eval(Vec_3<double> &f1, Vec_3<double> &f2,
    const Vec_2<double> &dR, double theta1, double theta2) const;

  void eval2(Vec_3<double> &f1, Vec_3<double> &f2, const Vec_2<double> &dR,
    double theta1, double theta2) const;

  void eval3(Vec_3<double> &f1, Vec_3<double> &f2,
    const Vec_2<double> &dR, double theta1, double theta2) const;

  void charge_pair(Vec_2<double> &f, Vec_2<double> &tau, const Vec_2<double> &dR,
    double epsilon, double prod1, double prod2, double prod3) const;

private:
  double eps_qh_qh;
  double eps_qh_qt;
  double eps_qt_qt;
  double r_cut_square;
  double kappa;
  double delta;
};

inline double ExtDipoleForce::cal_potential(const Vec_2<double>& dR,
                                            const Vec_2<double> &u1,
                                            const Vec_2<double> &u2) const {
  double r_hh = std::sqrt((dR + delta * (u1 - u2)).square());
  double Uhh = eps_qh_qh * std::exp(-kappa * r_hh) / r_hh;
  double r_ht = std::sqrt((dR + delta * (u1 + u2)).square());
  double Uht = eps_qh_qt * std::exp(-kappa * r_ht) / r_ht;
  double r_th = std::sqrt((dR - delta * (u1 + u2)).square());
  double Uth = eps_qh_qt * std::exp(-kappa * r_th) / r_th;
  double r_tt = std::sqrt((dR - delta * (u1 - u2)).square());
  double Utt = eps_qt_qt * std::exp(-kappa * r_tt) / r_tt;
  return Uhh + Uht + Uth + Utt;
}

inline void ExtDipoleForce::charge_pair(Vec_2<double>& f, Vec_2<double>& tau,
  const Vec_2<double>& dR, double epsilon, double prod1, double prod2, double prod3) const {
  double r_square = dR.square();
  double r = std::sqrt(r_square);
  double tmp = epsilon * (kappa * r + 1) / r_square * std::exp(-kappa * r);
  f += tmp / r * dR;
  tmp *= delta / r;
  tau.x += tmp * (prod1 - prod3);
  tau.y += -tmp * (prod2 - prod3);
}


template<class Force, class Par>
inline void ExtDipoleForce::operator()(Force & fi, Force & fj, Par & pi, Par & pj,
  const PBC_2 & bc, const WCAForce & fwca) const {
  Vec_2<double> dR(pi.x - pj.x, pi.y - pj.y);
  bc.nearest_dis(dR);
  double dR_square = dR.square();
  if (dR_square < r_cut_square) {
    //eval(fi, fj, dR, pi.theta, pj.theta);
    eval(fi, fj, dR, pi.u, pj.u);
    if (fwca.within_range(dR_square)) {
      fwca.eval(fi, fj, dR, dR_square);
    }
  }
}

class AlignForce {
public:
  AlignForce(double epsilon, double r_cut);

  void eval(Vec_3<double> &fi, Vec_3<double> &fj,
            double theta_i, double theta_j) const;

  template <class Force, class Par>
  void operator() (Force &fi, Force &fj, Par &pi, Par &pj,
    const PBC_2 &bc, const WCAForce &fwca) const;


private:
  double eps;
  double r_cut_square;
};


template<class Force, class Par>
inline void AlignForce::operator()(Force &fi, Force &fj, Par &pi, Par &pj,
                                   const PBC_2 &bc, const WCAForce &fwca) const {

}
#endif
