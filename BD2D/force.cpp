#include "force.h"
#include <iostream>

void WCAForce::ini(double eps, double sigma) {
  eps24 = eps * 24;
  double r_cut = std::pow(2.0, 1. / 6.);
  r_cut_square = r_cut * r_cut;
  set_sigma(sigma);
}

void WCAForce::set_sigma(double sigma) {
  sigma_inverse_square = 1 / (sigma * sigma);
}

DipoleForce::DipoleForce(double MM, double r_cut) {
  mm = MM;
  r_cut_square = r_cut * r_cut;
}

void DipoleForce::eval(Vec_3<double>& f1, Vec_3<double>& f2,
                       const Vec_2<double>& dR, double dR_square,
                       double theta1, double theta2) const {
  double inverse_dR = 1 / std::sqrt(dR_square);
  Vec_2<double> dR_hat(dR.x * inverse_dR, dR.y * inverse_dR);
  Vec_2<double> mu1(std::cos(theta1), std::sin(theta1));
  Vec_2<double> mu2(std::cos(theta2), std::sin(theta2));

  double mu1_cross_mu2 = - mu1.cross(mu2);
  double mu1_dot_dR_hat = mu1.dot(dR_hat);
  double mu2_dot_dR_hat = mu2.dot(dR_hat);

  Vec_2<double> f(3 * mm / (dR_square * dR_square) * (
    (mu1.dot(mu2) - 5 * mu1_dot_dR_hat * mu2_dot_dR_hat) * dR_hat +
    mu2_dot_dR_hat * mu1 + mu1_dot_dR_hat * mu2));
  f1 += f;
  f2 -= f;

  double tmp = mm / dR_square * inverse_dR;
  f1.z += tmp * (mu1_cross_mu2 + 3 * mu1.cross(dR_hat) * mu2_dot_dR_hat);
  f2.z += tmp * (-mu1_cross_mu2 + 3 * mu2.cross(dR_hat) * mu1_dot_dR_hat);
}

ExtDipoleForce::ExtDipoleForce(double eps, double qh, double qt,
                               double r_cut, double k, double d) {
  eps_qh_qh = eps * qh * qh;
  eps_qh_qt = eps * qh * qt;
  eps_qt_qt = eps * qt * qt;
  r_cut_square = r_cut * r_cut;
  kappa = k;
  delta = d;
  std::cout << "d = " << d << std::endl;
  std::cout << "kappa = " << k << std::endl;

}

void ExtDipoleForce::eval(Vec_3<double>& f1, Vec_3<double>& f2, const Vec_2<double> dR,
                          double theta1, double theta2) const {
  Vec_2<double> u1(std::cos(theta1), std::sin(theta1));
  Vec_2<double> u2(std::cos(theta2), std::sin(theta2));
  double prod1 = u1.cross(dR);
  double prod2 = u2.cross(dR);
  double prod3 = delta * u1.cross(u2);

  Vec_2<double> f;
  Vec_2<double> tau;
  charge_pair(f, tau, dR + delta * (u1 - u2), eps_qh_qh,
    prod1, prod2, prod3);
  charge_pair(f, tau, dR + delta * (u1 + u2), eps_qh_qt,
    prod1, -prod2, -prod3);
  charge_pair(f, tau, dR - delta * (u1 + u2), eps_qh_qt,
    -prod1, prod2, -prod3);
  charge_pair(f, tau, dR - delta * (u1 - u2), eps_qt_qt,
    -prod1, -prod2, prod3);
  f1 += f;
  f2 -= f;
  f1.z += tau.x;
  f2.z += tau.y;
}

void ExtDipoleForce::eval2(Vec_3<double>& f1, Vec_3<double>& f2, const Vec_2<double> dR,
                           double theta1, double theta2) const {
  double dx = 1e-10;
  Vec_2<double> u1(std::cos(theta1), std::sin(theta1));
  Vec_2<double> u2(std::cos(theta2), std::sin(theta2));
  Vec_2<double> u1_new(std::cos(theta1 + dx), std::sin(theta1 + dx));
  Vec_2<double> u2_new(std::cos(theta2 + dx), std::sin(theta2 + dx));
  Vec_2<double> dR2(dR.x + dx, dR.y);
  double u0 = cal_potential(dR, u1, u2);
  double ux = cal_potential(Vec_2<double>(dR.x + dx, dR.y), u1, u2);
  double uy = cal_potential(Vec_2<double>(dR.x, dR.y + dx), u1, u2);
  double u_t1 = cal_potential(dR, u1_new, u2);
  double u_t2 = cal_potential(dR, u1, u2_new);
  Vec_2<double> f((u0 - ux) / dx, (u0 - uy) / dx);
  f1 += f;
  f2 -= f;
  f1.z += (u0 - u_t1) / dx;
  f2.z += (u0 - u_t2) / dx;

}

void ExtDipoleForce::eval3(Vec_3<double>& f1, Vec_3<double>& f2,
  const Vec_2<double> dR, double theta1, double theta2) const {
  Vec_2<double> u1(std::cos(theta1), std::sin(theta1));
  Vec_2<double> u2(std::cos(theta2), std::sin(theta2));

  double u1_cross_dR = u1.cross(dR);
  double u2_cross_dR = u2.cross(dR);
  double u1_cross_u2 = u1.cross(u2);
  Vec_2<double> vec_r(dR + delta * (u1 - u2));
  double r_square = vec_r.square();
  double r = std::sqrt(r_square);
  double tmp = eps_qh_qh * (kappa * r + 1) / r_square * std::exp(-kappa * r);
  Vec_2<double> f(tmp / r * vec_r);
  f1 += f;
  f2 -= f;
  tmp *= delta / r;
  f1.z +=  tmp * (u1_cross_dR - delta * u1_cross_u2);
  f2.z += -tmp * (u2_cross_dR - delta * u1_cross_u2);

  vec_r = dR + delta * (u1 + u2);
  r_square = vec_r.square();
  r = std::sqrt(r_square);
  tmp = eps_qh_qt * (kappa * r + 1) / r_square * std::exp(-kappa * r);
  f = tmp / r * vec_r;
  f1 += f;
  f2 -= f;
  tmp *= delta / r;
  f1.z += tmp * (u1_cross_dR + delta * u1_cross_u2);
  f2.z += tmp * (u2_cross_dR - delta * u1_cross_u2);

  vec_r = dR - delta * (u1 + u2);
  r_square = vec_r.square();
  r = std::sqrt(r_square);
  tmp = eps_qh_qt * (kappa * r + 1) / r_square * std::exp(-kappa * r);
  f = tmp / r * vec_r;
  f1 += f;
  f2 -= f;
  tmp *= delta / r;
  f1.z += -tmp * (u1_cross_dR - delta * u1_cross_u2);
  f2.z += -tmp * (u2_cross_dR + delta * u1_cross_u2);

  vec_r = dR - delta * (u1 - u2);
  r_square = vec_r.square();
  r = std::sqrt(r_square);
  tmp = eps_qt_qt * (kappa * r + 1) / r_square * std::exp(-kappa * r);
  f = tmp / r * vec_r;
  f1 += f;
  f2 -= f;
  tmp *= delta / r;
  f1.z += -tmp * (u1_cross_dR + delta * u1_cross_u2);
  f2.z +=  tmp * (u2_cross_dR + delta * u1_cross_u2);
}
