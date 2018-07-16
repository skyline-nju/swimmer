#include "D2Q9.h"

lbm::D2Q9::Lattice::Lattice(int nx0, int ny0): nx(nx0), ny(ny0), size_scalar(nx * ny), size_pdf(size_scalar * ndir) {
  f0 = new double[size_scalar];
  f1 = new double[size_pdf - size_scalar];
  f2 = new double[size_pdf - size_scalar];
  rho = new double[size_scalar];
  ux = new double[size_scalar];
  uy = new double[size_scalar];
}

lbm::D2Q9::Lattice::~Lattice() {
  delete[] f0;
  delete[] f1;
  delete[] f2;
  delete[] rho;
  delete[] ux;
  delete[] uy;
}

void lbm::D2Q9::Lattice::init_pdf() const {
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      auto s_idx = scalar_idx(x, y);
      auto my_rho = rho[s_idx];
      auto my_ux = ux[s_idx];
      auto my_uy = uy[s_idx];
      auto uu = my_ux * my_ux + my_uy * my_uy;
      f0[scalar_idx(x, y)] = wi[0] * my_rho * (1.0 - 1.5 * uu);
      for (unsigned int i = 1; i < ndir; i++) {
        double ci_dot_u = dirx[i] * my_ux + diry[i] * my_uy;
        f1[pdfn_idx(x, y, i)] = wi[i] * my_rho * (1.0 + 3.0 * ci_dot_u + 4.5 * ci_dot_u * ci_dot_u - 1.5 * uu);
      }
    }
  }
}

void lbm::D2Q9::Lattice::compute_rho_u() const {
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      const auto s_idx = scalar_idx(x, y);
      double my_rho = f0[s_idx];
      double my_ux = 0;
      double my_uy = 0;
      for (int i = 1; i < ndir; i++) {
        const auto f_idx = pdfn_idx(x, y, i);
        my_rho += f1[f_idx];
        my_ux += dirx[i] * f1[f_idx];
        my_uy += diry[i] * f2[f_idx];
      }
      rho[s_idx] = my_rho;
      ux[s_idx] = my_ux / my_rho;
      uy[s_idx] = my_uy / my_rho;
    }
  }
}

void lbm::D2Q9::BGK::collide(double rho, double ux, double uy, double * f) const {
  const auto tw0r = omega_ * w0 * rho; // w[0]*rho/tau
  const auto twsr = omega_ * ws * rho; // w[1-4]*rho/tau
  const auto twdr = omega_ * wd * rho; // w[5-8]*rho/tau
  const auto omusq = 1.0 - 1.5*(ux * ux + uy * uy); // 1-(3/2)u.u
  const auto tux = 3.0 * ux;
  const auto tuy = 3.0 * uy;

  f[0] = omega_prime_ * f[0] + tw0r * (omusq);
  
  const auto c1_dot_3u = tux;
  f[1] = omega_prime_ * f[1] + twsr * (omusq + c1_dot_3u * (1 + 0.5 * c1_dot_3u));
  const auto c2_dot_3u = tuy;
  f[2] = omega_prime_ * f[2] + twsr * (omusq + c2_dot_3u * (1 + 0.5 * c2_dot_3u));
  const auto c3_dot_3u = -tux;
  f[3] = omega_prime_ * f[3] + twsr * (omusq + c3_dot_3u * (1 + 0.5 * c3_dot_3u));
  const auto c4_dot_3u = -tuy;
  f[4] = omega_prime_ * f[4] + twsr * (omusq + c4_dot_3u * (1 + 0.5 * c4_dot_3u));

  const auto c5_dot_3u = tux + tuy;
  f[5] = omega_prime_ * f[5] + twdr * (omusq + c5_dot_3u * (1 + 0.5 * c5_dot_3u));
  const auto c6_dot_3u = -tux + tuy;
  f[6] = omega_prime_ * f[6] + twdr * (omusq + c6_dot_3u * (1 + 0.5 * c6_dot_3u));
  const auto c7_dot_3u = -tux - tuy;
  f[7] = omega_prime_ * f[7] + twdr * (omusq + c7_dot_3u * (1 + 0.5 * c7_dot_3u));
  const auto c8_dot_3u = tux - tuy;
  f[8] = omega_prime_ * f[8] + twdr * (omusq + c8_dot_3u * (1 + 0.5 * c8_dot_3u));
}
