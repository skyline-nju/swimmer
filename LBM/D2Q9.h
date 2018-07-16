/**
 * @brief Lattice model D2Q9 for the lattice Boltzmann method
 * 
 * @file D2Q9.h
 * @author skyline-nju
 * @date 2018-06-29
 */
#pragma once

namespace lbm {

namespace D2Q9 {
// the lattice weights
const double w0 = 4.0 / 9.0;   // zero weight
const double ws = 1.0 / 9.0;   // adjacent weight
const double wd = 1.0 / 36.0;  // diagonal weight

// arrays of the lattice weights and direction components
// direction numbering scheme
// 6 2 5
// 3 0 1
// 7 4 8
const double wi[] = {w0, ws, ws, ws, ws, wd, wd, wd, wd};
const int dirx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int diry[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const int ndir = 9;

class BGK {
public:
  BGK(double tau): omega_(1 / tau), omega_prime_(1 - omega_) {}
  void collide(double rho, double ux, double uy, double *f) const;
private:
  double omega_;
  double omega_prime_;
};

class Lattice {
public:
  Lattice(int nx0, int ny0);
  ~Lattice();

  size_t scalar_idx(int x, int y) const { return nx * y + x;}

  // index for rest particle distribution function f0
  size_t pdf0_idx(int x, int y) const { return scalar_idx(x, y);}

  // index for non-rest particle distribution function f1-8
  size_t pdfn_idx(int x, int y, int dir) const { return nx * (ny * (dir-1) + y) + x;}

  void init_pdf() const;

  void compute_rho_u() const;

  template <class TCollide_model>
  void stream_collide(const TCollide_model &cm, bool flag_save);

  double get_rho(int x, int y) const { return rho[scalar_idx(x, y)];}

  double get_ux(int x, int y) const { return ux[scalar_idx(x, y)];}

  double get_uy(int x, int y) const { return uy[scalar_idx(x, y)];}

  int nx;
  int ny;
  size_t size_scalar;
  size_t size_pdf;
  double *f0;
  double *f1;
  double *f2;
  double *rho;
  double *ux;
  double *uy;
};

template<class TCollide_model>
void Lattice::stream_collide(const TCollide_model & cm, bool flag_save) {
  for (int y = 0; y < ny; y++) {
    int yp1 = y + 1;
    if (yp1 >= ny)
      yp1 = 0;
    int ym1 = y - 1;
    if (ym1 < 0)
      ym1 += ny;
    for (int x = 0; x < nx; x++) {
      int xp1 = x + 1;
      if (xp1 >= nx)
        xp1 = 0;
      int xm1 = x - 1;
      if (xm1 < 0)
        xm1 += nx;

      // load populations from adjacent nodes
      double f[] = {
        f0[pdf0_idx(x, y)],
        f1[pdfn_idx(xm1, y, 1)],
        f1[pdfn_idx(x, ym1, 2)],
        f1[pdfn_idx(xp1, y, 3)],
        f1[pdfn_idx(x, yp1, 4)],
        f1[pdfn_idx(xm1, ym1, 5)],
        f1[pdfn_idx(xp1, ym1, 6)],
        f1[pdfn_idx(xp1, yp1, 7)],
        f1[pdfn_idx(xm1, yp1, 8)]
      };

      // compute moments
      double my_rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
      double rhoinv = 1.0 / my_rho;
      double my_ux = rhoinv * (f[1] + f[5] + f[8] - (f[3] + f[6] + f[7]));
      double my_uy = rhoinv * (f[2] + f[5] + f[6] - (f[4] + f[7] + f[8]));

      // only write to memory when needed
      if (flag_save) {
        rho[scalar_idx(x, y)] = my_rho;
        ux[scalar_idx(x, y)] = my_ux;
        uy[scalar_idx(x, y)] = my_uy;
      }

      // collide
      cm.collide(my_rho, my_ux, my_uy, f);
      f0[pdf0_idx(x, y)] = f[0];
      f2[pdfn_idx(x, y, 1)] = f[1];
      f2[pdfn_idx(x, y, 2)] = f[2];
      f2[pdfn_idx(x, y, 3)] = f[3];
      f2[pdfn_idx(x, y, 4)] = f[4];
      f2[pdfn_idx(x, y, 5)] = f[5];
      f2[pdfn_idx(x, y, 6)] = f[6];
      f2[pdfn_idx(x, y, 7)] = f[7];
      f2[pdfn_idx(x, y, 8)] = f[8];
    }
  }
  // swap f1, f2
  double *tmp = f2;
  f2 = f1;
  f1 = tmp;
}

} // namespace D2Q9
  
} // namespace lbm