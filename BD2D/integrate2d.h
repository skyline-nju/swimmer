#ifndef INTEGRATE2D_H
#define INTEGRATE2D_H
#include "rand.h"
#include "comn.h"

class Int_base_2 {
public:
  Int_base_2(double h0, double Lx0, double Ly0): h(h0), Lx(Lx0), Ly(Ly) {};

  void PBC(double &x, double &y) const;

protected:
  double h;
  double Lx;
  double Ly;
};

inline void Int_base_2::PBC(double &x, double &y) const {
  if (x < 0) {
    x += Lx;
  } else if (x >= Lx) {
    x -= Lx;
  }
  if (y < 0) {
    y += Ly;
  } else if (y >= Ly) {
    y -= Ly;
  }
}

// Integration of Stochastic Differential Equations by the Euler Method
class Int_EM_2 : public Int_base_2 {
public:
  Int_EM_2(double h0, double Lx0, double Ly0);

  template <class Par>
  void int_T(Par &p, Ran *myran) const;

  template <class Par>
  void int_R(Par &p, Ran *myran) const;

  template <class Par>
  void int_SP(Par &p, Ran *myran, double Pe) const;

  template <class Par>
  void int_SP(Par &p, Ran *myran, double Pe, double tau) const;

  template <class Par>
  void int_SP_all(Par *p, int N, Ran *myran, double Pe) const;

  template <class Par>
  void int_SP_all(Par *p, int N, Ran *myran, double Pe, double tau) const;

protected:
  double sqrt_24h;
  double sqrt_72h;
  double trible_h;
};

template <class Par>
inline void Int_EM_2::int_T(Par &p, Ran *myran) const {
  p.x += p.fx * h + (myran->doub() - 0.5) * sqrt_24h;
  p.y += p.fy * h + (myran->doub() - 0.5) * sqrt_24h;
  p.fx = p.fy = 0;
  PBC(p.x, p.y);
}

template <class Par>
inline void Int_EM_2::int_R(Par &p, Ran *myran) const {
  p.theta += p.tau * trible_h + (myran->doub() - 0.5) * sqrt_72h;
  p.tau = 0;
}

template <class Par>
inline void Int_EM_2::int_SP(Par &p, Ran *myran, double Pe) const {
  int_R(p, myran);
  p.fx += std::cos(p.theta) * Pe;
  p.fy += std::sin(p.theta) * Pe;
  int_T(p, myran);
}

template <class Par>
inline void Int_EM_2::int_SP(Par &p, Ran *myran, double Pe, double tau) const {
  p.tau += tau;
  int_SP(p, myran, Pe);
}

template <class Par>
void Int_EM_2::int_SP_all(Par *p, int N, Ran *myran, double Pe) const {
  for (int i = 0; i < N; i++) {
    int_SP(p[i], myran, Pe);
  }
}

template <class Par>
void Int_EM_2::int_SP_all(Par *p, int N, Ran *myran, double Pe, double tau) const {
  for (int i = 0; i < N; i++) {
    int_SP(p[i], myran, Pe, tau);
  }
}

#endif