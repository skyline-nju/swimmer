#ifndef FORCE2D_H
#define FORCE2D_H
class F_base_2 {
public:
  F_base_2(double Lx0, double Ly0): Lx(Lx0), Ly(Ly0) {}

protected:
  double Lx;
  double Ly;
};

class F_WCA_2 : public F_base_2 {
public:
  F_WCA_2(double Lx0, double Ly0, double eps) :
    F_base_2(Lx0, Ly0), epsilon(eps) {
    
  }
protected:
  double epsilon;
};
#endif
