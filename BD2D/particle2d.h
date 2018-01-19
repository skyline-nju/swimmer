#ifndef PARTICLE2D_H
#define PARTICLE2D_H
#include <cmath>
class BP_2 {
public:
  BP_2() {}
  BP_2(double x0, double y0) : x(x0), y(y0) {}

  double x;
  double y;
};

class BP_with_ori_2 : public: PB_2{
public:
  BP_with_ori_2(): BP_2() {}
  BP_with_ori_2(double x0, double y0, double theta) : BP_2(x0, y0) {
    ux = std::cos(theta);
    uy = std::sin(theta);
}
  double ux;
  double uy;
}

#endif
