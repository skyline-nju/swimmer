#ifndef PARTICLE2D_H
#define PARTICLE2D_H
#include <cmath>
class BP_2d {
public:
  BP_2d() {}
  BP_2d(double x0, double y0) : x(x0), y(y0) {}

  double x;
  double y;
};

class BP_with_ori_2d : public: PB_2d{
public:
  BP_with_ori_2d(): BP_2d() {}
  BP_with_ori_2d(double x0, double y0, double theta) : BP_2d(x0, y0) {
    ux = std::cos(theta);
    uy = std::sin(theta);
}
  double ux;
  double uy;
}

#endif
