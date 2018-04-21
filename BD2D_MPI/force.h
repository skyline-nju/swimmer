#ifndef FORCE_H
#define FORCE_H
#include <vect.h>
#include <cmath>

class SpringForce {
public:
  SpringForce() { spring_const_ = sigma_ = sigma_square_ = 0; }
  explicit SpringForce(double k, double sigma = 1):
    spring_const_(k), sigma_(sigma), sigma_square_(sigma * sigma) {}
  void ini(double k, double sigma = 1);       
  void set_sigma(double sigma0);
  bool within_range(double dd) const { return dd < sigma_square_; }

  template<typename TPar1, typename TPar2, typename TBc>
  void eval(TPar1 &p1, TPar2 &p2, const TBc &bc) const;

  template<typename TPar1, typename TPar2, typename TBc>
  void operator()(TPar1 &p1, TPar2 &p2, const TBc &bc) const {
    eval(p1, p2, bc);
  }

private:
  double spring_const_;
  double sigma_;
  double sigma_square_;
};

inline void SpringForce::ini(double k, double sigma) {
  spring_const_ = k;
  set_sigma(sigma);
}

inline void SpringForce::set_sigma(double sigma0) {
  sigma_ = sigma0;
  sigma_square_ = sigma_ * sigma_;
}

template<typename TPar1, typename TPar2, typename TBc>
void SpringForce::eval(TPar1 & p1, TPar2 & p2, const TBc &bc) const {
  Vec_2<double> r12_vec(p2.x - p1.x, p2.y - p1.y);
  bc.nearest_dis(r12_vec);
  const auto r12_square = r12_vec.square();
  if (r12_square < sigma_square_) {
    const auto r12 = std::sqrt(r12_square);
    const auto f = (sigma_ - r12) * spring_const_;
    const auto fx = f * r12_vec.x / r12;
    const auto fy = f * r12_vec.y / r12;
    p1.fx -= fx;
    p1.fy -= fy;
    p2.fx += fx;
    p2.fy += fy;
  }
}


#endif
