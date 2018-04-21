#ifndef PARTICLE_H
#define PARTICLE_H
#include <rand.h>
#include <comn.h>

/**
 * \brief check whether the new particle is overlapped with existed particle.
 * \tparam TPar
 * \tparam TBc 
 * \param p_new 
 * \param p_arr 
 * \param sigma_square 
 * \param bc 
 * \return flag_overlap
 */
template <typename TPar, typename TBc>
bool check_overlap(const TPar &p_new, const std::vector<TPar> &p_arr,
    double sigma_square, const TBc &bc) {
  bool flag_overlap = false;
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    if (get_dis_square(p_new, *it, bc) < sigma_square) {
      flag_overlap = true;
      break;
    }
  }
  return flag_overlap;
}

/**
 * \brief Add new particle with random position, avoiding overlapped with existed
 * particles
 * \tparam TPar 
 * \tparam TBc 
 * \param p_arr 
 * \param n 
 * \param sigma_square 
 * \param myran 
 * \param bc 
 * \return flag
 */
template <typename TPar, typename TBc>
bool add_new_rand_2(std::vector<TPar> &p_arr, int n, double sigma_square,
                    Ran &myran, const TBc &bc) {
  bool flag = false;
  int count = 0;
  while (count < 10 * n) {
    TPar p_new(myran, bc);
    if (!check_overlap(p_new, p_arr, sigma_square, bc)) {
      p_arr.push_back(p_new);
      flag = true;
      break;
    } else {
      count++;
    }
  }
  return flag;
}

/**
 * \brief Create particle randomly without overlapping
 * \tparam TPar 
 * \tparam TBc 
 * \param p_arr 
 * \param n 
 * \param sigma 
 * \param myran 
 * \param bc 
 */
template <class TPar, typename TBc>
void create_rand_2(std::vector<TPar> &p_arr, int n, double sigma,
                   Ran &myran, const TBc &bc) {
  auto sigma_square = sigma * sigma;
  for (int i = 0; i < n; i++) {
    if (!add_new_rand_2(p_arr, n, sigma_square, myran, bc)) {
      std::cout << "Failed to add the " << i << "-th particles" << std::endl;
      exit(1);
    }
  }
  std::cout << "Create " << n << " particles with diameter " << sigma
    << " and random positions" << std::endl;
}

/**
 * \brief Class for Brownian particles in 2D
 */
class BrownPar_2 {
public:
  BrownPar_2():x(0), y(0), fx(0), fy(0) {}
  BrownPar_2(double x0, double y0) : x(x0), y(y0) { fx = fy = 0; }
  template <typename TBc>
  BrownPar_2(Ran &myran, const TBc &bc);

  double x;
  double y;
  double fx;
  double fy;
};

template<typename TBc>
BrownPar_2::BrownPar_2(Ran &myran, const TBc & bc) {
  x = y = 0;
  bc.new_rand_pos(x, y, myran);
  fx = fy = 0;
}

/**
 * \brief Class for Active Brownian particles in 2D
 */
class ActiveBrownPar_2{
public:
  ActiveBrownPar_2() { x = y = ux = uy = fx = fy = 0; }
  ActiveBrownPar_2(double x0, double y0, double ux0, double uy0) :
    x(x0), y(y0), ux(ux0), uy(uy0), fx(0), fy(0) {}
  template <typename TBc>
  ActiveBrownPar_2(Ran &myran, const TBc &bc);

  double x;
  double y;
  double ux;
  double uy;
  double fx;
  double fy;
};


template<typename TBc>
ActiveBrownPar_2::ActiveBrownPar_2(Ran &myran, const TBc & bc) {
  double x0, y0;
  bc.new_rand_pos(x0, y0, myran);
  x = x0;
  y = y0;
  const auto theta = 2 * PI * myran.doub();
  ux = std::cos(theta);
  uy = std::sin(theta);
  fx = fy = 0;
}
#endif
