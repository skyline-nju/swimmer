#ifndef PARTICLE2D_H
#define PARTICLE2D_H
#include <typeinfo>
#include <iostream>
#include "rand.h"

class BP_2 {
public:
  BP_2() { x = y = fx = fy = 0; }
  BP_2(double x0, double y0) : x(x0), y(y0), fx(0), fy(0) {}

  virtual void set_data(double x0, double y0, Ran *myran);
  
  double x;
  double y;
  double fx;
  double fy;
};

class BP_with_ori_2: public BP_2 {
public:
  BP_with_ori_2(): BP_2() { theta = tau = 0; }
  BP_with_ori_2(double X, double Y, double Theta):
    BP_2(X, Y), theta(Theta), tau(0) {}
  virtual void set_data(double x0, double y0, Ran *myran);
  double theta;
  double tau;
};

/****************************************************************************/

// Check whether two particles are overlapped
bool overlap_2(double dx, double dy, double sigma, double Lx, double Ly);

template <class Par>
bool overlap_2(double x, double y, Par *p, int i,
               double sigma, double Lx, double Ly) {
  if (i == 0)
    return false;
  for (int j = 0; j < i; j++) {
    double dx = x - p[j].x;
    double dy = y - p[j].y;
    if (overlap_2(dx, dy, sigma, Lx, Ly))
      return true;
  }
  return false;
}

// Add new particle with random position, avoiding overlapped with existed
// particle.
template <class Par>
int add_new_rand_2(int i, Par *p, int n, double sigma,
                    double Lx, double Ly, Ran *myran) {
  int count = 0;
  while (count < 10 * n) {
    double x = myran->doub() * Lx;
    double y = myran->doub() * Ly;
    if (!overlap_2(x, y, p, i, sigma, Lx, Ly)) {
      p[i].set_data(x, y, myran);

      return 1;
    } else {
      count++;
    }
  }
  return 0;
}

// Create particle randomly without overlapping
template <class Par>
int create_rand_2(Par *p, int n, double sigma,
                  double Lx, double Ly, Ran *myran) {
  for (int i = 0; i < n; i++) {
    if (!add_new_rand_2(i, p, n, sigma, Lx, Ly, myran)) {
      std::cout << "Failed to add the " << i << "-th particles" << std::endl;
      exit(1);
    }
  }
  std::cout << "Create " << n << " particles with diameter " << sigma 
            << " and random positions" << std::endl;
  return 1;
}
#endif


