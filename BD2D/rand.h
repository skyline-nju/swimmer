#ifndef RAND_H
#define RAND_H
#include <cmath>

// Ref: Numerical Recipes, The Art of Scientific Computing, 3rd, p366-367
// Generate a random number
struct Ran
{
  unsigned long long u, v, w;
  Ran(unsigned long long j) :v(4101842887655102017LL), w(1) {
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
  }
  inline unsigned long long int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }
  inline double doub() {
    return 5.42101086242752217E-20 * int64();
  }
  inline unsigned int int32() {
    return (unsigned int)int64();
  }
  template<typename T>
  void circle_point_picking(T &x, T &y);

  template<typename T>
  void sphere_point_picking(T &x, T &y, T &z);

  template<typename T>
  void hypersphere_point_picking(T &X);

  template<typename T>
  void shuffle(T *a, int n);
};

// Uniform distribution of points on the circumference of a unit circle
// Ref: http://mathworld.wolfram.com/CirclePointPicking.html
template<typename T>
void Ran::circle_point_picking(T &x, T &y) {
  double a, b, aa, bb, S;
  do {
    a = doub() * 2 - 1;
    b = doub() * 2 - 1;
    aa = a * a;
    bb = b * b;
    S = aa + bb;
  } while (S >= 1);
  x = (aa - bb) / S;
  y = 2 * a * b / S;
}

// Uniform distribution of points on the surface of a unit sphere
// Ref: http://mathworld.wolfram.com/SpherePointPicking.html
template<typename T>
void Ran::sphere_point_picking(T &x, T &y, T &z) {
  double a, b, S;
  do {
    a = doub() * 2 - 1;
    b = doub() * 2 - 1;
    S = a * a + b * b;
  } while (S >= 1);
  double R = std::sqrt(1 - S);
  x = 2 * a * R;
  y = 2 * b * R;
  z = 1 - 2 * S;
}

// Unifor distribution of points on the surface of a unit 4d sphere
// Ref: http://mathworld.wolfram.com/HyperspherePointPicking.html
template<typename T>
void Ran::hypersphere_point_picking(T &X) {
  double a, b, S1, S2;
  do {
    a = doub() * 2 - 1;
    b = doub() * 2 - 1;
    S1 = a * a + b * b;
  } while (S1 >= 1);
  X[0] = a;
  X[1] = b;
  do {
    a = doub() * 2 - 1;
    b = doub() * 2 - 1;
    S2 = a * a + b * b;
  } while (S2 >= 1);
  double Q = std::sqrt((1 - S1) / S2);
  X[2] = a * Q;
  X[3] = b * Q;
}

// For other cases of random point picking, ref.:
// http://mathworld.wolfram.com/topics/RandomPointPicking.html

// Shuffle a array randomly
template<class T>
void Ran::shuffle(T *a, int n) {
  for (int i = n - 1; i >= 0; i--) {
    // generate a random int j that 0 <= j <= i  
    int j = int(doub() * (i + 1));
    if (j > i)
      j = i;
    else if (j < 0)
      j = 0;
    T tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }
}
#endif


