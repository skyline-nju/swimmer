#ifndef VECT_H
#define VECT_H

template <typename T>
struct Vec_2 {
  T x, y;

  Vec_2(): x(0), y(0) {}
  Vec_2(T a) : x(a), y(a) {}
  Vec_2(const Vec_2<T> &a) : x(a.x), y(a.y) {}
  Vec_2(T x0, T y0) : x(x0), y(y0) {}
  
  inline double square() const { return x * x + y * y; }
};

template <typename T>
struct Vec_3 {
  T x, y, z;

  Vec_3() : x(0), y(0), z(0) {};
  Vec_3(T a) : x(a), y(a), z(a) {};
  Vec_3(const Vec_3<T> &a) : x(a.x), y(a.y), z(a.z) {};
  Vec_3(T x0, T y0, T z0) : x(x0), y(y0), z(z0) {}

  inline double square() const { return x * x + y * y + z * z; }

};

#endif
