#ifndef VECT_H
#define VECT_H
#include <cmath>

template <typename T>
struct Vec_2 {
  T x;
  T y;

  Vec_2(): x(0), y(0) {}
  Vec_2(T x0, T y0) : x(x0), y(y0) {}
  Vec_2(const Vec_2<T> &a) : x(a.x), y(a.y) {}

  Vec_2<T>& operator =(const Vec_2<T>& rhs);
  void operator +=(const Vec_2<T>& rhs);
  void operator -=(const Vec_2<T>& rhs);
  void operator *=(double rhs);
  void operator /=(double rhs);
  Vec_2<T> operator +(const Vec_2<T>& rhs) const;
  Vec_2<T> operator -(const Vec_2<T>& rhs) const;
  Vec_2<T> operator +(double rhs) const;
  Vec_2<T> operator -(double rhs) const;
  Vec_2<T> operator *(double rhs) const;
  Vec_2<T> operator /(double rhs) const;

  double dot(const Vec_2<T>& a) const;
  double cross(const Vec_2<T>& a) const;
  void normalize();
  double square() const;
  Vec_2 inverse() const;
  void rotate(double theta);
};


template <typename T>
Vec_2<T>& Vec_2<T>::operator=(const Vec_2<T>& rhs) {
  this->x = rhs.x;
  this->y = rhs.y;
  return *this;
}

template <typename T>
void Vec_2<T>::operator+=(const Vec_2<T>& rhs) {
  x += rhs.x;
  y += rhs.y;
}

template <typename T>
void Vec_2<T>::operator-=(const Vec_2<T>& rhs) {
  x -= rhs.x;
  y -= rhs.y;
}

template <typename T>
void Vec_2<T>::operator*=(double rhs) {
  x *= rhs;
  y *= rhs;
}

template <typename T>
void Vec_2<T>::operator/=(double rhs) {
  x /= rhs;
  y /= rhs;
}

template <typename T>
Vec_2<T> Vec_2<T>::operator+(double rhs) const {
  return Vec_2<T>(x + rhs, y + rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator+(const Vec_2<T>& rhs) const {
  return Vec_2<T>(x + rhs.x, y + rhs.y);
}

template <typename T>
Vec_2<T> operator+(double a, const Vec_2<T>& b) {
  return Vec_2<T>(a + b.x, a + b.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-(double rhs) const {
  return Vec_2<T>(x - rhs, y - rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-(const Vec_2<T>& rhs) const {
  return Vec_2<T>(x - rhs.x, y - rhs.y);
}

template <typename T>
Vec_2<T> operator-(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs - rhs.x, lhs - rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator*(double rhs) const {
  return Vec_2<T>(x * rhs, y * rhs);
}

template <typename T>
Vec_2<T> operator*(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs * rhs.x, lhs * rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator/(double rhs) const {
  return Vec_2<T>(x / rhs, y / rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::inverse() const {
  return Vec_2<T>(1 / x, 1 / y);
}

template <typename T>
double Vec_2<T>::square() const { return x * x + y * y; }

template <typename T>
double Vec_2<T>::dot(const Vec_2<T>& a) const {
  return x * a.x + y * a.y;
}

template <typename T>
double Vec_2<T>::cross(const Vec_2<T>& a) const {
  return x * a.y - y * a.x;
}

template <typename T>
void Vec_2<T>::normalize() {
  const auto one_over_r = 1 / std::sqrt(square());
  x *= one_over_r;
  y *= one_over_r;
}

template <typename T>
void Vec_2<T>::rotate(double theta) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  const auto x_new = x * c - y * s;
  const auto y_new = x * s + y * c;
  x = x_new;
  y = y_new;
}

template <typename T>
struct Vec_3 {
  T x;
  T y;
  T z;

  Vec_3(): x(0), y(0), z(0) {}
  Vec_3(T x0, T y0, T z0) : x(x0), y(y0), z(z0) {}
  Vec_3(const Vec_2<T> &a) : x(a.x), y(a.y), z(a.z) {}

  Vec_3<T>& operator =(const Vec_3<T>& rhs);
  void operator +=(const Vec_3<T>& rhs);
  void operator -=(const Vec_3<T>& rhs);
  void operator *=(double rhs);
  void operator /=(double rhs);
  Vec_3<T> operator +(const Vec_3<T>& rhs) const;
  Vec_3<T> operator -(const Vec_3<T>& rhs) const;
  Vec_3<T> operator +(double rhs) const;
  Vec_3<T> operator -(double rhs) const;
  Vec_3<T> operator *(double rhs) const;
  Vec_3<T> operator /(double rhs) const;

  double dot(const Vec_3<T>& a) const;
  double cross(const Vec_3<T>& a) const;
  void normalize();
  double square() const;
  Vec_3 inverse() const;
  void rotate(double theta);
};

template <typename T>
Vec_3<T>& Vec_3<T>::operator=(const Vec_3<T>& rhs) {
  this->x = rhs.x;
  this->y = rhs.y;
  this->z = rhs.z;
  return *this;
}

template <typename T>
void Vec_3<T>::operator+=(const Vec_3<T>& rhs) {
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
}

template <typename T>
void Vec_3<T>::operator-=(const Vec_3<T>& rhs) {
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
}

template <typename T>
void Vec_2<T>::operator*=(double rhs) {
  x *= rhs;
  y *= rhs;
}

template <typename T>
void Vec_2<T>::operator/=(double rhs) {
  x /= rhs;
  y /= rhs;
}

template <typename T>
Vec_2<T> Vec_2<T>::operator+(double rhs) const {
  return Vec_2<T>(x + rhs, y + rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator+(const Vec_2<T>& rhs) const {
  return Vec_2<T>(x + rhs.x, y + rhs.y);
}

template <typename T>
Vec_2<T> operator+(double a, const Vec_2<T>& b) {
  return Vec_2<T>(a + b.x, a + b.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-(double rhs) const {
  return Vec_2<T>(x - rhs, y - rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator-(const Vec_2<T>& rhs) const {
  return Vec_2<T>(x - rhs.x, y - rhs.y);
}

template <typename T>
Vec_2<T> operator-(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs - rhs.x, lhs - rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator*(double rhs) const {
  return Vec_2<T>(x * rhs, y * rhs);
}

template <typename T>
Vec_2<T> operator*(double lhs, const Vec_2<T>& rhs) {
  return Vec_2<T>(lhs * rhs.x, lhs * rhs.y);
}

template <typename T>
Vec_2<T> Vec_2<T>::operator/(double rhs) const {
  return Vec_2<T>(x / rhs, y / rhs);
}

template <typename T>
Vec_2<T> Vec_2<T>::inverse() const {
  return Vec_2<T>(1 / x, 1 / y);
}

template <typename T>
double Vec_2<T>::square() const { return x * x + y * y; }

template <typename T>
double Vec_2<T>::dot(const Vec_2<T>& a) const {
  return x * a.x + y * a.y;
}

template <typename T>
double Vec_2<T>::cross(const Vec_2<T>& a) const {
  return x * a.y - y * a.x;
}

template <typename T>
void Vec_2<T>::normalize() {
  const auto one_over_r = 1 / std::sqrt(square());
  x *= one_over_r;
  y *= one_over_r;
}

template <typename T>
void Vec_2<T>::rotate(double theta) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  const auto x_new = x * c - y * s;
  const auto y_new = x * s + y * c;
  x = x_new;
  y = y_new;
}

#endif


