#ifndef COMN_H
#define COMN_H
//#define USE_MPI

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>

const double PI = 3.14159265358979;

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

template <typename T1, typename T2>
void tangle_1(T1 &x, T1 x_min, T2 x_max, T2 len) {
  if (x < x_min) {
    x += len;
  } else if (x >= x_max) {
    x -= len;
  }
}

template <typename T>
void untangle_1(T &dx, T len, T half_len) {
  if (dx < -half_len) {
    dx += len;
  } else if (dx > half_len) {
    dx -= len;
  }
}

// create folder
void mkdir(const char *folder);

inline void mkdir(const std::string &folder) {mkdir(folder.c_str());}

// split string by a delimiter
std::vector<std::string> split(const std::string &str, const std::string &dlm);

template <class T>
void str_to_num(const std::string str, T &num) {
  std::stringstream ss;
  ss << str;
  ss >> num;
}

template <class T>
void num_to_str(const T &num, std::string str) {
  std::stringstream ss;
  ss << num;
  ss >> str;
}

// Calculate packing fraction from particle number in 2D domain
double cal_packing_fraction_2(int n, double Lx, double Ly, double sigma);

// Calculate particle number from packing fraction in 2D domain
int cal_particle_number_2(double phi, double Lx, double Ly, double sigma);

// Calculate the elapsed time to run a function
template <typename TFunc>
void cal_elapsed_time(TFunc f) {
  const auto t1 = std::chrono::system_clock::now();
  f();
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}

//#define SPATIAL_SORT
#ifdef SPATIAL_SORT
#include <CGAL/spatial_sort.h>
#include <CGAL/hilbert_sort.h>
struct MyLessX {
  template <class MyPoint>
  bool operator()(const MyPoint &p, const MyPoint &q) const {
    return p.x < q.x;
  }
};

struct MyLessY {
  template <class MyPoint>
  bool operator()(const MyPoint &p, const MyPoint &q) const {
    return p.y < q.y;
  }
};

template <class MyPoint>
struct MySpatialSortingTraits {
  typedef MyPoint Point_2;
  typedef MyLessX Less_x_2;
  typedef MyLessY Less_y_2;

  Less_x_2 less_x_2_object() const {
    return Less_x_2();
  }
  Less_y_2 less_y_2_object() const {
    return Less_y_2();
  }
};
#endif

#endif

