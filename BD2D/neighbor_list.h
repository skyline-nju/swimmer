#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H
#include "comn.h"
#include "cell_list.h"

template <class Par>
class Par_w_Pre_Pos : public Par {
public:
  Par_w_Pre_Pos(): Par(), x_pre(0), y_pre(0) {}

  void update_pre_pos() { x_pre = this->x; y_pre = this->y; }


  double x_pre;
  double y_pre;
};
class NeighborList_2 {
public:
  NeighborList_2(double Lx0, double Ly0, double r_cut0);

  NeighborList_2(double Lx0, double Ly0, double r_cut0,
                double r_buf_ratio, int nPar);
  ~NeighborList_2();

  void set_r_buf(double r_buf_ratio, int nPar);
  
  void apply_PBC(double &dx, double &dy) const;

  template <class Par1, class Par2>
  double cal_spatial_dis_square(const Par1 &p1, const Par2 &p2) const;

  template <class Par>
  double cal_traveled_dis_square(const Par &p0) const;

  template <class Par>
  bool should_refresh(const Par *p, int nPar);

  template <class Par>
  void create(Par *p, int nPar);

  template <class Par>
  void recreate(Par *p, int nPar);

  template <class Par>
  void renew_pair(const Par *p, int i, int j);

  template <class Par>
  void update(Par *par, int nPar);


  template <class BiFunc>
  void for_each_pair(BiFunc f2);

protected:
  double Lx;
  double Ly;
  double half_Lx;
  double half_Ly;
  double r_buf;
  double half_r_buf_square;
  double r_verlet;
  double r_verlet_square;
  double r_cut;
  
  int * mat; // matrix for neighbor list, nPar rows by max_neighbor columns
  int nrows; // total number of particles
  int ncols; // max neighbors of one particle
  int ntot;  // size of mat

  int last_update_time;
  int count_refresh;
  double sum_interval;
  double sum_interval_square;
  int min_interval;
  int max_interval;

  int max_neighbor;

  CellList_w_list *cl;
};


inline void NeighborList_2::apply_PBC(double &dx, double &dy) const {
  if (dx > half_Lx)
    dx -= Lx;
  else if (dx < -half_Lx)
    dx += Lx;
  if (dy > half_Ly)
    dy -= Ly;
  else if (dy < -half_Ly)
    dy += Ly;
}

template<class Par1, class Par2>
inline double NeighborList_2::cal_spatial_dis_square(
  const Par1 & p1, const Par2 & p2) const {
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;
  apply_PBC(dx, dy);
  return (dx * dx + dy * dy);
}

template<class Par>
inline double NeighborList_2::cal_traveled_dis_square(
  const Par & p0) const {
  double dx = p0.x - p0.x_pre;
  double dy = p0.y - p0.y_pre;
  apply_PBC(dx, dy);
  return (dx * dx + dy * dy);
}

template<class Par>
bool NeighborList_2::should_refresh(const Par * p, int nPar) {
  bool flag = false;
  for (int i = 0; i < nPar; i++) {
    if (cal_traveled_dis_square(p[i]) > half_r_buf_square) {
      flag = true;
      break;
    }
  }
  return flag;
}

template<class Par>
void NeighborList_2::create(Par * par, int nPar) {
  for (int ip = 0; ip < nPar; ip++) {
    mat[ip * ncols] = 0;
    par[ip].update_pre_pos();
  }
  cl->create(par, nPar);
  auto lambda = [this, par](int pi, int pj) {
    renew_pair(par, pi, pj);
  };
  cl->for_each_pair(lambda);

  std::cout << "create neighbor list" << std::endl;
}

template<class Par>
void NeighborList_2::recreate(Par * par, int nPar) {
  for (int ip = 0; ip < nPar; ip++) {
    mat[ip * ncols] = 0;
    par[ip].update_pre_pos();
  }
  cl->update(par, nPar);
  auto lambda = [this, par](int pi, int pj) {
    renew_pair(par, pi, pj);
  };
  cl->for_each_pair(lambda);
}


template<class Par>
inline void NeighborList_2::renew_pair(const Par * par, int pi, int pj) {
  if (cal_spatial_dis_square(par[pi], par[pj]) < r_verlet_square) {
    int im = pi * ncols;
    mat[im]++;
    mat[im + mat[im]] = pj;
  }
}

template <class Par>
void NeighborList_2::update(Par *par, int nPar) {
  if (min_interval > 0  && last_update_time < 0.75 * min_interval) {
    last_update_time++;
  } else if (should_refresh(par, nPar)) {
    recreate(par, nPar);
    count_refresh++;
    sum_interval += last_update_time;
    sum_interval_square += last_update_time * last_update_time;
    if (last_update_time > max_interval)
      max_interval = last_update_time;
    if (min_interval == 0 || last_update_time < min_interval)
      min_interval = last_update_time;
    last_update_time = 0;
  } else {
    last_update_time++;
  }
}

template<class BiFunc>
void NeighborList_2::for_each_pair(BiFunc f2) {
  // row from 0 to nPar - 1
  for (int row = 0; row < nrows; row++) {
    int i0 = row * ncols;
    int len = mat[i0];
    if (len > max_neighbor)
      max_neighbor = len;
    for (int col = 1; col <= len; col++) {
      f2(row, mat[i0 + col]);
      //std::cout << row << std::endl;
    }
  }
}

#endif
