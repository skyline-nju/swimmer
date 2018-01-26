#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H
#include "comn.h"
#include "cell_list.h"
#include <vector>

class NeighborList_2 {
public:
  NeighborList_2(double Lx0, double Ly0, double r_cut0);

  NeighborList_2(double Lx0, double Ly0, double r_cut0,
                double r_buf_ratio, int nPar);
  ~NeighborList_2();

  void set_r_buf(double r_buf_ratio, int nPar);
  
  void apply_PBC(double &dx, double &dy) const;

  template <class Par1, class Par2>
  double cal_dis_square(const Par1 &p1, const Par2 &p2) const;

  template <class Par>
  bool should_refresh(const std::vector<Par> &p);

  template <class Par>
  void create(std::vector<Par> &p);

  template <class Par, class BiFunc>
  void recreate(std::vector<Par> &p, BiFunc force,
                bool recreate_cell_list=false);

  template <class Par>
  void renew_pair(std::vector<Par> &p, int i, int j);

  template <class Par, class BiFunc>
  void renew_pair(std::vector<Par> &p, int i, int j, BiFunc force);

  template <class Par, class BiFunc>
  bool update(std::vector<Par> &p, BiFunc force);

  template <class BiFunc>
  void for_each_pair(BiFunc f2);

  template <class Par, class BiFunc>
  void cal_force(std::vector<Par> &p, BiFunc pair_force);

#ifdef SPATIAL_SORT
  template <class Par, class BiFunc>
  bool update(std::vector<Par> &p, BiFunc force,
              MySpatialSortingTraits<Par> &sst);

  template<class Par, class BiFunc>
  void cal_force(std::vector<Par> &par,
                 BiFunc pair_force,
                 MySpatialSortingTraits<Par> &sst);
#endif

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
  std::vector<Vec_2<double>> last_pos;
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
inline double NeighborList_2::cal_dis_square(const Par1& p1, const Par2 & p2) const {
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;
  apply_PBC(dx, dy);
  return (dx * dx + dy * dy);
}

template<class Par>
bool NeighborList_2::should_refresh(const std::vector<Par> &par) {
  bool flag = false;
  for (int i = 0; i < nrows; i++) {
    if (cal_dis_square(par[i], last_pos[i])> half_r_buf_square) {
      flag = true;
      break;
    }
  }
  return flag;
}

template<class Par>
void NeighborList_2::create(std::vector<Par> &par) {
  for (int ip = 0; ip < nrows; ip++) {
    mat[ip * ncols] = 0;
    last_pos.emplace_back(par[ip].x, par[ip].y);
  }
  cl->create(par);
  auto lambda = [this, &par](int pi, int pj) {
    renew_pair(par, pi, pj);
  };
  cl->for_each_pair(lambda);
}

template<class Par, class BiFunc>
void NeighborList_2::recreate(std::vector<Par> &par, BiFunc force,
                              bool recreate_cell_list) {
  for (int ip = 0; ip < nrows; ip++) {
    mat[ip * ncols] = 0;
    last_pos[ip].x = par[ip].x;
    last_pos[ip].y = par[ip].y;
  }
  if (recreate_cell_list) {
    cl->recreate(par);
  } else {
    cl->update(par);
  }
  auto lambda = [this, &par, &force](int pi, int pj) {
    renew_pair(par, pi, pj, force);
  };
  cl->for_each_pair(lambda);
}

template<class Par>
inline void NeighborList_2::renew_pair(std::vector<Par> &p, int pi, int pj) {
  if (cal_dis_square(p[pi], p[pj]) < r_verlet_square) {
    int im = pi * ncols;
    mat[im]++;
    mat[im + mat[im]] = pj;
  }
}

template<class Par, class BiFunc>
inline void NeighborList_2::renew_pair(std::vector<Par> &p,
                                       int pi, int pj, BiFunc force) {
  if (cal_dis_square(p[pi], p[pj]) < r_verlet_square) {
    force(pi, pj);
    int im = pi * ncols;
    mat[im]++;
    mat[im + mat[im]] = pj;
  }  
}

template <class Par, class BiFunc>
bool NeighborList_2::update(std::vector<Par> &par, BiFunc force) {
  if (min_interval > 0 && last_update_time < 0.75 * min_interval) {
    last_update_time++;
    return false;
  } else if (should_refresh(par)) {
    recreate(par, force);
    count_refresh++;
    sum_interval += last_update_time;
    sum_interval_square += last_update_time * last_update_time;
    if (last_update_time > max_interval)
      max_interval = last_update_time;
    if (min_interval == 0 || last_update_time < min_interval)
      min_interval = last_update_time;
    last_update_time = 0;
    return true;
  } else {
    last_update_time++;
    return false;
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
    }
  }
}
#ifdef SPATIAL_SORT
template <class Par, class BiFunc>
bool NeighborList_2::update(std::vector<Par> &par, BiFunc force,
  MySpatialSortingTraits<Par> &sst) {
  if (min_interval > 0 && last_update_time < 0.75 * min_interval) {
    last_update_time++;
    return false;
  } else if (should_refresh(par)) {
    if (count_refresh % 500 == 0) {
      CGAL::spatial_sort(par.begin(), par.end(), sst);
      recreate(par, force, true);
    } else {
      recreate(par, force, false);
    }
    count_refresh++;
    sum_interval += last_update_time;
    sum_interval_square += last_update_time * last_update_time;
    if (last_update_time > max_interval)
      max_interval = last_update_time;
    if (min_interval == 0 || last_update_time < min_interval)
      min_interval = last_update_time;
    last_update_time = 0;
    return true;
  } else {
    last_update_time++;
    return false;
  }
}

template<class Par, class BiFunc>
inline void NeighborList_2::cal_force(std::vector<Par> &par, BiFunc pair_force) {
  if (!update(par, pair_force))
    for_each_pair(pair_force);
}

template<class Par, class BiFunc>
inline void NeighborList_2::cal_force(std::vector<Par> &par,
                                      BiFunc pair_force,
                                      MySpatialSortingTraits<Par> &sst) {

  if (!update(par, pair_force, sst))
    for_each_pair(pair_force);
}
#endif

#endif
