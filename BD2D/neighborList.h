#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
#include <vector>
#include "cellList.h"
#include "vect.h"
#include "boundary.h"
#include "comn.h"

class NeighborList_2 {
public:
  NeighborList_2(double Lx0, double Ly0, double r_cut0,
                double r_buf_ratio, int nPar);

  ~NeighborList_2();

  void count_time();

  bool is_safe();

  template <class Par>
  bool should_refresh(const std::vector<Par> &p) const;

  template <class Par>
  void create(std::vector<Par> &p);

  template <class Par, class BiFunc>
  void recreate(std::vector<Par> &p, BiFunc force, bool recreate_cl = false);

  template <class Par>
  void renew_mat(std::vector<Par> &p, int i, int j);

  template <class Par, class BiFunc>
  void renew_mat(std::vector<Par> &p, int i, int j, BiFunc force);

  template <class Par, class BiFunc>
  bool update(std::vector<Par> &p, BiFunc force);

  template <class BiFunc>
  void for_each_pair(BiFunc f2);

  template <class Par, class BiFunc>
  void cal_force(std::vector<Par> &p, BiFunc pair_force);

#ifdef SPATIAL_SORT
  template <class Par, class BiFunc>
  bool update(std::vector<Par> &p, BiFunc force,
              const MySpatialSortingTraits<Par> &sst);

  template<class Par, class BiFunc>
  void cal_force(std::vector<Par> &par, BiFunc pair_force,
                 const MySpatialSortingTraits<Par> &sst);
#endif

protected:
  double half_r_buf_square;
  double r_verlet_square;

  std::vector<int> mat;
  Vec_2<int> mat_bins;

  int refresh_interval;
  int refresh_count;
  int safe_inteval;
  double dt;
  double dt_square;
  int dt_min;
  int dt_max;
  int max_neighbor;

  PBC_2 pbc_2;

  CellList_list_2 *cl;
  std::vector<Vec_2<double>> last_pos;
};

inline bool NeighborList_2::is_safe() {
  return refresh_count >= 20 && refresh_interval < safe_inteval;
  //return false;
}

template<class Par>
inline bool NeighborList_2::should_refresh(const std::vector<Par>& p_arr) const {
  bool flag = false;
  for (int i = 0; i < mat_bins.y; i++) {
    if (pbc_2.nearest_dis_square(p_arr[i], last_pos[i]) > half_r_buf_square) {
      flag = true;
      break;
    }
  }
  return flag;
}

template<class Par>
inline void NeighborList_2::create(std::vector<Par>& p_arr) {
  for (int i = 0; i < mat_bins.y; i++) {
    mat[i * mat_bins.x] = 0;
    last_pos.emplace_back(p_arr[i].x, p_arr[i].y);
  }
  cl->create(p_arr);
  auto lambda = [this, &p_arr](int i, int j) {
    renew_mat(p_arr, i, j);
  };
  cl->for_each_pair(lambda);
}

template<class Par, class BiFunc>
inline void NeighborList_2::recreate(std::vector<Par>& p_arr, BiFunc force,
                                  bool recreate_cl) {
  for (int i = 0; i < mat_bins.y; i++) {
    mat[i * mat_bins.x] = 0;
    last_pos[i].x = p_arr[i].x;
    last_pos[i].y = p_arr[i].y;
  }
  if (recreate_cl) {
    cl->recreate(p_arr);
  } else {
    cl->update(p_arr);
  }
  auto lambda = [&, this](int i, int j) {
    renew_mat(p_arr, i, j, force);
  };
  cl->for_each_pair(lambda);
}

template<class Par>
inline void NeighborList_2::renew_mat(std::vector<Par>& p, int i, int j) {
  if (pbc_2.nearest_dis_square(p[i], p[j]) < r_verlet_square) {
    int im = i * mat_bins.x;
    mat[im]++;
    mat[im + mat[im]] = j;
  }
}

template<class Par, class BiFunc>
inline void NeighborList_2::renew_mat(std::vector<Par>& p, int i, int j,
                                    BiFunc force) {
  Vec_2<double> dR(p[i].x - p[j].x, p[i].y - p[j].y);
  pbc_2.nearest_dis(dR);
  if (dR.square() < r_verlet_square) {
    force(i, j);
    int im = i * mat_bins.x;
    mat[im]++;
    mat[im + mat[im]] = j;
  }
}

template<class Par, class BiFunc>
bool NeighborList_2::update(std::vector<Par>& p_arr, BiFunc force) {
  bool flag = false;
  if (is_safe()) {
    refresh_interval++;
  } else if (should_refresh(p_arr)) {
    recreate(p_arr, force);
    count_time();
    flag = true;
  } else {
    refresh_interval++;
  }
  return flag;
}

template<class BiFunc>
void NeighborList_2::for_each_pair(BiFunc f2) {
  for (int iy = 0; iy < mat_bins.y; iy++) {
    int i0 = mat_bins.x  * iy;
    int len = mat[i0];
    if (len > max_neighbor)
      max_neighbor = len;
    for (int ix = 1; ix <= len; ix++) {
      f2(iy, mat[i0 + ix]);
    }
  }
}

template<class Par, class BiFunc>
inline void NeighborList_2::cal_force(std::vector<Par>& p_arr, BiFunc f_ij) {
  if (!update(p_arr, f_ij))
    for_each_pair(f_ij);
}

#ifdef SPATIAL_SORT
template<class Par, class BiFunc>
bool NeighborList_2::update(std::vector<Par>& p_arr, BiFunc force,
                          const MySpatialSortingTraits<Par>& sst) {
  bool flag = false;
  if (is_safe()) {
    refresh_interval++;
  } else if (should_refresh(p_arr)) {
    if (refresh_count % 500 == 0) {
      CGAL::spatial_sort(p_arr.begin(), p_arr.end(), sst);
      recreate(p_arr, force, true);
    } else {
      recreate(p_arr, force, false);
    }
    count_time();
    flag = true;
  } else {
    refresh_interval++;
  }
  return flag;
}

template<class Par, class BiFunc>
inline void NeighborList_2::cal_force(std::vector<Par>& p_arr, BiFunc f_ij,
                                    const MySpatialSortingTraits<Par>& sst) {
  if (!update(p_arr, f_ij, sst))
    for_each_pair(f_ij);
}

#endif

#endif