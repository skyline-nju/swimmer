#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
#include <list>
#include <iostream>
#include "vect.h"
#include "comn.h"

class CellListBase_2 {
public:
  CellListBase_2(double Lx, double Ly, double rcut);

  template <class Par>
  int get_ic(const Par &p) const;

  virtual bool has_member(int ic) const = 0;

  template <class UniFunc>
  void for_each_nearest_cell(UniFunc f) const;

protected:
  Vec_2<double> l_cell;
  Vec_2<int> bins;
  int ncells;
};

template <class Par>
int CellListBase_2::get_ic(const Par &p) const {
  return int(p.x / l_cell.x) + int(p.y / l_cell.y) * bins.x;
}

template<class UniFunc>
void CellListBase_2::for_each_nearest_cell(UniFunc f1) const {
  int idx_cell[5];
  for (int iy = 0; iy < bins.y; iy++) {
    int iy_up = iy + 1;
    if (iy_up == bins.y)
      iy_up = 0;
    int iy_times_nx = iy * bins.x;
    int iy_up_times_nx = iy_up * bins.x;

    for (int ix = 0; ix < bins.x; ix++) {
      int my_i = ix + iy_times_nx;
      if (has_member(my_i)) {
        int ix_rt = ix + 1;
        if (ix_rt == bins.x)
          ix_rt = 0;
        int ix_lt = ix - 1;
        if (ix_lt == -1)
          ix_lt = bins.x - 1;

        idx_cell[0] = my_i;
        idx_cell[1] = ix_rt + iy_times_nx;
        idx_cell[2] = ix_lt + iy_up_times_nx;
        idx_cell[3] = ix + iy_up_times_nx;
        idx_cell[4] = ix_rt + iy_up_times_nx;
        f1(idx_cell);
      }
    }
  }
}

class CellList_list_2 : public CellListBase_2 {
public:
  CellList_list_2(double Lx, double Ly, double rcut,
                  double r_buf_ratio = 0, int nPar = 0);

  bool has_member(int ic) const { return ! cell[ic].empty(); }

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij) const;

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij, int ic) const;

  template <class BiFunc>
  void for_each_pair(BiFunc f_ij, int ic1, int ic2) const;

  template <class Par>
  void create(const std::vector<Par> &p_arr);

  template <class Par>
  void recreate(const std::vector<Par> &p_arr);

  template <class Par>
  void update(const std::vector<Par> &p_arr);

  template <class Par, class BiFunc>
  void cal_force(const std::vector<Par> &p_arr, BiFunc f_ij);

#ifdef SPATIAL_SORT
  template <class Par, class BiFunc>
  void cal_force(std::vector<Par>, BiFunc pair_force,
                 MySpatialSortingTraits<Par> &sst);
#endif


protected:
  std::vector< std::list<int> > cell;
  std::vector<int> list_len;
  int count;
};


template<class BiFunc>
void CellList_list_2::for_each_pair(BiFunc f_ij) const {
  auto lambda = [this, f_ij](int *ic) {
    if (cell[ic[0]].size() > 1)
      for_each_pair(f_ij, ic[0]);
    if (!cell[ic[1]].empty())
      for_each_pair(f_ij, ic[0], ic[1]);
    if (!cell[ic[2]].empty())
      for_each_pair(f_ij, ic[0], ic[2]);
    if (!cell[ic[3]].empty())
      for_each_pair(f_ij, ic[0], ic[3]);
    if (!cell[ic[4]].empty())
      for_each_pair(f_ij, ic[0], ic[4]);
  };
  for_each_nearest_cell(lambda);
}

template<class BiFunc>
void CellList_list_2::for_each_pair(BiFunc f_ij, int ic) const {
  auto end2 = cell[ic].cend();
  auto end1 = std::prev(end2);
  for (auto it1 = cell[ic].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = std::next(it1); it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template<class BiFunc>
void CellList_list_2::for_each_pair(BiFunc f_ij, int ic1, int ic2) const {
  auto end1 = cell[ic1].cend();
  auto end2 = cell[ic2].cend();
  auto beg2 = cell[ic2].cbegin();
  for (auto it1 = cell[ic1].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = beg2; it2 != end2; ++it2) {
      f_ij(*it1, *it2);
    }
  }
}

template<class Par>
void CellList_list_2::create(const std::vector<Par>& p_arr) {
  int nPar = p_arr.size();
  for (int ip = 0; ip < nPar; ip++) {
    cell[get_ic(p_arr[ip])].push_back(ip);
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class Par>
void CellList_list_2::recreate(const std::vector<Par>& p_arr) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic].clear();
  }
  create(p_arr);
}

template<class Par>
void CellList_list_2::update(const std::vector<Par>& p_arr) {
  for (int iy = 0; iy < bins.y; iy++) {
    int iy_times_nx = iy * bins.x;

    for (int ix = 0; ix < bins.x; ix++) {
      int ic = ix + iy_times_nx;
      int depth = list_len[ic];
      if (depth) {
        int it_count = 0;
        for (auto it = cell[ic].begin(); it_count < depth; it_count++) {
          int ip = *it;
          int ic_new = get_ic(p_arr[ip]);
          if (ic_new == ic) {
            ++it;
          } else {
            it = cell[ic].erase(it);
            cell[ic_new].push_back(ip);
          }
        }
      }
    }
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class Par, class BiFunc>
inline void CellList_list_2::cal_force(const std::vector<Par>& p_arr,
                                       BiFunc f_ij) {
  update(p_arr);
  for_each_pair(f_ij);
}

template<class Par, class BiFunc>
inline void CellList_list_2::cal_force(std::vector<Par> p_arr, BiFunc pair_force,
                                       MySpatialSortingTraits<Par>& sst) {
  if (count % 1000 == 0) {
    CGAL::spatial_sort(p_arr.begin(), p_arr.end(), sst);
    recreate(p_arr);
  } else {
    update(p_arr);
  }
  for_each_pair(pair_force);
}

#endif
