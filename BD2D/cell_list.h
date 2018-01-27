#ifndef CELL_LIST_H
#define CELL_LIST_H
#include <vector>
#include <list>
#include <iostream>
#include <functional>

class CellListBase {
public:
  CellListBase(double Lx0, double Ly0, double rcut0);

  void PBCx(int &col) const;

  void PBCy(int &row) const;

  template <class Par>
  int check_pos(int col0, int row0, const Par &p) const;

  virtual bool has_member(int ic) const = 0;

  template <class UniFunc>
  void for_each_nearest_cell(UniFunc f) const;

  void set_r_buf(double r_buf, int nPar) {}

protected:
  double Lx;
  double Ly;
  double lx;
  double ly;
  int nrows;
  int ncols;
  int ncells;

};

inline void CellListBase::PBCx(int &col) const {
  if (col < 0)
    col += ncols;
  else if (col >= ncols)
    col -= ncols;
}

inline void CellListBase::PBCy(int &row) const {
  if (row < 0)
    row += nrows;
  else if (row >= nrows)
    row -= nrows;
}

template<class Par>
int CellListBase::check_pos(int col0, int row0, const Par & p) const {
  int col = int(p.x / lx);
  int row = int(p.y / ly);
  if (col == col0) {
    if (row == row0) {
      return -1;          // the position has not been changed.
    } else {
      PBCy(row);
      return col + ncols * row;
    }
  } else {
    PBCx(col);
    if (row != row0) {
      PBCy(row);
    }
    return col + ncols * row;
  }
}

template <class UniFunc>
void CellListBase::for_each_nearest_cell(UniFunc f) const {
  int idx_cell[5];
  for (int row = 0; row < nrows; row++) {
    int row_upper = row + 1;
    if (row_upper >= nrows)
      row_upper = 0;
    int row_times_ncols = row * ncols;
    int row_upper_times_ncols = row_upper * ncols;

    for (int col = 0; col < ncols; col++) {
      int i = col + row_times_ncols;
      if (has_member(i)) {
        int col_right = col + 1;
        if (col_right >= ncols)
          col_right = 0;
        int col_left = col - 1;
        if (col_left < 0)
          col_left = ncols - 1;

        idx_cell[0] = i;
        idx_cell[1] = col_right + row_times_ncols;
        idx_cell[2] = col_left + row_upper_times_ncols;
        idx_cell[3] = col + row_upper_times_ncols;
        idx_cell[4] = col_right + row_upper_times_ncols;
        f(idx_cell);
      }
    }
  }
}

/*****************************************************************************/

class CellList_w_list: public CellListBase{
public:
  CellList_w_list(double Lx0, double Ly0, double rcut0);
  
  template <class BiFunc>
  void for_each_pair(BiFunc f2) const;

  template <class BiFunc>
  void for_each_pair(BiFunc f2, int ic) const;

  template <class BiFunc>
  void for_each_pair(BiFunc f2, int ic1, int ic2) const;

  template <class Par>
  void create(const std::vector<Par> &p);

  template <class Par>
  void recreate(const std::vector<Par> &p);

  template <class Par>
  void update(const std::vector<Par> &par);

  bool has_member(int ic) const { return cell[ic].size() > 0; }

  template <class Par, class BiFunc>
  void cal_force(std::vector<Par> &par, BiFunc pair_force) {
    update(par, pair_force);
    for_each_pair(pair_force);
  }

protected:
  std::vector<std::list <int> > cell;
  std::vector<int> list_len;
};

template <class BiFunc>
void CellList_w_list::for_each_pair(BiFunc f2) const {
  auto f1 = [this, f2](int *idx) {
    for_each_pair(f2, idx[0]);
    if (cell[idx[1]].size() > 0)
      for_each_pair(f2, idx[0], idx[1]);
    if (cell[idx[2]].size() > 0)
      for_each_pair(f2, idx[0], idx[2]);
    if (cell[idx[3]].size() > 0)
      for_each_pair(f2, idx[0], idx[3]);
    if (cell[idx[4]].size() > 0)
      for_each_pair(f2, idx[0], idx[4]);
  };
  for_each_nearest_cell(f1);
}

template<class BiFunc>
void CellList_w_list::for_each_pair(BiFunc f2, int ic) const {
  if (cell[ic].size() > 1) {
    auto end2 = cell[ic].cend();
    auto end1 = std::prev(end2);
    for (auto it1 = cell[ic].cbegin(); it1 != end1; ++it1) {
      for (auto it2 = std::next(it1); it2 != end2; ++it2) {
        f2(*it1, *it2);
      }
    }
  }
}

template<class BiFunc>
void CellList_w_list::for_each_pair(BiFunc f2, int ic1, int ic2) const {
  auto end1 = cell[ic1].cend();
  auto end2 = cell[ic2].cend();
  auto beg2 = cell[ic2].cbegin();
  for (auto it1 = cell[ic1].cbegin(); it1 != end1; ++it1) {
    for (auto it2 = beg2; it2 != end2; ++it2) {
      f2(*it1, *it2);
    }
  }
}

template<class Par>
void CellList_w_list::create(const std::vector<Par> &par) {
  int nPar = par.size();
  for (int ip = 0; ip < nPar; ip++) {
    int idx_cell = int(par[ip].x / lx) + ncols * int(par[ip].y / ly);
    cell[idx_cell].push_back(ip);
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class Par>
void CellList_w_list::recreate(const std::vector<Par> &p) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic].clear();
  }
  create(p);
}

template<class Par>
void CellList_w_list::update(const std::vector<Par> &p) {
  for (int row = 0; row < nrows; row++) {
    int row_times_ncols = row * ncols;

    for (int col = 0; col < ncols; col++) {
      int ic = col + row_times_ncols;
      int depth = list_len[ic];
      if (depth == 1) {
        int ip = cell[ic].front();
        int ic_new = check_pos(col, row, p[ip]);
        if (ic_new >= 0) {
          cell[ic].erase(cell[ic].begin());
          cell[ic_new].push_back(ip);
        }
      } else if (depth > 1) {
        int it_count = 0;
        for (auto it = cell[ic].begin(); it_count < depth; it_count++) {
          int ip = *it;
          int ic_new = check_pos(col, row, p[ip]);
          if (ic_new < 0) {
            ++it;
          } else {
            it = cell[ic].erase(it);
            cell[ic_new].push_back(ip);
          }
        }
      }
    }
  }
  for (int i = 0; i < ncells; i++) {
    list_len[i] = cell[i].size();
  }
}


#endif
