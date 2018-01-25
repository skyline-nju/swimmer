#ifndef CELL_LIST_H
#define CELL_LIST_H
#include <vector>
#include <list>
#include <iostream>
#include <functional>

template <class Par>
class NodeWrapper :public Par {
public:
  NodeWrapper() : Par(), next(nullptr), cell_idx(-1) {}

  int get_cell_idx() { return cell_idx; }
  bool update_cell_idx(double lx, double ly, int ncols);

  template <class BiFunc>
  void for_each_pair(BiFunc f2);

  template <class BiFunc>
  void for_each_pair(NodeWrapper<Par> *head2, BiFunc f2);

  NodeWrapper * next;
  int cell_idx;

};

template <class Par>
inline bool NodeWrapper<Par>::update_cell_idx(double lx, double ly, int ncols) {
  int idx_new = int(this->x / lx) + int(this->y / ly) * ncols;
  if (idx_new == cell_idx) {
    return false;
  } else {
    cell_idx = idx_new;
    return true;
  }
}

template<class Par>
template<class BiFunc>
void NodeWrapper<Par>::for_each_pair(BiFunc f2) {
  NodeWrapper<Par> *node1 = this;
  NodeWrapper<Par> *node2;
  while (node1->next) {
    node2 = node1->next;
    do {
      f2(*node1, *node2);
      node2 = node2->next;
    } while (node2);
    node1 = node1->next;
  }
}

template<class Par>
template<class BiFunc>
void NodeWrapper<Par>::for_each_pair(NodeWrapper<Par>* head2, BiFunc f2) {
  if (head2) {
    NodeWrapper<Par> *node1 = this;
    NodeWrapper<Par> *node2;
    do {
      node2 = head2;
      do {
        f2(*node1, *node2);
        node2 = node2->next;
      } while (node2);
      node1 = node1->next;
    } while (node1);
  }
}

/*************************************************************************/
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

/*************************************************************************/
template <class Node>
class CellList_w_Node: public CellListBase {
public:
  CellList_w_Node(double Lx0, double Ly0, double rcut0);
  ~CellList_w_Node() {};

  template<class BiFunc>
  void for_each_pair(BiFunc f2) const;

  void link_nodes(Node *p, int nPar);

  void refresh(Node *p, int nPar);

  void update(Node *p, int nPar) { refresh(p, nPar); }

  bool has_member(int ic) const { return nullptr != cell[ic]; }

protected:
  std::vector<Node *> cell;
};

template<class Node>
CellList_w_Node<Node>::CellList_w_Node(double Lx0, double Ly0, double rcut0) :
                                       CellListBase(Lx0, Ly0, rcut0) {
  cell.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.push_back(nullptr);
  }
  std::cout << "initialize cell list with " << ncols << " columns by "
    << nrows << " rows" << std::endl;
  std::cout << "cell size is " << lx << " by " << ly << std::endl;
}

template<class Node>
template<class BiFunc>
void CellList_w_Node<Node>::for_each_pair(BiFunc f2) const {
  auto f1 = [this, f2](int *idx) {
    int i = idx[0];
    cell[i]->for_each_pair(f2);
    cell[i]->for_each_pair(cell[idx[1]], f2);
    cell[i]->for_each_pair(cell[idx[2]], f2);
    cell[i]->for_each_pair(cell[idx[3]], f2);
    cell[i]->for_each_pair(cell[idx[4]], f2);
  };
  for_each_nearest_cell(f1);
}

template<class Node>
void CellList_w_Node<Node>::link_nodes(Node *p, int nPar) {
  for (int i = 0; i < nPar; i++) {
    p[i].update_cell_idx(lx, ly, ncols);
    int i_cell = p[i].get_cell_idx();
    p[i].next = cell[i_cell];
    cell[i_cell] = &p[i];
  }
}

template<class Node>
void CellList_w_Node<Node>::refresh(Node *p, int nPar) {
  for (int i = 0; i < ncells; i++) {
    cell[i] = nullptr;
  }
  link_nodes(p, nPar);
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
  void create(const Par *par, int nPar);

  template <class Par>
  void recreate(const Par *p, int nPar);

  template <class Par>
  void update(const Par *p, int nPar);

  bool has_member(int ic) const { return cell[ic].size() > 0; }

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
void CellList_w_list::create(const Par * par, int nPar) {
  for (int ip = 0; ip < nPar; ip++) {
    int idx_cell = int(par[ip].x / lx) + ncols * int(par[ip].y / ly);
    cell[idx_cell].push_back(ip);
  }
  for (int ic = 0; ic < ncells; ic++) {
    list_len[ic] = cell[ic].size();
  }
}

template<class Par>
void CellList_w_list::recreate(const Par * p, int nPar) {
  for (int ic = 0; ic < ncells; ic++) {
    cell[ic].clear();
  }
  create(p, nPar);
}

template<class Par>
void CellList_w_list::update(const Par * p, int nPar) {
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
