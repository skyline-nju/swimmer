#ifndef CELL_LIST_H
#define CELL_LIST_H
#include <vector>
#include <list>

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
  int idx_new = int(x / lx) + int(y / ly) * ncols;
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

template <class Node>
class Cell_list_2 {
public:
  Cell_list_2(double Lx0, double Ly0, double lx0, double ly0);
  ~Cell_list_2() {};

  template<class BiFunc>
  void for_each_pair(BiFunc f2) const;

  void link_nodes(Node *p, int nPar);

  void refresh(Node *p, int nPar);

protected:
  std::vector<Node *> cell;
  double Lx;
  double Ly;
  double lx;
  double ly;
  int nrows;
  int ncols;
  int ncells;
};

template<class Node>
Cell_list_2<Node>::Cell_list_2(double Lx0, double Ly0, double lx0, double ly0):
  Lx(Lx0), Ly(Ly0) {
  ncols = int(Lx / lx0);
  nrows = int(Ly / ly0);
  lx = Lx / ncols;
  ly = Ly / nrows;
  ncells = nrows * ncols;
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
void Cell_list_2<Node>::for_each_pair(BiFunc f2) const {
  for (int row = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      int i = col + row * ncols;
      if (cell[i]) {
        int col_right = col + 1;
        if (col_right >= ncols)
          col_right = 0;
        int col_left = col - 1;
        if (col_left < 0)
          col_left = ncols - 1;
        int row_upper = row + 1;
        if (row_upper >= nrows)
          row_upper = 0;      
        // self --> right --> upper left --> upper -->upper right
        cell[i]->for_each_pair(f2);
        cell[i]->for_each_pair(cell[col_right + row * ncols], f2);
        cell[i]->for_each_pair(cell[col_left + row_upper * ncols], f2);
        cell[i]->for_each_pair(cell[col + row_upper * ncols], f2);
        cell[i]->for_each_pair(cell[col_right + row_upper * ncols], f2);
      }
    }
  }
}

template<class Node>
void Cell_list_2<Node>::link_nodes(Node *p, int nPar) {
  for (int i = 0; i < nPar; i++) {
    p[i].update_cell_idx(lx, ly, ncols);
    int i_cell = p[i].get_cell_idx();
    p[i].next = cell[i_cell];
    cell[i_cell] = &p[i];
  }
}

template<class Node>
void Cell_list_2<Node>::refresh(Node *p, int nPar) {
  for (int i = 0; i < ncells; i++) {
    cell[i] = nullptr;
  }
  link_nodes(p, nPar);
}

#endif

