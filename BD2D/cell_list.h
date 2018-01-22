#ifndef CELL_LIST_H
#define CELL_LIST_H

template <class Par>
class NodeWrapper :public Par {
public:
  NodeWrapper() : Par(), next(nullptr) {}


  NodeWrapper * next;
  int cell_idx;
};

template <class Par>
class Cell_list_2 {
public:
  Cell_list_2(double Lx0, double Ly0, double lx0, double ly0);
  ~Cell_list_2();
  NodeWrapper<Par> *cell;
  
protected:
  double Lx;
  double Ly;
  int nrows;
  int ncols;
  int ncells;
};


template<class Par>
Cell_list_2<Par>::Cell_list_2(double Lx0, double Ly0, double lx0, double ly0):
  Lx(Lx0), Ly(Ly0) {
  nrows = int(Lx / lx0);
  ncols = int(Ly / ly0);
  ncells = nrows * ncols;
  cell = new NodeWrapper<Par>[ncells];
  std::cout << "initialize cell list with " << ncols << " columns by "
    << nrows << " rows" << std::endl;
}

template<class Par>
Cell_list_2<Par>::~Cell_list_2() {
  delete[] cell;
}

#endif