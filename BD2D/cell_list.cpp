#include "cell_list.h"

CellListBase::CellListBase(double Lx0, double Ly0, double rcut0):
                           Lx(Lx0), Ly(Ly0) {
  ncols = int(Lx / rcut0);
  nrows = int(Ly / rcut0);
  lx = Lx / ncols;
  ly = Ly / nrows;
  ncells = ncols * nrows;

}
CellList_w_list::CellList_w_list(double Lx0, double Ly0, double rcut0): 
                                 CellListBase(Lx0, Ly0, rcut0) {
  cell.reserve(ncells);
  list_len.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.push_back(std::list<int>());
    list_len.push_back(0);
  }
}
