#include "cellList.h"

CellListBase_2::CellListBase_2(double Lx, double Ly, double rcut) {
  bins.x = int(Lx / rcut);
  bins.y = int(Ly / rcut);
  l_cell.x = Lx / bins.x;
  l_cell.y = Ly / bins.y;
  ncells = bins.x * bins.y;
}

CellList_list_2::CellList_list_2(double Lx, double Ly, double rcut,
  double r_buf_ratio, int nPar): CellListBase_2(Lx, Ly, rcut) {
  cell.reserve(ncells);
  list_len.reserve(ncells);
  for (int i = 0; i < ncells; i++) {
    cell.emplace_back();
    list_len.push_back(0);
  }
}
