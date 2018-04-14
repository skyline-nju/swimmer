#include "cellList2D.h"

Cell_list_base_2::Cell_list_base_2(const Vec_2<double>& l_domain,
                                   const Vec_2<double>& origin,
                                   double r_cut,
                                   const Vec_2<bool>& flag_comm):
                                   flag_comm_(flag_comm), r_cut_(r_cut) {
  if (!flag_comm.x && !flag_comm.y) {
    origin_ = origin;
    n_.x = int(l_domain.x / r_cut);
    n_.y = int(l_domain.y / r_cut);
    l_cell_.x = l_domain.x / n_.x;
    l_cell_.y = l_domain.y / n_.y;
    l_domain_ = l_domain;
  } else if (!flag_comm.x && flag_comm.y) {
    n_.x = int(l_domain.x / r_cut);
    n_.y = int(l_domain.y / r_cut) + 2;
    l_cell_.x = l_domain.x / n_.x;
    l_cell_.y = l_domain.y / (n_.y - 2);
    l_domain_.x = l_domain.x;
    l_domain_.y = l_domain.y + 2 * l_cell_.y;
    origin_.x = origin.x;
    origin_.y = origin.y - l_cell_.y;
  } else if (flag_comm.x && flag_comm.y) {
    n_.x = int(l_domain.x / r_cut) + 2;
    n_.y = int(l_domain.y / r_cut) + 2;
    l_cell_.x = l_domain.x / (n_.x - 2);
    l_cell_.y = l_domain.y / (n_.y - 2);
    l_domain_ = l_domain + 2 * l_cell_;
    origin_ = origin - l_cell_;
  }
  inverse_l_cell_ = l_cell_.inverse();
  ncells_ = n_.x * n_.y;
}

int Cell_list_base_2::get_row(int my_row, int drow) const {
  auto row = my_row + drow;
  if (row < 0) {
    row += n_.y;
  } else if (row >= n_.y) {
    row -= n_.y;
  }
  return row;
}

int Cell_list_base_2::get_col(int my_col, int dcol) const {
  auto col = my_col + dcol;
  if (col < 0) {
    col += n_.x;  
  } else if (col >= n_.x) {
    col -= n_.x;
  }
  return col;
}

