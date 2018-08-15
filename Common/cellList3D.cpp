#include "cellList3D.h"

void divide_cell(const Vec_3<double>& l, double r_cut,
                 Vec_3<int>& cell_num, Vec_3<double>& cell_len) {
  for (int dim = 0; dim < 3; dim++) {
    cell_num[dim] = l[dim] / r_cut;
    cell_len[dim] = l[dim] / cell_num[dim];
  }
}

const int CellListBase_3::cell_offset[13][3] = {
  { 0, 0, 1 },
  { 0, 1, 0 },
  { 0, 1, 1 },
  { 1, 0, 0 },
  { 1, 0, 1 },
  { 1, 1, 0 },
  { 1, 1, 1 },
  { 0, 1, -1 },
  { 1, 0, -1 },
  { 1, 1, -1 },
  { 1, -1, 0 },
  { 1, -1, 1 },
  { 1, -1, -1 }
};


CellListBase_3::CellListBase_3(const Vec3d & l,
                               double r_cut,
                               const Vec3d &gl_l,
                               const Vec3d & origin,
                               const Vec_3<bool>& flag_ext)
  : origin_(origin), l_(l), gl_l_(gl_l), flag_ext_(flag_ext) {
  for (int dim = 0; dim < 3; dim++) {
    n_[dim] = int(l[dim] / r_cut);
    inverse_lc_[dim] = n_[dim] / l[dim];
    if(flag_ext[dim]) {
      origin_[dim] -= l[dim] / n_[dim];
      n_[dim] += 2;
      l_[dim] += 2 * l[dim] / n_[dim];
    }
  }
  ncells_ = n_.x * n_.y * n_.z;
  nxny_ = n_.x * n_.y;
}

CellListBase_3::CellListBase_3(const Vec_3<int> &cell_size,
                               const Vec3d& lc,
                               const Vec3d &gl_l,
                               const Vec3d& origin,
                               const Vec_3<bool>& flag_ext)
  : n_(cell_size), origin_(origin), l_(lc * n_), gl_l_(gl_l), flag_ext_(flag_ext) {
  for (int dim = 0; dim < 3; dim++) {
    inverse_lc_[dim] = 1 / lc[dim];
    if(flag_ext[dim]) {
      origin_[dim] -= lc[dim];
      n_[dim] += 2;
      l_[dim] += 2. * lc[dim];
    }
  }
  ncells_ = n_.x * n_.y * n_.z;
  nxny_ = n_.x * n_.y;
}

void CellListBase_3::show_para() const{
  std::cout << "cell origin: " << origin_ << std::endl;
  std::cout << "cell size: " << n_ << std::endl;
}

Vec_3<double> CellListBase_3::get_offset(const Vec3d& pos) {
  Vec_3<double> offset{};
  Vec_3<double> dR = pos - origin_;
  for (int dim = 0; dim < 3; dim++) {
    if (dR[dim] > l_[dim])
      offset[dim] = -gl_l_[dim];
    else if (dR[dim] < -l_[dim])
      offset[dim] = gl_l_[dim];
    else
      offset[dim] = 0;
  }
  return offset;
}
