#include "cellList3D.h"

const int CellListBase_3::cell_offset[13][3] = {
  { 0,  0,  1 },
  { 0,  1,  0 },
  { 0,  1,  1 },
  { 1,  0,  0 },
  { 1,  0,  1 },
  { 1,  1,  0 },
  { 1,  1,  1 },
  { 0,  1, -1 },
  { 1,  0, -1 },
  { 1,  1, -1 },
  { 1, -1,  0 },
  { 1, -1,  1 },
  { 1, -1, -1 }
};


CellListBase_3::CellListBase_3(const Vec3d & l,
                               double r_cut,
                               const Vec3d &gl_l,
                               const Vec3d & origin,
                               const Vec_3<bool>& flag_ext)
  : n_(), origin_(origin), inverse_lc_(), l_(l), gl_l_(gl_l),
    flag_ext_(flag_ext) {
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
  : n_(cell_size), origin_(origin), inverse_lc_(),
    l_(lc * n_), gl_l_(gl_l), flag_ext_(flag_ext) {
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

void CellListBase_3::partition(const Vec_3<double>& l, double r_cut,
                               Vec_3<int>& cells_size, Vec_3<double>& cell_len) {
  for (int dim = 0; dim < 3; dim++) {
    cells_size[dim] = l[dim] / r_cut;
    cell_len[dim] = l[dim] / cells_size[dim];
  }
}

/**
 * @brief Get a vector to offset the periodic boundary condition
 * 
 * @param pos             Position of a particle
 * @return Vec_3<double>  
 */
Vec_3<double> CellListBase_3::get_offset(const Vec3d& pos) const {
  Vec_3<double> offset{};
  Vec_3<double> dR = pos - origin_;

  //! If canceling the annotation of the following line, the function can give
  //! right results.
  // std::cout << dR << std::endl;

  for (int dim = 0; dim < 3; dim++) {
    if (dR[dim] < 0) {
      offset[dim] = gl_l_[dim];
    } else if (dR[dim] > l_[dim]) {
      offset[dim] = -gl_l_[dim];
    }
  }
  return offset;
}
