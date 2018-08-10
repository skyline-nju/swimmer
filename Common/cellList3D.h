#pragma once

#include <vector>
#include "vect.h"
#include "node.h"

class CellListBase_3 {
public:
  typedef Vec_3<double> Vec3d;
  CellListBase_3(const Vec3d &l, double r_cut, const Vec3d &origin);
  
  int get_nx(double x) const {
    return int((x - origin_.x) * inverse_lc_.x);
  }
  int get_ny(double y) const {
    return int((y - origin_.y) * inverse_lc_.y);
  }
  int get_nz(double z) const {
    return int((z - origin_.z) * inverse_lc_.z);
  }
  template <typename TPar>
  int get_ic(const TPar &p) const {
    return get_nx(p.pos.x) + get_ny(p.pos.y) * n_.x + get_nz(p.pos.z) * nxny_;
  }

protected:
  int ncells_;
  int nxny_;
  Vec_3<int> n_;
  Vec_3<double> origin_;
  Vec_3<double> inverse_lc_;
};

inline CellListBase_3::CellListBase_3(const Vec3d &l, double r_cut, const Vec3d &origin)
  : origin_(origin) {
  n_.x = int(l.x / r_cut);
  n_.y = int(l.y / r_cut);
  n_.z = int(l.z / r_cut);
  inverse_lc_.x = n_.x / l.x;
  inverse_lc_.y = n_.y / l.y;
  inverse_lc_.z = n_.z / l.z;
  ncells_ = n_.x * n_.y * n_.z;
  nxny_ = n_.x * n_.y;
}

template <typename TNode>
class CellListNode_3: public CellListBase_3 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  
  CellListNode_3(const Vec3d &l, double r_cut, const Vec3d &origin = Vec_3<double>())
    : CellListBase_3(l, r_cut, origin), head_(ncells_) {}

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij, const Vec_3<int> &ic_beg, const Vec_3<int> &ic_end);

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij);

  void create(std::vector<TNode> &p_arr);

  void recreate(std::vector<TNode> &p_arr);

protected:
  std::vector<TNode*> head_;
};

template <typename TNode>
template <typename BiFunc>
void CellListNode_3<TNode>::for_each_pair(BiFunc f_ij,
                                          const Vec_3<int>& ic_beg,
                                          const Vec_3<int>& ic_end) {
  for(int z0 = ic_beg.z; z0 < ic_end.z; z0++) {
    auto z1 = z0 + 1;
    if (z1 >= n_.z)
      z1 -= n_.z;
    const auto z0_nxny = z0 * nxny_;
    const auto z1_nxny = z1 * nxny_;
    for (int y0 = ic_beg.y; y0 < ic_end.y; y0++) {
      auto y1 = y0 + 1;
      if (y1 >= n_.y)
        y1 -= n_.y;
      const auto y0_nx = y0 * n_.x;
      const auto y1_nx = y1 * n_.x;
      for (int x0 = ic_beg.x; x0 < ic_end.x; x0++) {
        auto x1 = x0 + 1;
        if (x1 >= n_.x)
          x1 -= n_.x;
        int ic[2][2][2] = {
          x0 + y0_nx + z0_nxny,
          x1 + y0_nx + z0_nxny,
          x0 + y1_nx + z0_nxny,
          x1 + y1_nx + z0_nxny,
          x0 + y0_nx + z1_nxny,
          x1 + y0_nx + z1_nxny,
          x0 + y1_nx + z1_nxny,
          x1 + y1_nx + z1_nxny
        };
        if (head_[ic[0][0][0]]) {
          for_each_node_pair(head_[ic[0][0][0]], f_ij);

          if (head_[ic[0][0][1]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[0][0][1]], f_ij);
            if (head_[ic[0][1][0]]) {
              for_each_node_pair(head_[ic[0][0][1]], head_[ic[0][1][0]], f_ij);
            }
            if (head_[ic[1][0][0]]) {
              for_each_node_pair(head_[ic[0][0][1]], head_[ic[1][0][0]], f_ij);
            }
            if (head_[ic[1][1][0]]) {
              for_each_node_pair(head_[ic[0][0][1]], head_[ic[1][1][0]], f_ij);
            }
          }

          if (head_[ic[0][1][0]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[0][1][0]], f_ij);
            if (head_[ic[1][0][0]]) {
              for_each_node_pair(head_[ic[0][1][0]], head_[ic[1][0][0]], f_ij);
            }
            if (head_[ic[1][0][1]]) {
              for_each_node_pair(head_[ic[0][1][0]], head_[ic[1][0][1]], f_ij);
            }
          }

          if (head_[ic[0][1][1]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[0][1][1]], f_ij);
            if (head_[ic[1][0][0]]) {
              for_each_node_pair(head_[ic[0][1][1]], head_[ic[1][0][0]], f_ij);
            }
          }

          if (head_[ic[1][0][0]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[1][0][0]], f_ij);
          }

          if (head_[ic[1][0][1]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[1][0][1]], f_ij);
          }

          if (head_[ic[1][1][0]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[1][1][0]], f_ij);
          }

          if (head_[ic[1][1][1]]) {
            for_each_node_pair(head_[ic[0][0][0]], head_[ic[1][1][1]], f_ij);
          }
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc>
void CellListNode_3<TNode>::for_each_pair(BiFunc f_ij) {
  for_each_pair(f_ij, Vec_3<int>(), n_);
}

template<typename TNode>
void CellListNode_3<TNode>::create(std::vector<TNode>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head_[ic]);
  }
}

template <typename TNode>
void CellListNode_3<TNode>::recreate(std::vector<TNode>& p_arr) {
  for (int ic=0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
  }
  create(p_arr);
}
