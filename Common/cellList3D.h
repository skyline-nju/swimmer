/**
 * @brief 3d linked-cell list
 * 
 * @file cellList3D.h
 * @author your name
 * @date 2018-08-17
 */
#pragma once
#include <vector>
#include "vect.h"
#include "node.h"

class CellListBase_3 {
public:
  typedef Vec_3<double> Vec3d;

  CellListBase_3(const Vec3d &l,
                 double r_cut,
                 const Vec3d &gl_l,
                 const Vec3d &origin,
                 const Vec_3<bool> &flag_ext);

  CellListBase_3(const Vec_3<int> &cell_size,
                 const Vec3d &lc,
                 const Vec3d &gl_l,
                 const Vec3d &origin,
                 const Vec_3<bool> &flag_ext);
  
  static void partition(const Vec_3<double> &l, double r_cut,
                        Vec_3<int> &cells_size, Vec_3<double> &cell_len);

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

  void show_para() const;

  int ncells() const { return ncells_; }
  const Vec3d& origin() const { return origin_; }
  const Vec3d& l() const { return l_; }
  const Vec_3<bool>& flag_ext() const { return flag_ext_; }
  const Vec_3<int>& cells_size() const { return n_; }

  Vec_3<double> get_offset(const Vec3d &pos);
protected:
  int ncells_;
  int nxny_;
  Vec_3<int> n_;
  Vec_3<double> origin_;
  Vec_3<double> inverse_lc_;
  Vec_3<double> l_;
  Vec_3<double> gl_l_;
  Vec_3<bool> flag_ext_;

  const static int cell_offset[13][3];
};

template <typename TNode>
class CellListNode_3: public CellListBase_3 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  typedef Vec_3<int> Vec3i;
  
  CellListNode_3(const Vec3d &l, double r_cut,
                 const Vec3d &gl_l,
                 const Vec3d &origin = Vec3d(),
                 const Vec_3<bool> &flag_ext = Vec_3<bool>())
                 : CellListBase_3(l, r_cut, gl_l, origin, flag_ext), head_(ncells_) {}

  CellListNode_3(const Vec_3<int> &cell_size,
                 const Vec3d &lc,
                 const Vec3d &gl_l,
                 const Vec3d &origin = Vec3d(),
                 const Vec_3<bool> &flag_ext = Vec_3<bool>())
                 : CellListBase_3(cell_size, lc, gl_l, origin, flag_ext), head_(ncells_) {}


  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2, const Vec3i &ic_beg, const Vec3i &ic_end) const;

  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2, const Vec3i &ic_beg, const Vec3i &ic_end) const;

  template <typename BiFunc>
  void for_each_pair_slow(BiFunc f_ij, const Vec3i &ic_beg, const Vec3i &ic_end) const;

  template <typename BiFunc1, typename BiFunc2>
  void for_each_pair(BiFunc1 f1, BiFunc2 f2) const;

  template <typename BiFunc, typename TriFunc>
  void for_each_pair_fast(BiFunc f1, TriFunc f2) const;

  void create(std::vector<TNode> &p_arr);

  template <typename T>
  void create(std::vector<TNode> &p_arr, T* par_num_arr);

  void recreate(std::vector<TNode> &p_arr);

  template <typename T>
  void recreate(std::vector<TNode> &p_arr, T* par_num_arr);

  template <typename UniFunc>
  void for_each_cell(UniFunc f, const Vec_3<int> &beg, const Vec_3<int> &end);

  void add_node(TNode &p) {
    const int idx = get_ic(p);
    p.append_at_front(&head_[idx]);
  }

  void clear(const Vec_3<int> &first, const Vec_3<int> &last);
  
  template <typename TPar>
  void replace(BiNode<TPar> &p_new, const BiNode<TPar> &p_old);

  int get_par_num(const Vec_3<int> &beg, const Vec_3<int> &end) const;

  int get_par_num() const;
protected:
  std::vector<TNode*> head_;
};

template <typename TNode>
template <typename BiFunc>
void CellListNode_3<TNode>::for_each_pair_slow(BiFunc f_ij,
                                               const Vec3i& ic_beg,
                                               const Vec3i& ic_end) const {
  for (int z0 = ic_beg.z; z0 < ic_end.z; z0++) {
    for (int y0 = ic_beg.y; y0 < ic_end.y; y0++) {
      for (int x0 = ic_beg.x; x0 < ic_end.x; x0++) {
        int i0 = x0 + y0 * n_.x + z0 * nxny_;
        if (head_[i0]) {
          for_each_node_pair(head_[i0], f_ij);
          for (int j = 0; j < 13; j++) {
            int z1 = z0 + cell_offset[j][0];
            int y1 = y0 + cell_offset[j][1];
            int x1 = x0 + cell_offset[j][2];
            if (z1 >= n_.z) {
              z1 = 0;
            }
            if (y1 >= n_.y) {
              y1 = 0;
            } else if (y1 < 0) {
              y1 += n_.y;
            }
            if (x1 >= n_.x) {
              x1 = 0;
            } else if (x1 < 0) {
              x1 += n_.x;
            }
            int i1 = x1 + y1 * n_.x + z1 * nxny_;
            if (head_[i1]) {
              for_each_node_pair(head_[i0], head_[i1], f_ij);
            }
          }
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_3<TNode>::for_each_pair(BiFunc1 f1, BiFunc2 f2, const Vec3i& ic_beg, const Vec3i& ic_end) const {
  int z[2];
  int y[2];
  int x[2];
  int i[8];
  for (z[0] = ic_beg.z; z[0] < ic_end.z; z[0]++) {
    z[1] = z[0] + 1;
    if (z[1] >= n_.z)
      z[1] -= n_.z;
    const int z_nx_ny[2] = { z[0] * nxny_, z[1] * nxny_ };
    for (y[0] = ic_beg.y; y[0] < ic_end.y; y[0]++) {
      y[1] = y[0] + 1;
      if (y[1] >= n_.y)
        y[1] -= n_.y;
      const int y_nx[2] = { y[0] * n_.x, y[1] * n_.y };
      for (x[0] = ic_beg.x; x[0] < ic_end.x; x[0]++) {
        x[1] = x[0] + 1;
        if (x[1] >= n_.x)
          x[1] -= n_.x;
        i[0] = x[0] + y_nx[0] + z_nx_ny[0];
        i[1] = x[1] + y_nx[0] + z_nx_ny[0];
        i[2] = x[0] + y_nx[1] + z_nx_ny[0];
        i[3] = x[1] + y_nx[1] + z_nx_ny[0];
        i[4] = x[0] + y_nx[0] + z_nx_ny[1];
        i[5] = x[1] + y_nx[0] + z_nx_ny[1];
        i[6] = x[0] + y_nx[1] + z_nx_ny[1];
        i[7] = x[1] + y_nx[1] + z_nx_ny[1];

        if (head_[i[0]]) {
          for_each_node_pair(head_[i[0]], f1);
          if (head_[i[1]]) {
            for_each_node_pair(head_[i[0]], head_[i[1]], f2);
          }
          if (head_[i[2]]) {
            for_each_node_pair(head_[i[0]], head_[i[2]], f2);
          }
          if (head_[i[3]]) {
            for_each_node_pair(head_[i[0]], head_[i[3]], f2);
          }
          if (head_[i[4]]) {
            for_each_node_pair(head_[i[0]], head_[i[4]], f2);
          }
          if (head_[i[5]]) {
            for_each_node_pair(head_[i[0]], head_[i[5]], f2);
          }
          if (head_[i[6]]) {
            for_each_node_pair(head_[i[0]], head_[i[6]], f2);
          }
          if (head_[i[7]]) {
            for_each_node_pair(head_[i[0]], head_[i[7]], f2);
          }
        }
        if (head_[i[1]]) {
          if (head_[i[2]]) {
            for_each_node_pair(head_[i[1]], head_[i[2]], f2);
          }
          if (head_[i[4]]) {
            for_each_node_pair(head_[i[1]], head_[i[4]], f2);
          }
          if (head_[i[6]]) {
            for_each_node_pair(head_[i[1]], head_[i[6]], f2);
          }
        }
        if (head_[i[2]]) {
          if (head_[i[4]]) {
            for_each_node_pair(head_[i[2]], head_[i[4]], f2);
          }
          if (head_[i[5]]) {
            for_each_node_pair(head_[i[2]], head_[i[5]], f2);
          }
        }
        if (head_[i[3]] && head_[i[4]]) {
          for_each_node_pair(head_[i[3]], head_[i[4]], f2);
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc, typename TriFunc>
void CellListNode_3<TNode>::for_each_pair_fast(BiFunc f1, TriFunc f2,
                                           const Vec3i& ic_beg, const Vec3i& ic_end) const {
  int z[2];
  int y[2];
  int x[2];
  int i[8];
  for (z[0] = ic_beg.z; z[0] < ic_end.z; z[0]++) {
    double lz = 0;
    z[1] = z[0] + 1;
    if (z[1] >= n_.z) {
      z[1] -= n_.z;
      lz = gl_l_.z;
    }
    const int z_nx_ny[2] = { z[0] * nxny_, z[1] * nxny_ };
    for (y[0] = ic_beg.y; y[0] < ic_end.y; y[0]++) {
      double ly = 0;
      y[1] = y[0] + 1;
      if (y[1] >= n_.y) {
        y[1] -= n_.y;
        ly = gl_l_.y;
      }
      const int y_nx[2] = { y[0] * n_.x, y[1] * n_.y };
      for (x[0] = ic_beg.x; x[0] < ic_end.x; x[0]++) {
        double lx = 0;
        x[1] = x[0] + 1;
        if (x[1] >= n_.x) {
          x[1] -= n_.x;
          lx = gl_l_.x;
        }
        i[0] = x[0] + y_nx[0] + z_nx_ny[0];
        i[1] = x[1] + y_nx[0] + z_nx_ny[0];
        i[2] = x[0] + y_nx[1] + z_nx_ny[0];
        i[3] = x[1] + y_nx[1] + z_nx_ny[0];
        i[4] = x[0] + y_nx[0] + z_nx_ny[1];
        i[5] = x[1] + y_nx[0] + z_nx_ny[1];
        i[6] = x[0] + y_nx[1] + z_nx_ny[1];
        i[7] = x[1] + y_nx[1] + z_nx_ny[1];

        if (head_[i[0]]) {
          for_each_node_pair(head_[i[0]], f1);
          if (head_[i[1]]) {
            for_each_node_pair(head_[i[0]], head_[i[1]], Vec3d(lx, 0., 0.), f2);
          }
          if (head_[i[2]]) {
            for_each_node_pair(head_[i[0]], head_[i[2]], Vec3d(0., ly, 0.), f2);
          }
          if (head_[i[3]]) {
            for_each_node_pair(head_[i[0]], head_[i[3]], Vec3d(lx, ly, 0.), f2);
          }
          if (head_[i[4]]) {
            for_each_node_pair(head_[i[0]], head_[i[4]], Vec3d(0., 0., lz), f2);
          }
          if (head_[i[5]]) {
            for_each_node_pair(head_[i[0]], head_[i[5]], Vec3d(lx, 0., lz), f2);
          }
          if (head_[i[6]]) {
            for_each_node_pair(head_[i[0]], head_[i[6]], Vec3d(0., ly, lz), f2);
          }
          if (head_[i[7]]) {
            for_each_node_pair(head_[i[0]], head_[i[7]], Vec3d(lx, ly, lz), f2);
          }
        }
        if (head_[i[1]]) {
          if (head_[i[2]]) {
            for_each_node_pair(head_[i[1]], head_[i[2]], Vec3d(-lx, ly, 0.), f2);
          }
          if (head_[i[4]]) {
            for_each_node_pair(head_[i[1]], head_[i[4]], Vec3d(-lx, 0., lz), f2);
          }
          if (head_[i[6]]) {
            for_each_node_pair(head_[i[1]], head_[i[6]], Vec3d(-lx, ly, lz), f2);
          }
        }
        if (head_[i[2]]) {
          if (head_[i[4]]) {
            for_each_node_pair(head_[i[2]], head_[i[4]], Vec3d(0., -ly, lz), f2);
          }
          if (head_[i[5]]) {
            for_each_node_pair(head_[i[2]], head_[i[5]], Vec3d(lx, -ly, lz), f2);
          }
        }
        if (head_[i[3]] && head_[i[4]]) {
          for_each_node_pair(head_[i[3]], head_[i[4]], Vec3d(-lx, -ly, lz), f2);
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc1, typename BiFunc2>
void CellListNode_3<TNode>::for_each_pair(BiFunc1 f1, BiFunc2 f2) const {
  Vec_3<int> beg(0, 0, 0);
  Vec_3<int> end(n_);
  if (flag_ext_.x)
    end.x -= 1;
  if (flag_ext_.y)
    end.y -= 1;
  if (flag_ext_.z)
    end.z -= 1;

  for_each_pair(f1, f2, beg, end);
}

template <typename TNode>
template <typename BiFunc, typename TriFunc>
void CellListNode_3<TNode>::for_each_pair_fast(BiFunc f1, TriFunc f2) const {
  Vec_3<int> end(n_);
  if (flag_ext_.x)
    end.x -= 1;
  if (flag_ext_.y)
    end.y -= 1;
  if (flag_ext_.z)
    end.z -= 1;

  for_each_pair_fast(f1, f2, Vec_3<int>(), end);
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
template <typename T>
void CellListNode_3<TNode>::create(std::vector<TNode>& p_arr, T* par_num_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head_[ic]);
    par_num_arr[ic] += 1;
  }
}

template <typename TNode>
void CellListNode_3<TNode>::recreate(std::vector<TNode>& p_arr) {
  for (int ic=0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
  }
  create(p_arr);
}

template <typename TNode>
template <typename T>
void CellListNode_3<TNode>::recreate(std::vector<TNode>& p_arr, T* par_num_arr) {
  for (int ic = 0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
    par_num_arr[ic] = 0;
  }
  create(p_arr);
}

template <typename TNode>
template <typename UniFunc>
void CellListNode_3<TNode>::for_each_cell(UniFunc f, const Vec_3<int>& beg, const Vec_3<int>& end) {
  for (int z = beg.z; z < end.z; z++) {
    const auto z_nx_ny = z * nxny_;
    for (int y = beg.y; y < end.y; y++) {
      const auto y_nx = y * n_.x;
      for (int x = beg.x; x < end.x; x++) {
        int ic = x + y_nx + z_nx_ny;
        f(&head_[ic]);
      }
    }
  }
}

template <typename TNode>
void CellListNode_3<TNode>::clear(const Vec_3<int>& first, const Vec_3<int>& last) {
  for (int z0 = first.z; z0 < last.z; z0++) {
    const auto z0_nxny = z0 * nxny_;
    for (int y0 = first.y; y0 < last.y; y0++) {
      const auto y0_nx = y0 * n_.x;
      for (int x0 = first.x; x0 < last.x; x0++) {
        int ic = x0 + y0_nx + z0_nxny;
        head_[ic] = nullptr;
      }
    }
  }
}

template <typename TNode>
template <typename TPar>
void CellListNode_3<TNode>::replace(BiNode<TPar>& p_new, const BiNode<TPar>& p_old) {
  p_new = p_old;
  if (p_new.next) {
    p_new.next->prev = &p_new;
  }
  if (p_new.prev) {
    p_new.prev->next = &p_new;
  } else {
    head_[get_ic(p_new)] = &p_new;
  }
}

template <typename TNode>
int CellListNode_3<TNode>::get_par_num(const Vec_3<int>& beg, const Vec_3<int>& end) const {
  int count = 0;
  for (int z0 = beg.z; z0 < end.z; z0++) {
    const auto z0_nxny = z0 * nxny_;
    for (int y0 = beg.y; y0 < end.y; y0++) {
      const auto y0_nx = y0 * n_.x;
      for (int x0 = beg.x; x0 < end.x; x0++) {
        int ic = x0 + y0_nx + z0_nxny;
        TNode * cur_node = head_[ic];
        while(cur_node) {
          count++;
          cur_node = cur_node->next;
        }
      }
    }
  }
  return count;
}

template <typename TNode>
int CellListNode_3<TNode>::get_par_num() const {
  return get_par_num(Vec_3<int>(), n_);
}
