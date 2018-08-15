#pragma once
#include <vector>
#include "vect.h"
#include "node.h"

void divide_cell(const Vec_3<double> &l, double r_cut,
                 Vec_3<int> &cell_num, Vec_3<double> &cell_len);

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

  int ncells() const { return ncells_; }

  void show_para() const;

  const Vec3d& origin() const { return origin_; }
  const Vec3d& l() const { return l_; }
  const Vec_3<bool>& flag_ext() const { return flag_ext_; }

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

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij, const Vec_3<int> &ic_beg,
                     const Vec_3<int> &ic_end) const;

  template <typename BiFunc>
  void for_each_pair2(BiFunc f_ij, const Vec_3<int> &ic_beg,
                      const Vec_3<int> &ic_end) const;

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij) const;

  void create(std::vector<TNode> &p_arr);

  template <typename T>
  void create(std::vector<TNode> &p_arr, T* par_num_arr);

  void recreate(std::vector<TNode> &p_arr);

  template <typename T>
  void recreate(std::vector<TNode> &p_arr, T* par_num_arr);

  template <typename UniFunc>
  void for_each_node(UniFunc f, const Vec_3<int> &first,
                     const Vec_3<int> &last) const;

  template <typename UniFunc>
  void for_each_node_clean(UniFunc f, const Vec_3<int> &first,
                           const Vec_3<int> &last);

  void add_node(TNode &p) { p.append_at_front(&head_[get_ic[p]]); }

  void clear(const Vec_3<int> &first, const Vec_3<int> &last);
  
  template <typename TPar>
  void replace(BiNode<TPar> &p_new, const BiNode<TPar> &p_old);
protected:
  std::vector<TNode*> head_;
};

template <typename TNode>
template <typename BiFunc>
void CellListNode_3<TNode>::for_each_pair(BiFunc f_ij,
                                          const Vec_3<int>& ic_beg,
                                          const Vec_3<int>& ic_end) const {
  auto cell_cell = [f_ij, this](int i1, int i2) {
    if (this->head_[i2]) {
      for_each_node_pair(this->head_[i1], this->head_[i2], f_ij);
    }
  };
  for (int z0 = ic_beg.z; z0 < ic_end.z; z0++) {
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
        int i0 = x0 + y0_nx + z0_nxny;
        int i1 = x1 + y0_nx + z0_nxny;
        int i2 = x0 + y1_nx + z0_nxny;
        int i3 = x1 + y1_nx + z0_nxny;
        int i4 = x0 + y0_nx + z1_nxny;
        int i5 = x1 + y0_nx + z1_nxny;
        int i6 = x0 + y1_nx + z1_nxny;
        int i7 = x1 + y1_nx + z1_nxny;
        if (head_[i0]) {
          for_each_node_pair(head_[i0], f_ij);
          cell_cell(i0, i1);
          cell_cell(i0, i2);
          cell_cell(i0, i3);
          cell_cell(i0, i4);
          cell_cell(i0, i5);
          cell_cell(i0, i6);
          cell_cell(i0, i7);
        }
        if (head_[i1]) {
          cell_cell(i1, i2);
          cell_cell(i1, i4);
          cell_cell(i1, i6);
        }
        if (head_[i2]) {
          cell_cell(i2, i4);
          cell_cell(i2, i5);
        }
        if (head_[i3]) {
          cell_cell(i3, i4);
        }
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc>
void CellListNode_3<TNode>::for_each_pair2(BiFunc f_ij,
                                           const Vec_3<int>& ic_beg,
                                           const Vec_3<int>& ic_end) const {
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
template <typename BiFunc>
void CellListNode_3<TNode>::for_each_pair(BiFunc f_ij) const {
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
void CellListNode_3<TNode>::for_each_node(UniFunc f, const Vec_3<int>& first,
                                          const Vec_3<int>& last) const {
  for (int z0 = first.z; z0 < last.z; z0++) {
    const auto z0_nxny = z0 * nxny_;
    for (int y0 = first.y; y0 < last.y; y0++) {
      const auto y0_nx = y0 * n_.x;
      for (int x0 = first.x; x0 < last.x; x0++) {
        int ic = x0 + y0_nx + z0_nxny;
        for_each_node(head_[ic], f);
      }
    }
  }
}

template <typename TNode>
template <typename UniFunc>
void CellListNode_3<TNode>::for_each_node_clean(UniFunc f,
                                                const Vec_3<int>& first,
                                                const Vec_3<int>& last) {
  for (int z0 = first.z; z0 < last.z; z0++) {
    const auto z0_nxny = z0 * nxny_;
    for (int y0 = first.y; y0 < last.y; y0++) {
      const auto y0_nx = y0 * n_.x;
      for (int x0 = first.x; x0 < last.x; x0++) {
        int ic = x0 + y0_nx + z0_nxny;
        for_each_node(head_[ic], f);
        head_[ic] = nullptr;
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
