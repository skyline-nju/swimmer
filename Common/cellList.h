/**
 * @brief Linked-cell list
 * @file cellList.h
 * @author skyline-nju
 * @date 2018-04-23
 */
#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
#include "vect.h"
#include "node.h"

/*************************************************************************//**
 * \brief Base class for 2d linked-cell list.
 ****************************************************************************/
class CellListBase_2 {
public:
  /**
   * \brief Constructor 
   * \param l Domain size, 2d vector.
   * \param r_cut Cutoff radius
   * \param origin Origin point of the damain, 2d vector.
   */
  CellListBase_2(const Vec_2<double> &l, double r_cut,
                 const Vec_2<double> &origin);
  /**
   * \brief Visit the nearest cells and itself.
   * \param my_row The row index of the cell.
   * \param my_col The column index of the cell.
   * \param f Template function: void(int).
   */
  template<typename UniFunc>
  void visit_nearby_cell(int my_row, int my_col, UniFunc f) const;

  int get_row(double y) const {
    return int((y - origin_.y) * inverse_lc_.y);
  }

  int get_col(double x) const {
    return int((x - origin_.x) * inverse_lc_.x);
  }

  template <typename TPar>
  int get_ic(const TPar &p) const {
    return get_col(p.x) + n_.x * get_row(p.y);
  }

  int get_upper_row(int my_row) const;

  int get_lower_row(int my_row) const;

  int get_left_col(int my_col) const;

  int get_right_col(int my_col) const;

  int get_row(int my_row, int row_offset) const;

  int get_col(int my_col, int col_offset) const;
protected:
  int ncells_;                //!< number of cells
  Vec_2<int> n_;              //!< number of rows, columns
  Vec_2<double> origin_;      //!< origin point of the domain
  Vec_2<double> inverse_lc_;  //!< inverse length of one cell
};

inline CellListBase_2::CellListBase_2(const Vec_2<double>& l, double r_cut,
                                      const Vec_2<double> &origin): origin_(origin) {
  n_.x = int(l.x / r_cut);
  n_.y = int(l.y / r_cut);
  inverse_lc_.x = n_.x / l.x;
  inverse_lc_.y = n_.y / l.y;
  ncells_ = n_.x * n_.y;
}

inline int CellListBase_2::get_upper_row(int my_row) const {
  auto row = my_row + 1;
  if (row == n_.y)
    row = 0;
  return row;
}

inline int CellListBase_2::get_lower_row(int my_row) const {
  return my_row == 0 ? n_.y - 1 : my_row - 1;
}

inline int CellListBase_2::get_left_col(int my_col) const {
  return my_col == 0 ? n_.x - 1 : my_col - 1;
}

inline int CellListBase_2::get_right_col(int my_col) const {
  auto col = my_col + 1;
  if (col == n_.x) {
    col = 0;
  }
  return col;
}

inline int CellListBase_2::get_row(int my_row, int row_offset) const {
  auto row = my_row + row_offset;
  if (row < 0) {
    row += n_.y;
  } else if (row >= n_.y) {
    row -= n_.y;
  }
  return row;
}

inline int CellListBase_2::get_col(int my_col, int col_offset) const {
  auto col = my_col + col_offset;
  if (col < 0) {
    col += n_.x;
  } else if (col >= n_.x) {
    col -= n_.x;
  }
  return col;
}

template<typename UniFunc>
void CellListBase_2::visit_nearby_cell(int my_row, int my_col, UniFunc f) const {
  for (int drow = -1; drow <= 1; drow++) {
    const auto row = get_row(my_row, drow);
    for (int dcol = -1; dcol <= 1; dcol++) {
      const auto col = get_col(my_col, dcol);
      f(col + row * n_.x);
    }
  }
}

/*************************************************************************//**
 * \brief Linked-cell list recording the index of the particles.
 * 
 * Use the STL container such as std::list or std::forward_list to record the
 * index of particles. Each cell is one such container.
 * 
 * \tparam TContainer Contain the indexs of particles located in one cell. 
 ***************************************************************************/
template <typename TContainer>
class CellListIdx_2: public CellListBase_2 {
public:
  /**
  * \brief Constructor
  * \param l Domain size, 2d vector.
  * \param r_cut Cutoff radius
  * \param origin Origin point of the damain, 2d vector, default zero.
  */
  CellListIdx_2(const Vec_2<double> &l, double r_cut,
                const Vec_2<double> &origin = Vec_2<double>())
    : CellListBase_2(l, r_cut, origin), head_(ncells_) {}

  /**
   * \brief Create cell list
   * \tparam TPar Template particle.
   * \param p_arr Array of particles.
   */
  template <typename TPar>
  void create(const std::vector<TPar> &p_arr);

  /**
   * \brief Visit nearby particles sorrounding one particle.
   * \tparam TPar   Template particle.
   * \tparam BiFunc Template function: void(TPar*, TPar*)
   * \param p       Pointer of the target particle
   * \param p_arr   Array of particles.
   * \param f_ij A  Function acting on the target particle and its' neighbor.
   */
  template <typename TPar, typename BiFunc>
  void for_nearby_par(const TPar* p, const std::vector<TPar> &p_arr,
                      BiFunc f_ij) const;
  
protected:
  std::vector<TContainer> head_;
};

template<typename TContainer>
template<typename TPar>
void CellListIdx_2<TContainer>::create(const std::vector<TPar> & p_arr) {
  for (unsigned int i = 0; i < p_arr.size(); i++) {
    int ic = get_ic(p_arr[i]);
    head_[ic].push_front(i);
  }
}

template<typename TContainer>
template<typename TPar, typename BiFunc>
void CellListIdx_2<TContainer>::for_nearby_par(const TPar * p,
                                               const std::vector<TPar>& p_arr,
                                               BiFunc f_ij) const {
  auto lambda = [this, p, &p_arr, f_ij](int ic) {
    if (!head_[ic].empty()) {
      const auto end = head_[ic].cend();
      for (auto it = head_[ic].cbegin(); it != end; ++it) {
        const TPar* pj = &p_arr[*it];
        if (p != pj) {
          f_ij(p, pj);
        }
      }
    }
  };
  visit_nearby_cell(get_row(p->y), get_col(p->x), lambda);
}

/*************************************************************************//**
 * \brief Linked-cell list constituted by bidirectional nodes.
 * 
 * Each cell contains one pointer of bidirectional node. Then first node has the
 * pointer pointing the second one, and so on. 
 * 
 * \tparam TNode Template node, which should contain two pointers: prev, next.
 ****************************************************************************/
template <typename TNode>
class CellListNode_2:public CellListBase_2 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  /**
  * \brief Constructor
  * \param l Domain size, 2d vector.
  * \param r_cut Cutoff radius
  * \param origin Origin point of the damain, 2d vector, defalut zero.
  */
  CellListNode_2(const Vec_2<double> &l, double r_cut,
                 const Vec_2<double> &origin = Vec_2<double>())
    : CellListBase_2(l, r_cut, origin), head_(ncells_) {}

  /**
  * \brief Visit particle pairs located in cell[row_beg: row_end, col_beg: col_end].
  * \tparam BiFunc   Template function: void(TNode*, TNode*)
  * \param f_ij      Function acting on two particles.
  * \param row_beg   Beginning row index.
  * \param row_end   Ending row index.
  * \param col_beg   Beginnig column index.
  * \param col_end   Ending column index.
  */
  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij, int row_beg, int row_end,
                     int col_beg, int col_end);

  /**
   * \brief Visit all particle pairs.
   * \tparam BiFunc Template function: void(TNode*, TNode*)
   * \param f_ij    Function acting on two particles.
   */
  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij) {
    for_each_pair(f_ij, 0, n_.y, 0, n_.x);
  }

 /**
  * @brief Visit nearby particles surronding one particle
  * 
  * @tparam BiFunc Template function: void(TPar*, TPar*)
  * @param p       Pointer of the target particle
  * @param f_ij    Function acting on the target particle and its' neighbor.
  * @return int    Index of the cell where the particle is.
  */
  template <typename BiFunc>
  int for_nearby_par(TNode *p, BiFunc f_ij) const;

  template <typename BiFunc>
  bool for_nearby_par(TNode& p, BiFunc u_ij, int* cell_idx = nullptr) const;

  /**
  * \brief Create cell list
  * \tparam TPar Template particle.
  * \param p_arr Array of particles.
  */
  void create(std::vector<TNode> &p_arr);

  /**
   * \brief Recreate cell list
   * \param p_arr Array of particles.
   */
  void recreate(std::vector<TNode> &p_arr);

  /**
   * \brief Update particles in cell ic.
   * \tparam UniFunc Template function: void(TNode *)
   * \param move     Function to update the coordinates of a particle.
   * \param ic       Index of the target cell.
   * \param head0    The head pointer of cell ic before updating nearby cells.
   */
  template <typename UniFunc>
  void update(UniFunc move, int ic, TNode* head0);

  /**
   * \brief Update particles in cells from first to last
   * \tparam UniFunc Template function: void(TNode *)
   * \param move     Function to update the coordinates of a particle.
   * \param first_ic Index of the first cell.
   * \param first    Iterator of the first head pointer.
   * \param last     Iterator of the last head pointer.
   */
  template <typename UniFunc>
  void update(UniFunc move, int first_ic, CIT first, CIT last);

  /**
   * \brief Update particles in an array of cells.
   * \tparam UniFunc Template function: void(TNode *)
   * \param move     Function to update the coordinates of a particle.
   * \param first_ic Index of the first cell.
   * \param head0    Vector of head pointer.
   */
  template <typename UniFunc>
  void update(UniFunc move, int first_ic, const std::vector<TNode*> &head0) {
    update(move, first_ic, head0.cbegin(), head0.cend());
  }

  /**
   * \brief Update all particles.
   * \tparam UniFunc Template function: void(TNode *)
   * \param move     Function to update the coordinates of a particle.
   */
  template <typename UniFunc>
  void update(UniFunc move) {
    std::vector<TNode*> head_cur(head_);
    update(move, 0, head_cur.cbegin(), head_cur.cend());
  }

  /**
   * \brief Update all particles row by row.
   * 
   * Update the cells row by row, so that less memories are requried.
   * \tparam UniFunc Template function: void(TNode *)
   * \param move     Function to update the coordinates of a particle.
   */
  template <typename UniFunc>
  void update_by_row(UniFunc move);

  void update(TNode *p, int old_cell_idx, int new_cell_idx);

protected:
  std::vector<TNode*> head_;
};

template <typename TNode>
template <typename BiFunc>
void CellListNode_2<TNode>::for_each_pair(BiFunc f_ij, int row_beg, int row_end,
                                          int col_beg, int col_end) {
  for (int row = row_beg; row < row_end; row++) {
    const auto row_upper = get_upper_row(row);
    const auto row_times_ncols = row * n_.x;
    const auto row_upper_times_ncols = row_upper * n_.x;
    for (int col = col_beg; col < col_end; col++) {
      const auto ic0 = col + row_times_ncols;
      const auto col_right = get_right_col(col);
      const auto ic1 = col_right + row_times_ncols;
      const auto ic2 = col + row_upper_times_ncols;
      if (head_[ic0]) {
        for_each_node_pair(head_[ic0], f_ij);
        const auto ic3 = col_right + row_upper_times_ncols;
        auto flag_c1_c2 = true;
        if (head_[ic1])
          for_each_node_pair(head_[ic0], head_[ic1], f_ij);
        else
          flag_c1_c2 = false;
        if (head_[ic2])
          for_each_node_pair(head_[ic0], head_[ic2], f_ij);
        else
          flag_c1_c2 = false;
        if (flag_c1_c2)
          for_each_node_pair(head_[ic1], head_[ic2], f_ij);
        if (head_[ic3])
          for_each_node_pair(head_[ic0], head_[ic3], f_ij);
      } else if (head_[ic1] && head_[ic2]) {
        for_each_node_pair(head_[ic1], head_[ic2], f_ij);
      }
    }
  }
}

template <typename TNode>
template <typename BiFunc>
int CellListNode_2<TNode>::for_nearby_par(TNode* p, BiFunc f_ij) const {
  int my_row = get_row(p->y);
  int my_col = get_col(p->x);
  int cell_idx = my_col + my_row * n_.x;
  visit_nearby_cell(my_row, my_col, [this, p, f_ij](int ic) {
    if (head_[ic]) {
      TNode *cur_node = head_[ic];
      do {
        if (p != cur_node)
          f_ij(p, cur_node);
        cur_node = cur_node->next;
      } while (cur_node);
    }
  });
  return cell_idx;
}

template <typename TNode>
template <typename BiFunc>
bool CellListNode_2<TNode>::for_nearby_par(TNode& p, BiFunc u_ij, int* ptr_ic0) const {
  int my_row = get_row(p.y);
  int my_col = get_col(p.x);
  if (ptr_ic0) {
    *ptr_ic0 = my_col + my_row * n_.x;
  }

  bool success = true;
  for (int drow = -1; drow <= 1; drow++) {
    const auto row = get_row(my_row, drow);
    for (int dcol = -1; dcol <= 1; dcol++) {
      const auto col = get_col(my_col, dcol);
      int ic = col + row * n_.x;
      if (head_[ic]) {
        TNode* cur_node = head_[ic];
        do {
          if (&p != cur_node) {
            if (!u_ij(p, *cur_node)) {
              success = false;
              break;
            }
          }
          cur_node = cur_node->next;
        } while (cur_node);
        
      }
      if (!success) {
        break;
      }
    }
    if (!success) {
      break;
    }
  }
  return success;
}

template<typename TNode>
void CellListNode_2<TNode>::create(std::vector<TNode>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head_[ic]);
  }
}

template<typename TNode>
void CellListNode_2<TNode>::recreate(std::vector<TNode>& p_arr) {
  for (int ic = 0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
  }
  create(p_arr);
}

template <typename TNode>
template <typename UniFunc>
void CellListNode_2<TNode>::update(UniFunc move, int ic, TNode* head0) {
  auto cur_node = head0;
  const auto flag = cur_node == head_[ic] ? true : false;
  do {
    move(cur_node);
    int ic_new = get_ic(*cur_node);
    if (ic_new != ic) {
      if (flag) {
        cur_node->break_away(&head_[ic]);
      } else {
        cur_node->break_away();
      }
      auto tmp_node = cur_node;
      cur_node = cur_node->next;
      tmp_node->append_at_front(&head_[ic_new]);
    } else {
      cur_node = cur_node->next;
    }
  } while (cur_node);
}

template <typename TNode>
template <typename UniFunc>
void CellListNode_2<TNode>::update(UniFunc move, int first_ic,
                                   CIT first, CIT last) {
  int ic = first_ic;
  for (auto it = first; it != last; ++it) {
    if (*it) {
      update(move, ic, *it);
    }
    ic++;
  }
}

template <typename TNode>
template <typename UniFunc>
void CellListNode_2<TNode>::update_by_row(UniFunc move) {
  std::vector<TNode*> pre_row(head_.cend() - n_.x, head_.cend());
  auto first = head_.cbegin();
  auto last = first + n_.x;
  std::vector<TNode*> cur_row(first, last);
  first += n_.x;
  last += n_.x;
  std::vector<TNode*> nxt_row(first, last);
  update(move, 0, cur_row);
  cur_row.swap(nxt_row);

  for (int row = 1; row < n_.y - 2; row++) {
    first += n_.x;
    last += n_.x;
    nxt_row.assign(first, last);
    update(move, row * n_.x, cur_row.cbegin(), cur_row.cend());
    cur_row.swap(nxt_row);
  }
  update(move, (n_.y - 2) * n_.x, cur_row);
  update(move, (n_.y - 1) * n_.x, pre_row);
}

template <typename TNode>
void CellListNode_2<TNode>::update(TNode *p, int old_cell_idx, int new_cell_idx) {
  if (new_cell_idx != old_cell_idx) {
    p->break_away(&head_[old_cell_idx]);
    p->append_at_front(&head_[new_cell_idx]);
  }

}
#endif
