#ifndef CELLLIST2D_H
#define CELLLIST2D_H
#include <vector>
#include "vect.h"
#include "node.h"

/**
 * Base class for linked-cell list for 2D 
 */
class Cell_list_base_2 {
public:
  Cell_list_base_2(const Vec_2<double> &l_domain,
                   const Vec_2<double> &origin,
                   double r_cut,
                   const Vec_2<bool> &flag_comm = Vec_2<bool>());
  template<typename TBc>
  Cell_list_base_2(const TBc& bc, double r_cut,
                   const Vec_2<bool>& flag_comm = Vec_2<bool>());

  int get_row(double y) const {return int((y - origin_.y) * inverse_l_cell_.y);}

  int get_col(double x) const {return int((x - origin_.x) * inverse_l_cell_.x);}

  template<typename TPar>
  int get_ic(const TPar &p) const {return get_col(p.x) + get_row(p.y) * n_.x;}

  int get_upper_row(int my_row) const;

  int get_lower_row(int my_row) const;

  int get_left_col(int my_col) const;

  int get_right_col(int my_col) const;

  int get_row(int my_row, int drow) const;

  int get_col(int my_col, int dcol) const;

  template<typename UniFunc>
  void visit_nearby_cell(int my_row, int my_col, UniFunc f) const;

protected:
  Vec_2<bool> flag_comm_;
  int ncells_;
  double r_cut_;
  Vec_2<int> n_;
  Vec_2<double> origin_;
  Vec_2<double> l_domain_;
  Vec_2<double> l_cell_;
  Vec_2<double> inverse_l_cell_;
};

/**
 * \brief  Class for cell list for 2D. 
 * \tparam TNode Template for node, such as BiNode<TPar>
 */
template <typename TNode>
class Cell_list_2:public Cell_list_base_2 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  template <typename TBc>
  Cell_list_2(const TBc &bc, double r_cut,
              const Vec_2<bool> &flag_comm=Vec_2<bool>());

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij, int row_beg = 0, int row_end = -1, int col_beg = 0, int col_end = -1);

  template <typename BiFunc>
  void for_nearby_par(const TNode *p, BiFunc f_ij) const;

  void create(std::vector<TNode> &p_arr);

  void recreate(std::vector<TNode> &p_arr);

  /**
   * \brief Update the location of particles in the ic-th cell.
   * \tparam UniFunc Template function with one input para
   * \param move Function dicatating how a particle moves. 
   * \param ic The index of the current cell.
   * \param head0 The head of ic-th cell before updating it and surrounding cells.
   */
  template <typename UniFunc>
  void update(UniFunc move, int ic, TNode* head0);

  /**
   * \brief Update the cell from the first to the last iterator.
   * \tparam UniFunc Template function with one input para
   * \param move Function dicatating how a particle moves.
   * \param first_ic The index of the starting cell.
   * \param first Iterator pointing to the original head of the first cell. 
   * \param last  Iterator pointing to the original head of the last cell.
   */
  template <typename UniFunc>
  void update(UniFunc move, int first_ic, CIT first, CIT last);

  template <typename UniFunc>
  void update(UniFunc move, int first_ic, const std::vector<TNode*> &head0) {
    update(move, first_ic, head0.cbegin(), head0.cend());
  }

  /**
   * \brief Update the whole cells.
   * \tparam UniFunc
   * \param move 
   */
  template <typename UniFunc>
  void update(UniFunc move);

  /**
   * \brief Update the whole cells by row by row.
   * \tparam UniFunc 
   * \param move 
   */
  template <typename UniFunc>
  void update_by_row(UniFunc move);

protected:
  std::vector<TNode*> head_;
};

/**
 * \brief  Class for cell list for 2D. 
 * \tparam TContainer STL container such as std::list<int>.
 */
template <typename TContainer> 
class Cell_list_idx_2: public Cell_list_base_2 {
public:
  template <typename TBc>
  Cell_list_idx_2(const TBc &bc, double r_cut,
                  const Vec_2<bool> &flag_comm = Vec_2<bool>());
  template <typename TPar>
  void create(const TPar &p_arr);

  template <typename TPar, typename BiFunc>
  void for_nearby_par(const TPar* p, const std::vector<TPar> &p_arr,
                      BiFunc f_ij) const;

protected:
  std::vector<TContainer> head_;
};

template <typename TBc>
Cell_list_base_2::Cell_list_base_2(const TBc& bc, double r_cut,
                                   const Vec_2<bool>& flag_comm):
    Cell_list_base_2(Vec_2<double>(bc.get_lx(), bc.get_ly()),
                     Vec_2<double>(bc.get_xmin(), bc.get_ymin()),
                     r_cut, flag_comm) {}

template <typename UniFunc>
void Cell_list_base_2::visit_nearby_cell(int my_row, int my_col, UniFunc f) const {
  for (int drow = -1; drow <= 1; drow++) {
    const auto row = get_row(my_row, drow);
    for (int dcol = -1; dcol <= 1; dcol++) {
      const auto col = get_col(my_col, dcol);
      f(col + row * n_.x);
    }
  }
}

inline int Cell_list_base_2::get_upper_row(int my_row) const {
  auto row = my_row + 1;
  if (row == n_.y)
    row = 0;
  return row;
}

inline int Cell_list_base_2::get_lower_row(int my_row) const {
  return my_row == 0? n_.y - 1: my_row - 1;
}

inline int Cell_list_base_2::get_left_col(int my_col) const {
  return my_col == 0? n_.x - 1: my_col - 1;
}

inline int Cell_list_base_2::get_right_col(int my_col) const {
  auto col = my_col + 1;
  if (col == n_.x) {
    col = 0;
  }
  return col;
}

template <typename TNode>
template <typename TBc>
Cell_list_2<TNode>::Cell_list_2(const TBc& bc, double r_cut,
                                const Vec_2<bool> &flag_comm):
                                Cell_list_base_2(bc, r_cut, flag_comm) {
  head_.reserve(ncells_);
  for (int i = 0; i < ncells_; i++) {
    head_.emplace_back();
  }
}

template<typename TNode>
void Cell_list_2<TNode>::create(std::vector<TNode>& p_arr) {
  auto end = p_arr.end();
  for (auto it = p_arr.begin(); it != end; ++it) {
    int ic = get_ic(*it);
    (*it).append_at_front(&head_[ic]);
  }
}

template<typename TNode>
void Cell_list_2<TNode>::recreate(std::vector<TNode>& p_arr) {
  for (int ic = 0; ic < ncells_; ic++) {
    head_[ic] = nullptr;
  }
  create(p_arr);
}

template<typename TNode>
template<typename BiFunc>
void Cell_list_2<TNode>::for_each_pair(BiFunc f_ij,
                                       int row_beg, int row_end,
                                       int col_beg, int col_end) {
  if (row_end == -1)
    row_end = n_.y;
  if (col_end == -1)
    col_end = n_.x;
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

template<typename TNode>
template<typename BiFunc>
void Cell_list_2<TNode>::for_nearby_par(const TNode * p, BiFunc f_ij) const {
  visit_nearby_cell(get_row(p->y), get_col(p->x), [this, p, f_ij](int ic) {
    if(head_[ic]) {
      const TNode *cur_node = head_[ic];
      do {
        if (p != cur_node)
          f_ij(p, cur_node);
        cur_node = cur_node->next;
      } while(cur_node);
    }
  });
}


template <typename TContainer>
template <typename TBc>
Cell_list_idx_2<TContainer>::Cell_list_idx_2(const TBc& bc, double r_cut,
                                             const Vec_2<bool>& flag_comm):
                                             Cell_list_base_2(bc, r_cut, flag_comm) {
  head_.reserve(ncells_);
  for (int i = 0; i < ncells_; i++) {
    head_.emplace_back();
  }
}

template <typename TContainer>
template <typename TPar>
void Cell_list_idx_2<TContainer>::create(const TPar& p_arr) {
  for (unsigned int i = 0; i < p_arr.size(); i++) {
    int ic = get_ic(p_arr[i]);
    head_[ic].push_front(i);
  }
}

template <typename TContainer>
template <typename TPar, typename BiFunc>
void Cell_list_idx_2<TContainer>::for_nearby_par(const TPar *p,
                                                 const std::vector<TPar> &p_arr,
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

template <typename TNode>
template <typename UniFunc>
void Cell_list_2<TNode>::update(UniFunc move, int ic, TNode* head0) {
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

template<typename TNode>
template<typename UniFunc>
void Cell_list_2<TNode>::update(UniFunc move, int first_ic, CIT first, CIT last) {
  int ic = first_ic;
  for (auto it = first; it != last; ++it) {
    if (*it) {
      update(move, ic, *it);
    }
    ic++;
  }
}

template<typename TNode>
template<typename UniFunc>
void Cell_list_2<TNode>::update(UniFunc move) {
  std::vector<TNode*> head_cur(head_);
  update(move, 0, head_cur.cbegin(), head_cur.cend());
}

template<typename TNode>
template<typename UniFunc>
void Cell_list_2<TNode>::update_by_row(UniFunc move) {
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

#endif

