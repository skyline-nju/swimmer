#ifndef CELLLIST2D_H
#define CELLLIST2D_H
#include <vector>
//#include <iterator>
#include "vect.h"
#include "node.h"
#include "cmdline.h"

class Cell_list_base_2 {
public:
  Cell_list_base_2(const Vec_2<double> &l_domain,
                   const Vec_2<double> &origin,
                   double r_cut,
                   const Vec_2<bool> &flag_comm);
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


template <typename TNode>
class Cell_list_2 {
public:
  typedef typename std::vector<TNode*>::iterator IT;
  typedef typename std::vector<TNode*>::const_iterator CIT;
  Cell_list_2(double lx, double ly,  // NOLINT
              double x0 = 0, double y0 = 0, double r_cut = 1,
              bool comm_x = false, bool comm_y = false) {
    set_para(lx, ly, x0, y0, r_cut, comm_x, comm_y); 
  }
  explicit Cell_list_2(const cmdline::parser &cmd);

  void set_para(double lx, double ly, double x0, double y0, double r_cut,
                bool comm_x, bool comm_y);
  int get_ic(const TNode &p) const;

  int get_row(int my_row, int drow) const;

  int get_col(int my_col, int dcol) const;

  template <typename BiFunc>
  void for_each_pair(BiFunc f_ij, int row_beg = 0, int row_end = -1, int col_beg = 0, int col_end = -1);

  template <typename BiFunc>
  void for_nearby_par(const TNode *p, BiFunc f_ij) const;

  void create(std::vector<TNode> &p_arr);

  void recreate(std::vector<TNode> &p_arr);


  template <typename UniFunc>
  void update(UniFunc move, int first_ic, CIT first, CIT last);

  template <typename UniFunc>
  void update(UniFunc move, int ic, TNode* head0);

  template <typename UniFunc>
  void update(UniFunc move, int first_ic, const std::vector<TNode*> &head0) {
    update(move, first_ic, head0.cbegin(), head0.cend());
  }

  template <typename UniFunc>
  void update(UniFunc move);

  template <typename UniFunc>
  void update_by_row(UniFunc move);

protected:
  bool comm_x_;
  bool comm_y_;
  int ncells_;
  double r_cut_;
  Vec_2<double> origin_;
  Vec_2<double> l_box_;
  Vec_2<double> inverse_lc_;
  Vec_2<double> l_cell_;
  Vec_2<int> n_bins_;
  std::vector<TNode*> head_;
};

template<typename TNode>  // NOLINT
Cell_list_2<TNode>::Cell_list_2(const cmdline::parser & cmd) {
  const auto lx = cmd.get<double>("Lx");
  const auto ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : lx;
  const auto x0 = 0;
  const auto y0 = 0;
  const auto r_cut = cmd.get<double>("sigma");
  const auto comm_x = false;
  const auto comm_y = false;
  set_para(lx, ly, x0, y0, r_cut, comm_x, comm_y);
}

template<typename TNode>
void Cell_list_2<TNode>::set_para(double lx, double ly,
                                  double x0, double y0,
                                  double r_cut,
                                  bool comm_x, bool comm_y) {
  r_cut_ = r_cut;
  comm_x_ = comm_x;
  comm_y_ = comm_y;
  if (!comm_x && !comm_y) {
    origin_.x = 0;
    origin_.y = 0;
    n_bins_.x = int(lx / r_cut);
    n_bins_.y = int(ly / r_cut);
    l_cell_.x = lx / n_bins_.x;
    l_cell_.y = ly / n_bins_.y;
    l_box_.x = lx;
    l_box_.y = ly;
  } else if (!comm_x && comm_y) {
    n_bins_.x = int(lx / r_cut);
    n_bins_.y = int(ly / r_cut) + 2;
    l_cell_.x = lx / n_bins_.x;
    l_cell_.y = ly / (n_bins_.y - 2);
    l_box_.x = lx;
    l_box_.y = ly + 2 * l_cell_.y;
    origin_.x = x0;
    origin_.y = y0 - l_cell_.y;
  } else if (comm_x && comm_y) {
    n_bins_.x = int(lx / r_cut) + 2;
    n_bins_.y = int(ly / r_cut) + 2;
    l_cell_.x = lx / (n_bins_.x - 2);
    l_cell_.y = ly / (n_bins_.y - 2);
    l_box_.x = lx + 2 * l_cell_.x;
    l_box_.y = ly + 2 * l_cell_.y;
    origin_.x = x0 - l_cell_.x;
    origin_.y = y0 - l_cell_.y;
  }
  inverse_lc_.x = 1 / l_cell_.x;
  inverse_lc_.y = 1 / l_cell_.y;
  ncells_ = n_bins_.x * n_bins_.y;
  head_.reserve(ncells_);
  for (int i = 0; i < ncells_; i++) {
    head_.emplace_back();
  }
}

template<typename TNode>
int Cell_list_2<TNode>::get_ic(const TNode & p) const {
  return int((p.x - origin_.x) * inverse_lc_.x)
    + int((p.y - origin_.y) * inverse_lc_.y) * n_bins_.x;
}

template <typename TNode>
int Cell_list_2<TNode>::get_row(int my_row, int drow) const {
  int row = my_row + drow;
  if(row < 0) {
    row += n_bins_.y;
  } else if (row >= n_bins_.y) {
    row -= n_bins_.y;
  }
  return row;
}

template <typename TNode>
int Cell_list_2<TNode>::get_col(int my_col, int dcol) const {
  int col = my_col + dcol;
  if(col < 0) {
    col += n_bins_.x;
  } else if (col >= n_bins_.x) {
    col -= n_bins_.x;
  }
  return col;  
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
void Cell_list_2<TNode>::for_each_pair(BiFunc f_ij, int row_beg, int row_end,
  int col_beg, int col_end) {
  if (row_end == -1)
    row_end = n_bins_.y;
  if (col_end == -1)
    col_end = n_bins_.x;
  for (int row = row_beg; row < row_end; row++) {
    int row_upper = row + 1;
    if (row_upper == n_bins_.y)
      row_upper = 0;
    const auto row_times_ncols = row * n_bins_.x;
    const auto row_upper_times_ncols = row_upper * n_bins_.x;
    for (int col = col_beg; col < col_end; col++) {
      int ic0 = col + row_times_ncols;
      int col_right = col + 1;
      if (col_right == n_bins_.x)
        col_right = 0;
      int ic1 = col_right + row_times_ncols;
      int ic2 = col + row_upper_times_ncols;
      if (head_[ic0]) {
        for_each_node_pair(head_[ic0], f_ij);
        int ic3 = col_right + row_upper_times_ncols;
        bool flag_c1_c2 = true;
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
  auto f_j = [p, f_ij](const TNode * p_j) {
    if (p != p_j) {
      f_ij(p, p_j);  
    }
  };
  const int my_ic = get_ic(*p);
  const int my_col = my_ic % n_bins_.x;
  const int my_row = my_ic / n_bins_.x;
  for (int drow = -1; drow <= 1; drow++) {
    auto row = get_row(my_row, drow);
    for (int dcol = -1; dcol <= 1; dcol++) {
      auto col = get_col(my_col, dcol);
      for_each_par(head_[col + row * n_bins_.x], f_j);
    }
  }
}

template<typename TNode>
template<typename UniFunc>
void Cell_list_2<TNode>::update_by_row(UniFunc move) {
  int nx = n_bins_.x;
  std::vector<TNode*> pre_row(head_.cend() - nx, head_.cend());
  auto first = head_.cbegin();
  auto last = first + nx;
  std::vector<TNode*> cur_row(first, last);
  first += nx;
  last += nx;
  std::vector<TNode*> nxt_row(first, last);
  update(move, 0, cur_row);
  cur_row.swap(nxt_row);

  for (int row = 1; row < n_bins_.y - 2; row++) {
    first += nx;
    last += nx;
    nxt_row.assign(first, last);
    update(move, row * nx, cur_row.cbegin(), cur_row.cend());
    cur_row.swap(nxt_row);
  }
  update(move, (n_bins_.y - 2) * nx, cur_row);
  update(move, (n_bins_.y - 1) * nx, pre_row);
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
#endif

