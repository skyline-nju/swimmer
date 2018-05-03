#include "latticeWetting.h"
void lattice::find_neighbors(int ic, std::list<int>& neighbor,
                             const Vec_2<int> &l_domain,
                             const std::vector<uint8_t>& cell,
                             const Vec_2<bool> &flag_bc) {
  auto lambda = [&neighbor, &l_domain, &cell](int x, int y) {
    const auto i = x + y * l_domain.x;
    if (cell[i])
      neighbor.push_back(i);
  };

  int y = ic / l_domain.x;
  int x = ic % l_domain.x;

  int x_l = x - 1;
  if (x_l < 0) {
    if (flag_bc.x) {
      x_l += l_domain.x;
      lambda(x_l, y);
    }
  } else {
    lambda(x_l, y);
  }
  int x_r = x + 1;
  if (x_r >= l_domain.x) {
    if (flag_bc.x) {
      x_r = 0;
      lambda(x_r, y);
    }
  } else {
    lambda(x_r, y);
  }

  int y_d = y - 1;
  if (y_d < 0) {
    if (flag_bc.y) {
      y_d += l_domain.y;
      lambda(x, y_d);
    }
  } else {
    lambda(x, y_d);
  }
  int y_u = y + 1;
  if (y_u >= l_domain.y) {
    if (flag_bc.y) {
      y_u = 0;
      lambda(x, y_u);
    }
  } else {
    lambda(x, y_u);
  }
}

Vec_2<int> lattice::Cluster::domain_len;

lattice::Cluster::Cluster(int ic, int size) {
  idx_cell_.reserve(size);
  idx_cell_.push_back(ic);
  x_min_ = x_max_ = ic % domain_len.x;
  y_min_ = y_max_ = ic / domain_len.x;
}

void lattice::Cluster::expand(const std::vector<uint8_t>& cell,
                              std::list<int>& neighbor,
                              std::vector<bool>& flag_visited,
                              std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end(); ++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      find_neighbors(k, new_neighbor_idx, domain_len, cell,
                     Vec_2<bool>(false, true));
      if (!new_neighbor_idx.empty()) {
        neighbor.splice(neighbor.end(), new_neighbor_idx);
      }
    }
    if (!flag_clustered[k]) {
      flag_clustered[k] = true;
      push_back(k);
    }
  }
  idx_cell_.reserve(idx_cell_.size());
}

void lattice::Cluster::check_wetting(std::vector<char>& flag_w) const {
  if (x_min_ == 0) {
    for (auto i: idx_cell_) {
      flag_w[i] = 'L';
    }
  } else if (x_max_ == domain_len.x - 1) {
    for (auto i : idx_cell_) {
      flag_w[i] = 'R';
    }
  }
}

void lattice::Cluster::find_clusters(std::vector<Cluster>& c_arr,
                                     const std::vector<uint8_t>& cell) {
  const auto n_cells = cell.size();
  c_arr.reserve(n_cells / 2);
  std::vector<bool> flag_clustered(n_cells, false);
  std::vector<bool> flag_visited(n_cells, false);
  for (size_t i = 0; i < n_cells; i++) {
    if (cell[i] && (!flag_visited[i])) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      find_neighbors(i, neighbor, domain_len, cell,
                     Vec_2<bool>(false, true));
      if (!neighbor.empty()) {
        c_arr.emplace_back(i, n_cells);
        flag_clustered[i] = true;
        c_arr.back().expand(cell, neighbor, flag_visited,
                            flag_clustered);
      }
    }
  }
}

void lattice::Cluster::find_clusters(std::vector<Cluster>& c_arr,
                                     const std::vector<uint8_t>& cell,
                                     std::vector<char>& flag_wetting) {
  const auto n_cells = cell.size();
  c_arr.reserve(n_cells / 2);
  std::vector<bool> flag_clustered(n_cells, false);
  std::vector<bool> flag_visited(n_cells, false);
  for (size_t i = 0; i < n_cells; i++) {
    if (cell[i] && (!flag_visited[i])) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      find_neighbors(i, neighbor, domain_len, cell, Vec_2<bool>(false, true));
      if (!neighbor.empty()) {
        c_arr.emplace_back(i, n_cells);
        flag_clustered[i] = true;
        c_arr.back().expand(cell, neighbor, flag_visited,
                            flag_clustered);
        c_arr.back().check_wetting(flag_wetting);
      }
    }
  }
  for (size_t i = 0; i < n_cells; i++) {
    if (cell[i]) {
      int x = i % domain_len.x;
      if (x == 0)
        flag_wetting[i] = 'L';
      else if (x == domain_len.x - 1)
        flag_wetting[i] = 'R';
    }
  }
}

void lattice::Cluster::push_back(int ic) {
  idx_cell_.push_back(ic);
  int x = ic % domain_len.x;
  int y = ic / domain_len.y;
  if (x_min_ > x)
    x_min_ = x;
  else if (x_max_ < x)
    x_max_ = x;
  if (y_min_ > y)
    y_min_ = y;
  else if (y_max_ < y)
    y_max_ = y;
}

void lattice::cal_non_empty_frac(const std::vector<uint8_t>& cell,
                                std::vector<double>& packing_frac) {
  std::vector<uint16_t> thickness_profile;
  std::vector<uint16_t> num_profile;
  cal_non_empty_profile(cell, thickness_profile, num_profile, packing_frac);
}

