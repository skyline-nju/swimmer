#ifndef WETTING_H
#define WETTING_H
#include "boundary.h"
#include "cellList2D.h"
#include "exporter.h"
#include <list>
#include <forward_list>
#include <cassert>

template <typename TNode, typename TCellList, typename TBc>
void find_neighbors(int i, double eps_square, const TNode *p0,
                    const TCellList &cl, const TBc &bc, std::list<int> &neighbor) {
  auto lambda = [p0, eps_square, &bc, &neighbor](const TNode* pi, const TNode* pj) {
    if (get_dis_square(*pi, *pj, bc) < eps_square) {
      neighbor.push_back(pj - p0);
    }
  };
  cl.for_nearby_par(p0 + i, lambda);
}

template <typename TPar, typename TCellList, typename TBc>
void find_neighbors(int i, double eps_square, const std::vector<TPar> &p_arr,
                    const TCellList &cl, const TBc &bc, std::list<int> &neighbor) {
  const auto p0 = &p_arr[0];
  auto lambda = [p0, eps_square, &bc, &neighbor](const TPar *pi, const TPar *pj) {
    if(get_dis_square(*pi, *pj, bc) < eps_square) {
      neighbor.push_back(pj - p0);
    }
  };
  cl.for_nearby_par(&p_arr[i], p_arr, lambda);
}

class Cluster {
public:
  Cluster()=default;

  explicit Cluster(int idx, int size) {
    idx_arr_.reserve(size);
    idx_arr_.push_back(idx);
  }

  explicit Cluster(int idx, int size, double x) {
    idx_arr_.reserve(size);
    idx_arr_.push_back(idx);
  }

  template <typename TNode, typename TCellList, typename TBc>
  void expand(double eps_square, unsigned int min_pts,
              const std::vector<TNode>& p_arr,
              const TCellList &cl,
              const TBc &bc,
              std::list<int> &neighbor,
              std::vector<bool> &flag_visited,
              std::vector<bool> &flag_clustered);

  template <typename UniFunc>
  void for_each(UniFunc f_i) const;

  int get_n() const {return idx_arr_.size();}

  template <typename T>
  void mark(std::vector<T> &flag_arr, T flag) const {
    for_each([flag, &flag_arr](auto i){flag_arr[i] = flag;});
  }

protected:
  std::vector<int> idx_arr_;
};

class Cluster_w_xlim : public Cluster {
public:
  Cluster_w_xlim() {x_min_=x_max_=0;}

  explicit Cluster_w_xlim(int idx, int size, double x):
      Cluster(idx, size), x_min_(x), x_max_(x) { }

  template <typename TNode, typename TCellList, typename TBc>
  void expand(double eps_square, int min_pts,
              const std::vector<TNode> &p_arr,
              const TCellList &cl,
              const TBc &bc, std::list<int> &neighbor,
              std::vector<bool> &flag_visited,
              std::vector<bool> &flag_clustered);

  void mark_wetting(double x_left, double x_right,
                    std::vector<char>& flag_w) const;
protected:
  double x_min_;
  double x_max_;
};

inline void Cluster_w_xlim::mark_wetting(double x_left, double x_right,
                                         std::vector<char>& flag_w) const {
  if (x_min_ < x_left)
    mark(flag_w, 'L');
  else if (x_max_ >= x_right)
    mark(flag_w, 'R');
}

template <typename TNode, typename TCellList, typename TBc>
void Cluster::expand(double eps_square, unsigned int min_pts,
                     const std::vector<TNode>& p_arr,
                     const TCellList& cl,
                     const TBc& bc,
                     std::list<int> &neighbor,
                     std::vector<bool>& flag_visited,
                     std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end();++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      find_neighbors(k, eps_square, p_arr, cl, bc, neighbor);
      if (new_neighbor_idx.size() >= min_pts) {
        neighbor.splice(neighbor.end(), new_neighbor_idx);
      } 
    }
    if (!flag_clustered[k]){
      flag_clustered[k] = true;
      idx_arr_.push_back(k);
    }
  }
  idx_arr_.shrink_to_fit();
}

template <typename UniFunc>
void Cluster::for_each(UniFunc f_i) const {
  const auto end = idx_arr_.cend();
  for (auto it = idx_arr_.cbegin(); it != end; ++it)
    f_i(*it);
}

template <typename TNode, typename TCellList, typename TBc>
void Cluster_w_xlim::expand(double eps_square, int min_pts,
                            const std::vector<TNode>& p_arr,
                            const TCellList& cl,
                            const TBc& bc, std::list<int>& neighbor,
                            std::vector<bool>& flag_visited,
                            std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end();++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      find_neighbors(k, eps_square, p_arr, cl, bc, neighbor);
      if (new_neighbor_idx.size() >= min_pts) {
        neighbor.splice(neighbor.end(), new_neighbor_idx);
      } 
    }
    if (!flag_clustered[k]){
      flag_clustered[k] = true;
      idx_arr_.push_back(k);
      // only valid for non-periodic boundary condition in the x direction
      if (p_arr[k].x < x_min_) {
        x_min_ = p_arr[k].x;
      }
      if (p_arr[k].x > x_max_) {
        x_max_ = p_arr[k].x;
      }
    }
  }
  idx_arr_.shrink_to_fit();
}

template <typename TPar, typename TBc, typename TCluster>
void dbscan(std::vector<TCluster> &c_arr,
            std::vector<bool> &flag_clustered,
            double eps, unsigned int min_pts,
            const std::vector<TPar> &p_arr,
            const TBc &bc) {
  const auto eps_square = eps * eps;
  const auto n_par = p_arr.size();
  c_arr.reserve(n_par);
  std::vector<bool> flag_visited(n_par, false);
  Cell_list_idx_2<std::forward_list<int>> cl(bc, eps);
  cl.create(p_arr);
  for (unsigned int i = 0; i < n_par; i++) {
    if (!flag_visited[i]) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      find_neighbors(i, eps_square, p_arr, cl, bc, neighbor);
      if (neighbor.size() >= min_pts) {
        c_arr.emplace_back(i, n_par, p_arr[i].x);
        flag_clustered[i] = true;
        c_arr.back().expand(eps_square, min_pts, p_arr, cl, bc,
                            neighbor, flag_visited, flag_clustered);
      }
    }
  }
}

template <typename TPar, typename TBc>
void dbscan_wall(std::vector<Cluster_w_xlim> &c_arr, 
                 std::vector<bool> &flag_clustered,
                 std::vector<char> &flag_wetting,
                 double eps, unsigned int min_pts, double height_min,
                 const std::vector<TPar> &p_arr,
                 const TBc &bc) {
  const auto eps_square = eps * eps;
  double x_left = height_min;
  double x_right = bc.get_xmin() + bc.get_lx() - height_min;
  const auto n_par = p_arr.size();
  c_arr.reserve(n_par);
  std::vector<bool> flag_visited(n_par, false);
  Cell_list_idx_2<std::forward_list<int>> cl(bc, eps);
  cl.create(p_arr);
  for (unsigned int i = 0; i < n_par; i++) {
    if (!flag_visited[i]) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      find_neighbors(i, eps_square, p_arr, cl, bc, neighbor);
      if (neighbor.size() >= min_pts) {
        c_arr.emplace_back(i, n_par, p_arr[i].x);
        flag_clustered[i] = true;
        c_arr.back().expand(eps_square, min_pts, p_arr, cl, bc,
                            neighbor, flag_visited, flag_clustered);
        c_arr.back().mark_wetting(x_left, x_right, flag_wetting);
      } else if (p_arr[i].x <= x_left) {
        flag_wetting[i] = 'L';
      } else if (p_arr[i].x >= x_right) {
        flag_wetting[i] = 'R';
      }
    }
  }
}

template <typename TCluster>
void get_cluster_size(const std::vector<TCluster> &cluster,
                      std::vector<int> &size_arr) {
  const auto end = cluster.cend();
  for (auto it = cluster.cbegin(); it != end; ++it) {
    (*it).mark(size_arr, (*it).get_n());
  }
}

template <typename T>
unsigned int get_n_flag(const std::vector<T> &flag_arr, T flag) {
  unsigned int n = 0;
  const auto end = flag_arr.cend();
  for (auto it = flag_arr.cbegin(); it != end; ++it) {
    if (*it == flag)
      n++;
  }
  return n;
}

/* Cal the thickness profile and corresponding particle number */

template <typename TPar, typename TBc>
void cal_wetting_profile(std::vector<float> &thickness_profile,
                         std::vector<unsigned short> &num_profile,
                         std::vector<double> &packing_frac,
                         const std::vector<TPar> &p_arr,
                         const std::vector<char> &flag_wetting,
                         const TBc &bc) {
  const auto nrow = int(bc.get_ly());
  const auto x_min = bc.get_xmin();
  const auto x_max = x_min + bc.get_lx();
  const auto y_min = bc.get_ymin();
  const auto n = p_arr.size();
  for (unsigned int i = 0; i < n; i++) {
    if (flag_wetting[i] == 'L') {
      const auto row = int(p_arr[i].y - y_min);
     double dx = p_arr[i].x - x_min;
      if (thickness_profile[row] < dx) {
        thickness_profile[row] = dx;
      }
      if (num_profile[row] == 0) {
        packing_frac[0] += 1;
      }
      num_profile[row]++;
    } else if (flag_wetting[i] == 'R') {
      const auto row = int(p_arr[i].y - y_min) + nrow;
      double dx = x_max - p_arr[i].x;
      if (thickness_profile[row] < dx) {
        thickness_profile[row] = dx;
      }
      if (num_profile[row] == 0) {
        packing_frac[1] += 1;
      }
      num_profile[row]++;
    }
  }
  packing_frac[0] /= nrow;
  packing_frac[1] /= nrow;
}

/* Export the data about wetting transition */
class WetDataExporter : public BaseExporter_2 {
public:
  WetDataExporter(const cmdline::parser &cmd);
  ~WetDataExporter();

  void dump_frame(int i_step,
                  const std::vector<float> &thickness_profile,
                  const std::vector<unsigned short> &num_profile,
                  const std::vector<double> &packing_frac);

  template <typename TPar, typename TBc>
  void write_frame(int i_step, const std::vector<TPar> &p_arr,
                   const std::vector<char> &flag_wetting, const TBc &bc);

private:
  void set_chunk_and_deflate();
  int deflate_level_;

  int ncid_;
  int row_id_;
  int time_id_;
  int wetting_frac_id_;
  int thickness_profile_id_;
  int num_profile_id_;

  size_t row_len_;
  size_t frame_len_;
  size_t time_idx_[1];
};


template<typename TPar, typename TBc>
void WetDataExporter::write_frame(int i_step, const std::vector<TPar>& p_arr,
                                  const std::vector<char>& flag_wetting,
                                  const TBc & bc) {
  std::vector<float> thickness_profile(row_len_ * 2, 0);
  std::vector<unsigned short> num_profile(row_len_ * 2, 0);
  std::vector<double> packing_frac(2, 0);
  cal_wetting_profile(thickness_profile, num_profile, packing_frac,
                      p_arr, flag_wetting, bc);
  dump_frame(i_step, thickness_profile, num_profile, packing_frac);
}
#endif
