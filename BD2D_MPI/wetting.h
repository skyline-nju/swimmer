#ifndef WETTING_H
#define WETTING_H
#include "cellList2D.h"
#include "boundary.h"
#include <list>

template <typename TNode, typename TCellList, typename TBc>
void get_neighbors(int i, double eps_square, const std::vector<TNode> &p_arr,
    const TCellList &cl, const TBc &bc, std::list<int> &neighbor) {
  auto lambda = [eps_square, &p_arr, &bc, &neighbor](const TNode *pi, const TNode *pj) {
    if (get_dis_square(*pi, *pj, bc) < eps_square) {
      neighbor.push_back(pj - &p_arr[0]);
    }  
  };
  cl.for_nearby_par(&p_arr[i], lambda);
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
  void expand(double eps_square, int min_pts, const std::vector<TNode> &p_arr,
              const TCellList &cl, const TBc &bc, std::list<int> &neighbor,
              std::vector<bool> &flag_visited, std::vector<bool> &flag_clustered);

  template <typename UniFunc>
  void for_each(UniFunc f_i) const {
    const auto end = idx_arr_.cend();
    for (auto it = idx_arr_.cbegin(); it != end; ++it)
      f_i(*it);
  }

  int get_n() const {return idx_arr_.size();}

protected:
  std::vector<int> idx_arr_;
};

class Cluster_w_xlim : public Cluster {
public:
  Cluster_w_xlim() {x_min_=x_max_=0;}

  explicit Cluster_w_xlim(int idx, int size, double x):
      Cluster(idx, size), x_min_(x), x_max_(x) { }

  template <typename TNode, typename TCellList, typename TBc>
  void expand(double eps_square, int min_pts, const std::vector<TNode> &p_arr,
              const TCellList &cl, const TBc &bc, std::list<int> &neighbor,
              std::vector<bool> &flag_visited, std::vector<bool> &flag_clustered);
  double get_x_min() const {return x_min_;}
  double get_x_max() const {return x_max_;}
protected:
  double x_min_;
  double x_max_;
};

template <typename TNode, typename TCellList, typename TBc>
void Cluster::expand(double eps_square, int min_pts,
                     const std::vector<TNode>& p_arr, const TCellList& cl,
                     const TBc& bc, std::list<int> &neighbor,
                     std::vector<bool>& flag_visited, std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end();++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      get_neighbors(k, eps_square, p_arr, cl, bc, new_neighbor_idx);
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


 template <typename TNode, typename TCellList, typename TBc>
void Cluster_w_xlim::expand(double eps_square, int min_pts,
                            const std::vector<TNode>& p_arr, const TCellList& cl,
                            const TBc& bc, std::list<int>& neighbor,
                            std::vector<bool>& flag_visited, std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end();++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      get_neighbors(k, eps_square, p_arr, cl, bc, new_neighbor_idx);
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

template <typename TNode, typename TCellList, typename TBc, typename TCluster>
void dbscan(std::vector<TCluster> &c_arr,
            std::vector<bool> &flag_clustered,
            double eps, int min_pts,
            const std::vector<TNode> &p_arr,
            const TBc &bc,
            const TCellList &cl) {
  const auto eps_square = eps * eps;
  const auto n_par = p_arr.size();
  c_arr.reserve(n_par);
  std::vector<bool> flag_visited(n_par, false);
  for (int i = 0; i < n_par; i++) {
    if (!flag_visited[i]) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      get_neighbors(i, eps_square, p_arr, cl, bc, neighbor);
      if (neighbor.size()>=min_pts) {
        c_arr.emplace_back(i, n_par, p_arr[i].x);
        flag_clustered[i] = true;
        c_arr.back().expand(eps_square, min_pts, p_arr, cl, bc,
                            neighbor, flag_visited, flag_clustered);
      }
    }
  }
}

#endif
