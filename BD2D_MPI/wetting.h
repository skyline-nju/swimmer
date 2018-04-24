/**
 * @brief Cal wetting profile, packing fraction, correlation length, etc...
 * 
 * Study the wetting transition of active Brownian colloids / run-and-tumble swimmers.
 * 
 * @file wetting.h
 * @author skyline-nju
 * @date 2018-04-24
 */
#ifndef WETTING_H
#define WETTING_H
#include "cellList.h"
#include <list>
#include <forward_list>

/**
 * \brief Find neighbors of particle i.
 * \tparam TPar         Particle template
 * \tparam TCellList    Cell list template like CellListIdx_2<Tcontainer>
 * \tparam BiFunc       Template binary function.
 * \param i             Index of the target particle.
 * \param eps_square    Square of epsilon.
 * \param p_arr         Particle array.
 * \param cl            Cell list.
 * \param f_dis_square  Function to calculate square distance of two particles.
 * \param neighbor      Neighbor list of i particle.
 */
template <typename TPar, typename TCellList, typename BiFunc>
void find_neighbors(int i, double eps_square, const std::vector<TPar> &p_arr,
                    const TCellList &cl, BiFunc f_dis_square,
                    std::list<int> &neighbor) {
  const auto p0 = &p_arr[0];
  auto add_neighbor = [=, &neighbor](const TPar *pi, const TPar *pj) {
    if (f_dis_square(*pi, *pj) < eps_square) {
      neighbor.push_back(pj - p0);
    }
  };
  cl.for_nearby_par(&p_arr[i], p_arr, add_neighbor);
}

/************************************************************************//**
 * \brief Class for cluster
 ***************************************************************************/
class Cluster {
public:
  Cluster() = default;

  explicit Cluster(int idx, int size) {
    idx_arr_.reserve(size);
    idx_arr_.push_back(idx);
  }

  explicit Cluster(int idx, int size, double x) {
    idx_arr_.reserve(size);
    idx_arr_.push_back(idx);
  }

  /**
   * \brief Expand the cluster.
   * \tparam TPar          Template for particle.
   * \tparam TCellList     Template for cell list.
   * \tparam BiFunc        Template for binary function
   * \param eps_square     Square of epsilon.
   * \param min_pts        Min neighbors that a core particle should have at least.
   * \param p_arr          Array of particls.
   * \param cl             Cell list.
   * \param f_dis_square   Function to calculate the square distance between two particles.
   * \param neighbor       Neighbor list.
   * \param flag_visited   Array to record whether one particle has been visited.
   * \param flag_clustered Array to record whether one particle have been clustered.
   */
  template <typename TPar, typename TCellList, typename BiFunc>
  void expand(double eps_square, int min_pts,
              const std::vector<TPar>& p_arr,
              const TCellList &cl,
              BiFunc f_dis_square,
              std::list<int> &neighbor,
              std::vector<bool> &flag_visited,
              std::vector<bool> &flag_clustered);
  /**
   * \brief Visit each particle in the cluster.
   * \tparam UniFunc Template for unitary function
   * \param f_i      A Unitary function.
   */
  template <typename UniFunc>
  void for_each(UniFunc f_i) const;

  /**
   * \brief Get the number of particle in the cluster.
   * \return Number of particle in the cluster.
   */
  int get_n() const { return idx_arr_.size(); } 

  /**
   * \brief Mark the cluster particles with a given flag.
   * \tparam T        Template for the type of flag.
   * \param flag_arr  Array of flag.
   * \param flag      Value of the flag to mark.
   */
  template <typename T>
  void mark(std::vector<T> &flag_arr, T flag) const {
    for_each([flag, &flag_arr](int i) {flag_arr[i] = flag; });
  }

protected:
  std::vector<int> idx_arr_;  //!< index of particls of this cluster.
};

/*************************************************************************//**
 * \brief Cluster with x_min and x_max.
 ****************************************************************************/
class Cluster_w_xlim : public Cluster {
public:
  Cluster_w_xlim() { x_min_ = x_max_ = 0; }

  explicit Cluster_w_xlim(int idx, int size, double x) :
    Cluster(idx, size), x_min_(x), x_max_(x) {}

  /**
   * \brief Expand the cluster and update the x_min and x_max of the cluster.
   * \tparam TPar          Template for particle.
   * \tparam TCellList     Template for cell list.
   * \tparam BiFunc        Template for binary function
   * \param eps_square     Square of epsilon.
   * \param min_pts        Min neighbors that a core particle should have at least.
   * \param p_arr          Array of particls.
   * \param cl             Cell list.
   * \param f_dis_square   Function to calculate the square distance between two particles.
   * \param neighbor       Neighbor list.
   * \param flag_visited   Array to record whether one particle has been visited.
   * \param flag_clustered Array to record whether one particle have been clustered.
   */
  template <typename TPar, typename TCellList, typename BiFunc>
  void expand(double eps_square, int min_pts,
              const std::vector<TPar> &p_arr,
              const TCellList &cl,
              BiFunc f_dis_square,
              std::list<int> &neighbor,
              std::vector<bool> &flag_visited,
              std::vector<bool> &flag_clustered);

  /**
   * \brief Mark the cluster particle as 'L', 'V' or "R".
   * 
   * 'V' for particles in vapor states, while 'L', 'R' denote the particles are
   * wetted on the left and right walls, respectively.
   * \param x_left  Leftmost threshold for wetting.
   * \param x_right Rightmost threshold for wetting.
   * \param flag_w  Array of flags for wetting.
   */
  void mark_wetting(double x_left, double x_right,
                    std::vector<char>& flag_w) const;
protected:
  double x_min_;  //!< Min x of the cluster.
  double x_max_;  //!< Max x of the cluster.
};

inline void Cluster_w_xlim::mark_wetting(double x_left, double x_right,
                                         std::vector<char>& flag_w) const {
  if (x_min_ < x_left)
    mark(flag_w, 'L');
  else if (x_max_ >= x_right)
    mark(flag_w, 'R');
}

template <typename TPar, typename TCellList, typename BiFunc>
void Cluster::expand(double eps_square, int min_pts,
                     const std::vector<TPar>& p_arr,
                     const TCellList& cl,
                     BiFunc f_dis_square,
                     std::list<int> &neighbor,
                     std::vector<bool>& flag_visited,
                     std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end(); ++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      find_neighbors(k, eps_square, p_arr, cl, f_dis_square, neighbor);
      if (new_neighbor_idx.size() >= min_pts) {
        neighbor.splice(neighbor.end(), new_neighbor_idx);
      }
    }
    if (!flag_clustered[k]) {
      flag_clustered[k] = true;
      idx_arr_.push_back(k);
    }
  }
  idx_arr_.reserve(idx_arr_.size());
}

template <typename UniFunc>
void Cluster::for_each(UniFunc f_i) const {
  const auto end = idx_arr_.cend();
  for (auto it = idx_arr_.cbegin(); it != end; ++it)
    f_i(*it);
}

template <typename TNode, typename TCellList, typename BiFunc>
void Cluster_w_xlim::expand(double eps_square, int min_pts,
                            const std::vector<TNode>& p_arr,
                            const TCellList& cl,
                            BiFunc f_dis_square,
                            std::list<int>& neighbor,
                            std::vector<bool>& flag_visited,
                            std::vector<bool>& flag_clustered) {
  for (auto it = neighbor.begin(); it != neighbor.end(); ++it) {
    auto k = *it;
    if (!flag_visited[k]) {
      flag_visited[k] = true;
      std::list<int> new_neighbor_idx;
      find_neighbors(k, eps_square, p_arr, cl, f_dis_square, neighbor);
      if (new_neighbor_idx.size() >= min_pts) {
        neighbor.splice(neighbor.end(), new_neighbor_idx);
      }
    }
    if (!flag_clustered[k]) {
      flag_clustered[k] = true;
      idx_arr_.push_back(k);
      //! only valid for non-periodic boundary condition in the x direction
      if (p_arr[k].x < x_min_) {
        x_min_ = p_arr[k].x;
      }
      if (p_arr[k].x > x_max_) {
        x_max_ = p_arr[k].x;
      }
    }
  }
  idx_arr_.reserve(idx_arr_.size());
}

/**
 * \brief DBSCAN for the periodic boundary condition.
 * \tparam TPar           Template for particle
 * \tparam BiFunc         Template binary function
 * \tparam TCluster       Template for cluster
 * \param c_arr           Array of cluster
 * \param flag_clustered  Array of flag for clustering.
 * \param eps             Threshold distance between two neighbors.
 * \param min_pts         Min neighbors that a core particle should have.
 * \param p_arr           Particle array.
 * \param f_dis_square    Func to cal the square distance of two particles.
 * \param l               Domain lengths (Vec_2)
 * \param origin          Origin point of the damain.
 */
template <typename TPar, typename BiFunc, typename TCluster>
void dbscan(std::vector<TCluster> &c_arr,
            std::vector<bool> &flag_clustered,
            double eps, unsigned int min_pts,
            const std::vector<TPar> &p_arr,
            BiFunc f_dis_square,
            const Vec_2<double> &l,
            const Vec_2<double> &origin) {
  const auto eps_square = eps * eps;
  const auto n_par = p_arr.size();
  c_arr.reserve(n_par);
  flag_clustered.reserve(n_par);
  for (size_t i = 0; i < n_par; i++)
    flag_clustered.push_back(false);
  std::vector<bool> flag_visited(n_par, false);
  CellListIdx_2<std::forward_list<int>> cl(l, eps, origin);
  cl.create(p_arr);
  for (unsigned int i = 0; i < n_par; i++) {
    if (!flag_visited[i]) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      find_neighbors(i, eps_square, p_arr, cl, f_dis_square, neighbor);
      if (neighbor.size() >= min_pts) {
        c_arr.emplace_back(i, n_par, p_arr[i].x);
        flag_clustered[i] = true;
        c_arr.back().expand(eps_square, min_pts, p_arr, cl, f_dis_square,
                            neighbor, flag_visited, flag_clustered);
      }
    }
  }
}

/**
 * \brief DBSCAN for two walls at x=0 and x=L.
 * \tparam TPar               Template for particle
 * \tparam BiFunc             Template binary function
 * \param c_arr               Array of Cluster_w_xlim 
 * \param flag_clustered      Array of flag for clustering.
 * \param flag_wetting        Array of flag for wetting.
 * \param eps                 Threshold distance between two neighbors.
 * \param min_pts             Min neighbors that a core particle should have.
 * \param height_min          Threshold distance from the particle to the wall.
 * \param p_arr               Particle array.
 * \param f_dis_square        Func to cal the square distance of two particles.
 * \param l                   Domain lengths (Vec_2)
 * \param origin              Origin point of the damain.
 */
template <typename TPar, typename BiFunc>
void dbscan_wall(std::vector<Cluster_w_xlim> &c_arr,
                 std::vector<bool> &flag_clustered,
                 std::vector<char> &flag_wetting,
                 double eps, unsigned int min_pts, double height_min,
                 const std::vector<TPar> &p_arr, BiFunc f_dis_square,
                 const Vec_2<double> &l, const Vec_2<double> &origin) {
  const auto eps_square = eps * eps;
  const auto x_left = height_min;
  const auto x_right = origin.x + l.x - height_min;
  const auto n_par = p_arr.size();
  c_arr.reserve(n_par);
  flag_clustered.reserve(n_par);
  flag_wetting.reserve(n_par);
  for (size_t i = 0; i < n_par; i++) {
    flag_clustered.push_back(false);
    flag_wetting.push_back('V');
  }
  std::vector<bool> flag_visited(n_par, false);
  CellListIdx_2<std::forward_list<int>> cl(l, eps, origin);
  cl.create(p_arr);
  for (unsigned int i = 0; i < n_par; i++) {
    if (!flag_visited[i]) {
      flag_visited[i] = true;
      std::list<int> neighbor;
      find_neighbors(i, eps_square, p_arr, cl, f_dis_square, neighbor);
      if (neighbor.size() >= min_pts) {
        c_arr.emplace_back(i, n_par, p_arr[i].x);
        flag_clustered[i] = true;
        c_arr.back().expand(eps_square, min_pts, p_arr, cl, f_dis_square,
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

/**
 * \brief Get the size of each cluster.
 * \tparam TCluster Tempalte for cluster.
 * \param cluster   Array of clusters.
 * \param size_arr  Array of cluster sizes.
 */
template <typename TCluster>
void get_cluster_size(const std::vector<TCluster> &cluster,
                      std::vector<int> &size_arr) {
  const auto end = cluster.cend();
  for (auto it = cluster.cbegin(); it != end; ++it) {
    (*it).mark(size_arr, (*it).get_n());
  }
}


/**
 * \brief Cal the wetting profile.
 * \tparam TPar               Template for particle
 * \param thickness_profile   Thickness profile for wetting layer.
 * \param num_profile         Number profile for wetting layer.
 * \param packing_frac        Packing fraction for wetting layer.
 * \param p_arr               Particle array.
 * \param flag_wetting        Array of flags for wetting.
 * \param l                   Domain lengths.
 * \param origin              Origin point of the domain.
 */
template <typename TPar>
void cal_wetting_profile(std::vector<float> &thickness_profile,
                         std::vector<unsigned short> &num_profile,
                         std::vector<double> &packing_frac,
                         const std::vector<TPar> &p_arr,
                         const std::vector<char> &flag_wetting,
                         const Vec_2<double> &l,
                         const Vec_2<double> &origin) {
  const auto nrow = int(l.y);
  const auto x_min = origin.x;
  const auto x_max = origin.x + l.x;
  const auto y_min = origin.y;
  const auto n = p_arr.size();
  for (unsigned int i = 0; i < n; i++) {
    if (flag_wetting[i] == 'L') {
      const auto row = int(p_arr[i].y - y_min);
      const double dx = p_arr[i].x - x_min;
      if (thickness_profile[row] < dx) {
        thickness_profile[row] = dx;
      }
      if (num_profile[row] == 0) {
        packing_frac[0] += 1;
      }
      num_profile[row]++;
    } else if (flag_wetting[i] == 'R') {
      const auto row = int(p_arr[i].y - y_min) + nrow;
      const double dx = x_max - p_arr[i].x;
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

#endif
