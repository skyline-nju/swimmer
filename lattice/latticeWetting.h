#ifndef LATTICEWETTING_H
#define LATTICEWETTING_H
#include <list>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include "vect.h"
namespace lattice {

void find_neighbors(int ic, std::list<int> &neighbor,
                    const Vec_2<int> &l_domain,
                    const std::vector<uint8_t> &cell,
                    const Vec_2<bool> &flag_bc);

class Cluster {
public:
  explicit Cluster(int ic, int size);

  void expand(const std::vector<uint8_t>& cell,
              std::list<int> &neighbor,
              std::vector<bool> &flag_visited,
              std::vector<bool> &flag_clustered);

  void check_wetting(std::vector<char> &flag_w) const;

  static void find_clusters(std::vector<Cluster> &c_arr,
                            const std::vector<uint8_t> &cell);

  static void find_clusters(std::vector<Cluster> &c_arr,
                            const std::vector<uint8_t> &cell,
                            std::vector<char> &flag_wetting);

  template <typename TPar, typename T1, typename T2>
  static void cal_wetting_profile(const std::vector<uint8_t> &cell,
                                  const std::vector<TPar> &p_arr,
                                  std::vector<char> &type_arr,
                                  std::vector<T1> &thickness_profile,
                                  std::vector<T2> &num_profile,
                                  std::vector<double> &pack_frac);

  void push_back(int ic);

  static Vec_2<int> domain_len;
protected:
  std::vector<int> idx_cell_;
  int x_min_;
  int x_max_;
  int y_min_;
  int y_max_;

};

void cal_non_empty_frac(const std::vector<uint8_t> &cell,
                        std::vector<double> &packing_frac);

template <typename T1, typename T2>
void cal_non_empty_profile(const std::vector<uint8_t> &cell,
                           std::vector<T1> &thickness_profile,
                           std::vector<T2> &num_profile,
                           std::vector<double> &packing_frac) {
  thickness_profile.assign(Cluster::domain_len.y * 2, 0);
  num_profile.assign(Cluster::domain_len.y * 2, 0);
  packing_frac.assign(2, 0);
  for (auto row = 0; row < Cluster::domain_len.y; row++) {
    const size_t row_lx = row * Cluster::domain_len.x;
    const auto i_lt = row * 2;
    const auto i_rt = i_lt + 1;
    for (auto col = 0; col < Cluster::domain_len.x; col++) {
      const auto ic = col + row_lx;
      if (cell[ic]) {
        ++thickness_profile[i_lt];
        num_profile[i_lt] += cell[ic];
      } else {
        if (thickness_profile[i_lt]) {
          packing_frac[0] += 1;
        }
        break;
      }
    }
    for (auto col = Cluster::domain_len.x - 1; col >= 0; col--) {
      const auto ic = col + row_lx;
      if (cell[ic]) {
        ++thickness_profile[i_rt];
        num_profile[i_rt] += cell[ic];
      } else {
        if (thickness_profile[i_rt]) {
          packing_frac[1] += 1;
        }
        break;
      }
    }
  }
  packing_frac[0] /= Cluster::domain_len.y;
  packing_frac[1] /= Cluster::domain_len.y;
}

template <typename TPar, typename T1, typename T2>
void Cluster::cal_wetting_profile(const std::vector<uint8_t>& cell,
                                  const std::vector<TPar>& p_arr,
                                  std::vector<char> &type_arr,
                                  std::vector<T1>& thickness_profile,
                                  std::vector<T2>& num_profile,
                                  std::vector<double>& pack_frac) {
  const auto n_par = p_arr.size();

  std::vector<Cluster> cluster_arr;
  std::vector<char> wetting_cell(cell.size(), 'V');
  find_clusters(cluster_arr, cell, wetting_cell);

  type_arr.reserve(n_par);
  for (size_t i = 0; i < n_par; i++) {
    auto ic = p_arr[i].x + domain_len.x * p_arr[i].y;
    type_arr.push_back(wetting_cell[ic]);
  }

  thickness_profile.assign(domain_len.y * 2, 0);
  num_profile.assign(domain_len.y * 2, 0);
  pack_frac.assign(4, 0);

  for (auto row = 0; row < domain_len.y; row++) {
    const size_t row_lx = row * domain_len.x;
    const auto i_lt = row * 2;
    const auto i_rt = i_lt + 1;
    for (auto col = 0; col < domain_len.x; col++) {
      const auto ic = col + row_lx;
      const auto cell_type = wetting_cell[ic];
      if (cell_type == 'L') {
        thickness_profile[i_lt] = col + 1;
        num_profile[i_lt] += cell[ic];
      } else if (cell_type == 'R') {
        thickness_profile[i_rt] = domain_len.x - i_rt;
        num_profile[i_rt] += cell[ic];
      }
    }
  }

  for (int i = 0; i < domain_len.y; i++) {
    if (thickness_profile[i * 2]) {
      pack_frac[0] += 1;
    }
    if (thickness_profile[i * 2 + 1]) {
      pack_frac[1] += 1;
    }
  }
  pack_frac[0] /= domain_len.y;
  pack_frac[1] /= domain_len.y;

  std::vector<double> pack_frac2;
  cal_non_empty_frac(cell, pack_frac2);
  pack_frac[2] = pack_frac2[0];
  pack_frac[3] = pack_frac2[1];
}

}
#endif