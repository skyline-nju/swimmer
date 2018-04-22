#ifndef UNIDOMAIN_H
#define UNIDOMAIN_H
#include "rand.h"
#include "cmdline.h"
#include "cellList.h"
#include "force.h"
#include "integrate.h"
#include "wet_profile.h"
#include "particle.h"
#include "netCDFExporter.h"
#include <functional>

class UniDomain_2 {
public:
  typedef unsigned long long ull;

  explicit UniDomain_2(const cmdline::parser &cmd);

  template <typename TPar>
  void ini_rand(std::vector<TPar> &p_arr, CellListNode_2<TPar> **cl);

  template <typename TPar, typename BiFunc>
  void output(int i_step, const std::vector<TPar> &p_arr,
              BiFunc f_dis_square);

  template <typename TPar>
  void run(const cmdline::parser &cmd, std::vector<TPar> &p_arr,
           CellListNode_2<TPar> *cl);
private:
  Vec_2<int> flag_bc_;
  Vec_2<double> l_;
  Vec_2<double> half_l_;

  double sigma_;
  int n_par_;
  Ran myran_;

  // exporters
  LogExporter_2 *log_;
  XyExporter *xy_;
  NcParExporter_2 *nc_;
  ProfileExporter *profile_;
};

template <typename TPar>
void UniDomain_2::ini_rand(std::vector<TPar>& p_arr,
                           CellListNode_2<TPar>** cl) {
  p_arr.reserve(n_par_);
  Vec_2<double> l(l_);
  Vec_2<double> origin;
  if (flag_bc_.x == 2) {
    l.x -= sigma_;
    origin.x += 0.5 * sigma_;
  }
  if (flag_bc_.y == 2) {
    l.y -= sigma_;
    origin.y += 0.5 * sigma_;
  }
  create_rand_2(p_arr, n_par_, sigma_, myran_, l, origin, flag_bc_);
  *cl = new CellListNode_2<TPar>(l_, sigma_);
  (*cl)->create(p_arr);
}

template <typename TPar, typename BiFunc>
void UniDomain_2::output(int i_step, const std::vector<TPar>& p_arr,
                         BiFunc f_dis_square) {
  if (log_)
    log_->record(i_step);
  if (profile_) {
    std::vector<Cluster_w_xlim> c_arr;
    std::vector<bool> flag_clustered;
    std::vector<char> flag_wetting;
    auto cal_cluster = [this, f_dis_square, &p_arr](std::vector<Cluster_w_xlim> &cluster,
                                                    std::vector<bool> &flag_c,
                                                    std::vector<char> &flag_w) {
      dbscan_wall(cluster, flag_c, flag_w, profile_->get_eps(),
                  profile_->get_min_pts(), profile_->get_height_min(), p_arr,
                  f_dis_square, l_, Vec_2<double>());
    };
    if (profile_->need_export(i_step)) {
      cal_cluster(c_arr, flag_clustered, flag_wetting);
      auto lambda = [this, &flag_wetting, &p_arr](std::vector<float> &thickness,
                                                  std::vector<unsigned short> &num,
                                                  std::vector<double> &frac) {
        cal_wetting_profile(thickness, num, frac, p_arr, flag_wetting, l_, Vec_2<double>());
      };
      profile_->write_frame(i_step, lambda);
    }
    if (xy_ && xy_->need_export(i_step)) {
      if (c_arr.empty())
        cal_cluster(c_arr, flag_clustered, flag_wetting);
      xy_->write_frame(i_step, p_arr, flag_wetting);
    }
    if (nc_ && nc_->need_export(i_step)) {
      if (c_arr.empty())
        cal_cluster(c_arr, flag_clustered, flag_wetting);
      nc_->write_frame(i_step, p_arr, flag_wetting);
    }
  } else {
    if (xy_ && xy_->need_export(i_step))
      xy_->write_frame(i_step, p_arr);
    if (nc_ && nc_->need_export(i_step))
      nc_->write_frame(i_step, p_arr);
  }
}

template <typename TPar>
void UniDomain_2::run(const cmdline::parser& cmd, std::vector<TPar>& p_arr,
                      CellListNode_2<TPar>* cl) {
  const Run_and_tumble move(cmd.get<double>("h"), cmd.get<double>("alpha"));
  const SpringForce f_spring(cmd.get<double>("spring_const"));
  const auto n_step = cmd.get<int>("n_step");

  auto wrap = [this](TPar &p) {
    wrap_1(p.y, 0, l_.y, l_.y);
  };

  auto unwrap = [this](Vec_2<double> &vec) {
    unwrap_1(vec.y, l_.y, half_l_.y);
  };

  auto f_dis_square = [this](const TPar &pi, const TPar &pj) {
    const double dx = pi.x - pj.x;
    double dy = pi.y - pj.y;
    unwrap_1(dy, l_.y, half_l_.y);
    return dx * dx + dy * dy;
  };

  const auto k_wall = cmd.get<double>("k_wall");
  const auto xl = 0.5 * sigma_;
  const auto xr = l_.x - 0.5 * sigma_;
  auto wall_force = [this, xl, xr, k_wall](TPar &p) {
    if (p.x < xl) {
      p.fx += (xl - p.x) * k_wall;
    } else if (p.x > xr) {
      p.fx += (xr - p.x) * k_wall;
    }
  };

  auto integ = [this, &p_arr, cl, &move, wrap, wall_force] () {
    for (int i = 0; i < n_par_; i++) {
      move(p_arr[i], myran_, wrap, wall_force);
    }
    cl->recreate(p_arr);
  };

  auto pair_force = [&p_arr, cl, &f_spring, unwrap]() {
    cl->for_each_pair([&f_spring, unwrap](TPar *p1, TPar *p2) {
      f_spring.eval2(*p1, *p2, unwrap);
    });
  };

  const auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    pair_force();
    integ();
    output(i, p_arr, f_dis_square);
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}
#endif
