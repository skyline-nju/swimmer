#ifndef UNIDOMAIN_H
#define UNIDOMAIN_H
#include <functional>
#include "rand.h"
#include "cmdline.h"
#include "cellList.h"
#include "force.h"
#include "integrate.h"
#include "wetting.h"
#include "particle.h"
#include "netCDFExporter.h"

class UniDomain_2 {
public:
  typedef unsigned long long ull;

  explicit UniDomain_2(const cmdline::parser &cmd);

  ~UniDomain_2();

  template <typename TPar>
  void ini_rand(std::vector<TPar> &p_arr, CellListNode_2<TPar> **cl);

  template <typename TPar, typename BiFunc>
  void output(int i_step, const std::vector<TPar> &p_arr,
              BiFunc f_dis_square);

  template <typename TPar>
  void run(const cmdline::parser &cmd, std::vector<TPar> &p_arr,
           CellListNode_2<TPar> *cl);
private:
  Vec_2<int> flag_bc_;       //!< boundary condition flag
  Vec_2<double> l_;          //!< domain lengths
  Vec_2<double> half_l_;     //!< half of domain lengths

  double sigma_;             //!< particle diameter
  int n_par_;                //!< particle number
  Ran myran_;                //!< random number generator

  LogExporter_2 *log_;       //!< log exporter
  XyExporter *xy_;           //!< snapshot exporter with xyz format
  NcParExporter_2 *nc_;      //!< snapshto exporter with netcdf format
  ProfileExporter *profile_; //!< profile exporter
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
    auto cal_cluster = [this, f_dis_square, &p_arr](auto &cluster, auto &flag_c, auto &flag_w) {
      dbscan_wall(cluster, flag_c, flag_w, profile_->get_eps(),
                  profile_->get_min_pts(), profile_->get_height_min(), p_arr,
                  f_dis_square, l_, Vec_2<double>());
    };
    if (profile_->need_export(i_step)) {
      cal_cluster(c_arr, flag_clustered, flag_wetting);
      auto lambda = [this, &flag_wetting, &p_arr](auto &thickness, auto &num, auto &frac) {
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
  const Run_and_tumble run_tumble(cmd.get<double>("h"), cmd.get<double>("alpha"));
  const SpringForce f_spring(cmd.get<double>("spring_const"));
  const auto n_step = cmd.get<int>("n_step");

  auto tangle = [this](TPar &p) {
    tangle_1(p.y, 0, l_.y, l_.y);
  };

  auto untangle = [this](Vec_2<double> &vec) {
    untangle_1(vec.y, l_.y, half_l_.y);
  };

  auto cal_dis_square = [this](const TPar &pi, const TPar &pj) {
    const double dx = pi.x - pj.x;
    double dy = pi.y - pj.y;
    untangle_1(dy, l_.y, half_l_.y);
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

  auto integ = [this, &p_arr, cl, &run_tumble, tangle, wall_force] () {
    for (int i = 0; i < n_par_; i++) {
      run_tumble(p_arr[i], myran_, tangle, wall_force);
    }
    cl->recreate(p_arr);
  };

  auto pair_force = [&p_arr, cl, &f_spring, untangle]() {
    cl->for_each_pair([&f_spring, untangle](TPar *p1, TPar *p2) {
      f_spring.eval(*p1, *p2, untangle);
    });
  };

  const auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    pair_force();
    integ();
    output(i, p_arr, cal_dis_square);
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}
#endif
