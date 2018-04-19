#ifndef SINGLE_DOMAIN_2_H
#define SINGLE_DOMAIN_2_H
#include <functional>
#include "cellList2D.h"
#include "rand.h"
#include "particle.h"
#include "integrate.h"
#include "force.h"
#include "netCDFExporter.h"
#include "wetting.h"

template <typename TNode, typename TBc>
class Single_domain_2 {
public:
  typedef unsigned long long int ull;
  explicit Single_domain_2(const cmdline::parser &cmd);
  void ini_rand(double sigma = 1);

  template <typename TForce>
  void cal_force(const TForce &force);

  template <typename TIntegrate>
  void integrate(const TIntegrate &integ);

  template <typename TIntegrate>
  void integrate2(const TIntegrate &integ);

  template <typename TIntegrate>
  void integrate3(const TIntegrate &integ);

  void run(const cmdline::parser &cmd);
private:
  Ran myran_;
  TBc bc_;
  Cell_list_2<TNode> cell_list_;
  std::vector<TNode> p_arr_;

  int n_par_;
  double packing_frac_;
};

template<typename TNode, typename TBc>
Single_domain_2<TNode, TBc>::Single_domain_2(const cmdline::parser & cmd):
    myran_(cmd.get<unsigned long long>("seed")), bc_(cmd),
    cell_list_(bc_, cmd.get<double>("sigma")) {
  packing_frac_ = cmd.get<double>("phi");
  const auto lx = cmd.get<double>("Lx");
  const auto ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : lx;
  n_par_ = cal_particle_number_2(packing_frac_, lx, ly, 1);
  p_arr_.reserve(n_par_);
}

template<typename TNode, typename TBc>
void Single_domain_2<TNode, TBc>::ini_rand(double sigma) {
  create_rand_2(p_arr_, n_par_, sigma, myran_, bc_);
  cell_list_.create(p_arr_);
}

template<typename TNode, typename TBc>
template<typename TForce>
void Single_domain_2<TNode, TBc>::cal_force(const TForce &force) {
  auto f_ij = [this, &force](TNode *pi, TNode *pj) {
    force(*pi, *pj, bc_);
  };
  cell_list_.for_each_pair(f_ij);
}

template<typename TNode, typename TBc>
template<typename TIntegrate>
void Single_domain_2<TNode, TBc>::integrate(const TIntegrate &integ) {
  for (int i = 0; i < n_par_; i++) {
    integ(p_arr_[i], bc_, myran_);
  }
  cell_list_.recreate(p_arr_);
}

template<typename TNode, typename TBc>
template<typename TIntegrate>
void Single_domain_2<TNode, TBc>::integrate2(const TIntegrate &integ) {
  auto move = [this, &integ](TNode *p) {
    integ(*p, bc_, myran_);
  };
  cell_list_.update(move);
}

template<typename TNode, typename TBc>
template<typename TIntegrate>
void Single_domain_2<TNode, TBc>::integrate3(const TIntegrate &integ) {
  auto move = [this, &integ](TNode *p) {
    integ(*p, bc_, myran_);
  };
  cell_list_.update_by_row(move);
}

template<typename TNode, typename TBc>
void Single_domain_2<TNode, TBc>::run(const cmdline::parser & cmd) {
  const Run_and_tumble move(cmd.get<double>("h"), cmd.get<double>("alpha"));
  const SpringForce f_spring(cmd.get<double>("spring_const"));
  const auto n_step = cmd.get<int>("n_step");
  std::function<void()> integ;
  const auto mode = cmd.get<int>("int_mode");
  if (mode == 0)
    integ = [this, &move]() {integrate(move); };
  else if (mode == 1)
    integ = [this, &move]() {integrate2(move); };
  else
    integ = [this, &move]() {integrate3(move); };

  std::vector<std::ofstream> fout;
  LogExporter_2 *log = nullptr;
  XyExporter *xy = nullptr;
  NcParExporter_2 *nc = nullptr;
  WetDataExporter *profile = nullptr;
  const auto output_on = cmd.exist("output") ? true: false;
  if (output_on) {
    xy = new XyExporter(cmd);
    log = new LogExporter_2(cmd);
    nc = new NcParExporter_2(cmd);
    profile = new WetDataExporter(cmd);
  }
  double eps = 1.1;
  int min_pts = 2;
  double height_min = 0.55;
  const auto t1 = std::chrono::system_clock::now();
  for (int i = 1; i <= n_step; i++) {
    cal_force(f_spring);
    integ();
    if (output_on) {
      log->record(i);
      xy->write_frame(i, p_arr_);
      std::vector<bool> flag_clustered(n_par_);
      std::vector<char> flag_wetting(n_par_, 'V');
      std::vector<Cluster_w_xlim> c_arr;
      if (profile->need_export(i)) {
        dbscan_wall(c_arr, flag_clustered, flag_wetting,
                    eps, min_pts, height_min, p_arr_, bc_);
        profile->write_frame(i, p_arr_, flag_wetting, bc_);
      }
      if (nc->need_export(i)) {
        if (flag_clustered.size() == 0) {
          dbscan_wall(c_arr, flag_clustered, flag_wetting,
                      eps, min_pts, height_min, p_arr_, bc_);
        }
        nc->write_frame(i, p_arr_, flag_wetting);
      }
    }
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
  std::cout << std::endl;

  delete log;
  delete xy;
  delete nc;
  delete profile;
}

#endif
