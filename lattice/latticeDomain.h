/**
 * @brief lattice domain
 * 
 * @file latticeDomain.h
 * @author skyline-nju
 * @date 2018-04-27
 */
#ifndef LATTICLEDOMAIN_H
#define LATTICLEDOMAIN_H

#include <cstdint>
#include <vector>
#include <functional>
#include "latticeExporter.h"
#include "vect.h"
#include "comn.h"
#include "cmdline.h"
#include "latticeParticle.h"

namespace lattice {

/**
 * \brief setting the sequence of particles in one mc move
 * 
 * \tparam TFunc1       std::function
 * \tparam TFunc2       Function template
 * \tparam TPar         Template for particles
 * \tparam TRan         Template for random number generator
 * \param mc_move       Functional for one mv move
 * \param trival_move   Functional for updating one particle
 * \param p_arr         Particle array
 * \param ptr_arr       Array of particle pointers
 * \param myran         Random number generator
 * \param mode          Squence to perform all trival moves.
 */
template <typename TFunc1, typename TFunc2, typename TPar, typename TRan>
void set_sampling_seque(TFunc1 &mc_move, TFunc2 trival_move,
                        std::vector<TPar> &p_arr, std::vector<TPar*> &ptr_arr,
                        TRan &myran, int mode) {
  switch (mode) {
  case(0): {
    std::cout << "choose particles randomly\n";
    mc_move = [&myran, &p_arr, trival_move]() {
      const unsigned int n = p_arr.size();
      for (unsigned int count = 0; count < n; count++) {
        auto i = int(myran.doub() * n);
        trival_move(p_arr[i], myran);
      }
    };
    break;
  }
  case(1): {
    std::cout << "generate a random squence in advance\n";
    mc_move = [&myran, &p_arr, trival_move]() {
      const unsigned int n = p_arr.size();
      auto idx = new int[n];
      for (unsigned int i = 0; i < n; i++)
        idx[i] = int(myran.doub() * n);
      for (unsigned int i = 0; i < n; i++)
        trival_move(p_arr[idx[i]], myran);
      delete[] idx;
    };
    break;
  }
  case(2): {
    std::cout << "shuffle the particles\n";
    mc_move = [&myran, &p_arr, trival_move]() {
      shuffle(p_arr, myran);
      for (auto &i : p_arr) {
        trival_move(i, myran);
      }
    };
    break;
  }
  case(3): {
    std::cout << "shuffle the pointers\n";
    ptr_arr.reserve(p_arr.size());
    for (auto &i : p_arr) {
      ptr_arr.push_back(&i);
    }
    mc_move = [&myran, &ptr_arr, trival_move]() {
      shuffle(ptr_arr, myran);
      for (auto &i : ptr_arr) {
        trival_move(*i, myran);
      }
    };
    break;
  }
  default: {
    std::cout << "Error, sequence mode should be one of 0, 1, 2, 3\n";
    exit(1);
  }
  }
}

/*************************************************************************//**
 * @brief Base class for 2D lattice domain
 * 
 ***************************************************************************/
class UniDomain_2 {
public:
  /**
   * @brief Construct a new UniDomain_2 object
   * 
   * @tparam TPar     Template of particle
   * @tparam TRan     Template of random number generator
   * @param cmd       cmdline parser
   * @param p_arr     particle array
   * @param myran     random number generator
   * @param is_rt     is run-and-tumble case?
   */
  template<typename TPar, typename TRan>
  explicit UniDomain_2(const cmdline::parser &cmd, std::vector<TPar> &p_arr,
                       TRan &myran, bool is_rt);
  /**
   * @brief Destroy the UniDomain_2 object
   * 
   */
  ~UniDomain_2() {delete log_; delete snap_; delete pf_;}

protected:
  template <typename TPar>
  size_t get_ic(const TPar &p) const { return p.x + p.y * l_.x; }

  size_t get_ic(int ix, int iy) const { return ix + iy * l_.x; }

  template <typename T>
  void tangle_x(T &x) const { tangle_1(x, 0, l_.x, l_.x); }

  template <typename T>
  void tangle_y(T &y) const { tangle_1(y, 0, l_.y, l_.y); }

  /**
   * @brief Base function to run the simulation.
   * 
   * @tparam TPar     Template of particle
   * @tparam TRan     Template of random number generator
   * @tparam TMove    Template function of trival moves
   * @param p_arr     Array of particles
   * @param myran     Random number generator
   * @param f_move    Trival moves
   * @param n_step    Total steps to run 
   * @param seq_mode  For each step, the sequences to update particles
   */
  template<typename TPar, typename TRan, typename TMove>
  void base_run(std::vector<TPar>& p_arr, TRan& myran, TMove f_move,
                int n_step, int seq_mode);

  Vec_2<int> l_;              //!> domain lengths
  std::vector<uint8_t> cell_; //!> number of particles at each cell
  const int max_capacity_;    //!> the max number of particle that one cell can hold

  LogExporter *log_;          //!> log exporter
  SnapExporter_2 * snap_;     //!> snapshot exporter
  ProfileExporter *pf_;       //!> profile exporter
};
/**
 * @brief Construct a new UniDomain_2 object
 * 
 * @tparam TPar     Template of particle
 * @tparam TRan     Template of random number generator
 * @param cmd       cmdline parser
 * @param p_arr     particle array
 * @param myran     random number generator
 * @param is_rt     is run-and-tumble case?
 */
template <typename TPar, typename TRan>
UniDomain_2::UniDomain_2(const cmdline::parser& cmd, std::vector<TPar>& p_arr,
                         TRan &myran, bool is_rt)
  : max_capacity_(cmd.get<int>("n_max")),
    log_(nullptr), snap_(nullptr), pf_(nullptr) {
  l_.x = cmd.get<int>("Lx");
  l_.y = cmd.exist("Ly") ? cmd.get<int>("Ly") : l_.x;
  const auto phi = cmd.get<double>("pack_frac");
  const size_t ncells = l_.x * l_.y;
  unsigned int n_par = int(ncells * phi);
  cell_.reserve(ncells);
  p_arr.reserve(n_par);
  for (size_t i = 0; i < ncells; i++) {
    cell_.push_back(0);
  }
  while (p_arr.size() < n_par) {
    TPar p(myran, l_);
    const auto ic = get_ic(p);
    if (cell_[ic] == 0) {
      p_arr.push_back(p);
      cell_[ic] += 1;
    }
  }
  set_output_2(cmd, &log_, &snap_, &pf_, is_rt);
}

/**
 * @brief Base function to run the simulation.
 * 
 * @tparam TPar     Template of particle
 * @tparam TRan     Template of random number generator
 * @tparam TMove    Template function of trival moves
 * @param p_arr     Array of particles
 * @param myran     Random number generator
 * @param f_move    Trival moves
 * @param n_step    Total steps to run 
 * @param seq_mode  For each step, the sequences to update particles
 */
template <typename TPar, typename TRan, typename TMove>
void UniDomain_2::base_run(std::vector<TPar>& p_arr, TRan& myran, TMove f_move,
                           int n_step, int seq_mode) {

  std::vector<TPar*> ptr_arr;
  std::function<void()> mc_move;
  set_sampling_seque(mc_move, f_move, p_arr, ptr_arr, myran, seq_mode);

  const auto t1 = std::chrono::system_clock::now();
  for (auto i = 1; i <= n_step; i++) {
    mc_move();
    output_2(i, p_arr, cell_, log_, snap_, pf_);
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}

/*************************************************************************//**
 * @brief 2D domain where lattice particles perform run-and-tumble motions
 * 
 ***************************************************************************/
class UniDomain_RT_2:public UniDomain_2 {
public:
  template <typename TPar, typename TRan>
  UniDomain_RT_2(const cmdline::parser &cmd, std::vector<TPar> &p_arr,
                 TRan &myran);

  template <typename TPar, typename TRan>
  void run(std::vector<TPar> &p_arr, TRan &myran, int n_step, int seq_mode);

  template <typename TPar, typename TRan>
  void run_wx(std::vector<TPar> &p_arr, TRan &myran, int n_step, int seq_mode);

private:
  template <typename TPar>
  void jump_forward(TPar &p);

  template <typename TPar>
  void jump_forward_wall_x(TPar &p);

  const double alpha_;  //!< tumbling rate  
};

/**
 * @brief Construct a new UniDomain_RT_2 object
 * 
 * @tparam TPar     Template of particle
 * @tparam TRan     Template of random number generator
 * @param cmd       cmdline parser
 * @param p_arr     particle array
 * @param myran     random number generator
 */
template<typename TPar, typename TRan>
UniDomain_RT_2::UniDomain_RT_2(const cmdline::parser &cmd,
                               std::vector<TPar>& p_arr,
                               TRan & myran)
  : UniDomain_2(cmd, p_arr, myran, true), alpha_(cmd.get<double>("alpha")) {
  std::cout << "Run-and-tumble lattice particle\n";
}
/**
* @brief Run simulation with periodic boundary condition
*
* @tparam TPar       Template of particle
* @tparam TRan       Template of random number generator
* @param p_arr       Array of particles
* @param myran       Random number generator
* @param n_step      Total steps to run
* @param seq_mode    Updating sequences for each time step
*/
template <typename TPar, typename TRan>
void UniDomain_RT_2::run(std::vector<TPar>& p_arr, TRan& myran,
                         int n_step, int seq_mode) {
  auto trival_move = [this](TPar &p, TRan &ran) {
    jump_forward(p);
    if (ran.doub() < alpha_)
      p.tumble(ran);
  };
  base_run(p_arr, myran, trival_move, n_step, seq_mode);
}
/**
* @brief Run simulation with two walls at xmin and xmax
*
* @tparam TPar       Template of particle
* @tparam TRan       Template of random number generator
* @param p_arr       Array of particles
* @param myran       Random number generator
* @param n_step      Total steps to run
* @param seq_mode    Updating sequences for each time step
*/
template<typename TPar, typename TRan>
void UniDomain_RT_2::run_wx(std::vector<TPar>& p_arr, TRan & myran,
                            int n_step, int seq_mode) {
  auto trival_move = [this](TPar &p, TRan &ran) {
    jump_forward_wall_x(p);
    if (ran.doub() < alpha_)
      p.tumble(ran);
  };
  base_run(p_arr, myran, trival_move, n_step, seq_mode);
}

/**
 * @brief Jump forward by one size under PBC
 * 
 * @tparam TPar  Template of particle
 * @param p      The particle to jump
 */
template <typename TPar>
void UniDomain_RT_2::jump_forward(TPar& p) {
  int x_new = p.x + p.get_ux();
  tangle_x(x_new);
  int y_new = p.y + p.get_uy();
  tangle_y(y_new);
  const auto ic_new = get_ic(x_new, y_new);
  if (cell_[ic_new] < max_capacity_) {
    const auto ic = get_ic(p);
    --cell_[ic];
    p.x = x_new;
    p.y = y_new;
    ++cell_[ic_new];
  }
}

/**
 * @brief Jump forward by one size under PBC
 * 
 * @tparam TPar  Template of particle
 * @param p      The particle to jump
 */
template <typename TPar>
void UniDomain_RT_2::jump_forward_wall_x(TPar& p) {
  const int x_new = p.x + p.get_ux();
  if (x_new >= 0 && x_new < l_.x) {
    int y_new = p.y + p.get_uy();
    tangle_y(y_new);
    const auto ic_new = get_ic(x_new, y_new);
    if (cell_[ic_new] < max_capacity_) {
      const auto ic = get_ic(p);
      --cell_[ic];
      p.x = x_new;
      p.y = y_new;
      ++cell_[ic_new];
    }
  }
}

/*************************************************************************//**
 * @brief 2D domain where lattice particles perform active Brownian motions
 * 
 ***************************************************************************/
class UniDomain_AB_2: public UniDomain_2 {
public:
  template <typename TPar, typename TRan>
  UniDomain_AB_2(const cmdline::parser &cmd, std::vector<TPar> &p_arr,
                 TRan &myran);

  template <typename TPar, typename TRan>
  void run(std::vector<TPar> &p_arr, TRan &myran, int n_step, int seq_mode);

  template <typename TPar, typename TRan>
  void run_wx(std::vector<TPar> &p_arr, TRan &myran, int n_step, int seq_mode);

protected:
  template <typename Tx, typename Tu>
  void rot_diffuse(Par_2<Tx, Tu> &p, double rand_value) const;

  template <typename Tx>
  void rot_diffuse(Par_s_2<Tx> &p, double rand_value) const;

  template <typename Tx, typename Tu>
  void trans_diffuse(Par_2<Tx, Tu> &p, double rand_value);

  template <typename Tx>
  void trans_diffuse(Par_s_2<Tx> &p, double rand_value);

  template <typename Tx, typename Tu>
  void trans_diffuse_wx(Par_2<Tx, Tu> &p, double rand_value);

  template <typename Tx>
  void trans_diffuse_wx(Par_s_2<Tx> &p, double rand_value);

private:
  const double D_rot_;
  const double D_rot_2_;
  double prob_slice_[3];
};

/**
 * @brief Construct a new UniDomain_AB_2 object
 * 
 * @tparam TPar     Template of particle
 * @tparam TRan     Template of random number generator
 * @param cmd       cmdline parser
 * @param p_arr     particle array
 * @param myran     random number generator
 */
template <typename TPar, typename TRan>
UniDomain_AB_2::UniDomain_AB_2(const cmdline::parser& cmd,
                               std::vector<TPar>& p_arr, TRan& myran)
  : UniDomain_2(cmd, p_arr, myran, false),
    D_rot_(cmd.get<double>("D_rot")), D_rot_2_(D_rot_ * 2), prob_slice_{} {
  const auto nu_f = cmd.get<double>("nu_f");
  const auto nu_b = cmd.get<double>("nu_b");
  const auto nu_t = cmd.get<double>("nu_t");
  const auto rate_tot = nu_f + nu_b + 2 * nu_t;
  prob_slice_[0] = nu_f / rate_tot;
  prob_slice_[1] = prob_slice_[0] + nu_b / rate_tot;
  prob_slice_[2] = prob_slice_[1] + nu_t / rate_tot;
  std::cout << "active Brownain on-lattice particles\n";
}
/**
 * @brief Rotation of on-lattice particles
 *
 * @tparam Tx         Template type of x, y
 * @tparam Tu         Template type of ux, uy
 * @param p           One particle
 * @param rand_value  Random number range from 0 to 1
 */
template <typename Tx, typename Tu>
void UniDomain_AB_2::rot_diffuse(Par_2<Tx, Tu>& p, double rand_value) const {
  auto tmp = p.ux;
  if (rand_value < D_rot_2_) {
    if (rand_value < D_rot_) {
      p.ux = p.uy;
      p.uy = -tmp;
    } else {
      p.ux = -p.uy;
      p.uy = tmp;
    }
  }
}
/**
* @brief Rotation of on-lattice particles
*
* @tparam Tx         Template type of x, y
* @param p           One particle
* @param rand_value  Random number range from 0 to 1
*/
template <typename Tx>
void UniDomain_AB_2::rot_diffuse(Par_s_2<Tx>& p, double rand_value) const {
  if (rand_value < D_rot_2_) {
    if (rand_value < D_rot_) {
      p.s -= 1;
      if (p.s < 0)
        p.s = 3;
    } else {
      p.s += 1;
      if (p.s > 3)
        p.s = 0;
    }
  }
}

  /**
 * @brief Translational motion of on-lattive active particles under PBC
 *
 * @tparam Tx             Template type of x, y
 * @tparam Tu             Template type of ux, uy
 * @param p               One particle
 * @param rand_value      Random number range from 0 to 1
 */
template <typename Tx, typename Tu>
void UniDomain_AB_2::trans_diffuse(Par_2<Tx, Tu>& p, double rand_value) {
  int ux, uy;
  if (rand_value < prob_slice_[0]) {
    ux = p.ux;
    uy = p.uy;
  } else if (rand_value < prob_slice_[1]) {
    ux = -p.ux;
    uy = -p.uy;
  }else if (rand_value < prob_slice_[2]) {
    ux = -p.uy;
    uy = p.ux;
  }else {
    ux = p.uy;
    uy = -p.ux;
  }
  int x_new = p.x + ux;
  tangle_x(x_new);
  int y_new = p.y + uy;
  tangle_y(y_new);
  const auto ic_new = get_ic(x_new, y_new);
  if (cell_[ic_new] < max_capacity_) {
    const auto ic = get_ic(p);
    cell_[ic] -= 1;
    p.x = x_new;
    p.y = y_new;
    cell_[ic_new] += 1;
  }
}

/**
* @brief Translational motion of on-lattive active particles under PBC
*
* @tparam Tx             Template type of x, y
* @param p               One particle
* @param rand_value      Random number range from 0 to 1
*/
template <typename Tx>
void UniDomain_AB_2::trans_diffuse(Par_s_2<Tx>& p, double rand_value) {
  int s;
  if (rand_value < prob_slice_[0]) {
    s = p.s;
  } else if (rand_value < prob_slice_[1]) {
    s = p.s + 2;
    if (s > 3)
      s -= 4;
  } else if (rand_value < prob_slice_[2]) {
    s = p.s + 1;
    if (s == 4)
      s = 0;
  } else {
    s = p.s - 1;
    if (s < 0)
      s = 3;
  }
  int x_new = p.x + state_2[s][0];
  tangle_x(x_new);
  int y_new = p.y + state_2[s][1];
  tangle_y(y_new);
  const auto ic_new = get_ic(x_new, y_new);
  if (cell_[ic_new] < max_capacity_) {
    const auto ic = get_ic(p);
    cell_[ic] -= 1;
    p.x = x_new;
    p.y = y_new;
    cell_[ic_new] += 1;
  }
}

  /**
 * @brief Translational motion of on-lattive active particles with walls
 * at xmin and xmax
 *
 * @tparam Tx            Template type of x, y
 * @tparam Tu            Template type of ux, uy
 * @param p              One particle
 * @param rand_value     Random number range from 0 to 1
 */
template <typename Tx, typename Tu>
void UniDomain_AB_2::trans_diffuse_wx(Par_2<Tx, Tu>& p, double rand_value) {
  int ux, uy;
  if (rand_value < prob_slice_[0]) {
    ux = p.ux;
    uy = p.uy;
  } else if (rand_value < prob_slice_[1]) {
    ux = -p.ux;
    uy = -p.uy;
  } else if (rand_value < prob_slice_[2]) {
    ux = -p.uy;
    uy = p.ux;
  } else {
    ux = p.uy;
    uy = -p.ux;
  }
  const int x_new = p.x + ux;
  if (x_new >= 0 && x_new < l_.x) {
    int y_new = p.y + uy;
    tangle_y(y_new);
    const auto ic_new = get_ic(x_new, y_new);
    if (cell_[ic_new] < max_capacity_) {
      const auto ic = get_ic(p);
      --cell_[ic];
      p.x = x_new;
      p.y = y_new;
      ++cell_[ic_new];
    }
  }
}

/**
* @brief Translational motion of on-lattive active particles with walls
* at xmin and xmax
*
* @tparam Tx            Template type of x, y
* @param p              One particle
* @param rand_value     Random number range from 0 to 1
*/
template <typename Tx>
void UniDomain_AB_2::trans_diffuse_wx(Par_s_2<Tx>& p, double rand_value) {
  int s;
  if (rand_value < prob_slice_[0]) {
    s = p.s;
  } else if (rand_value < prob_slice_[1]) {
    s = p.s + 2;
    if (s > 3)
      s -= 4;
  } else if (rand_value < prob_slice_[2]) {
    s = p.s + 1;
    if (s == 4)
      s = 0;
  } else {
    s = p.s - 1;
    if (s < 0)
      s = 3;
  }
  const int x_new = p.x + state_2[s][0];
  if (x_new >= 0 && x_new < l_.x) {
    int y_new = p.y + state_2[s][1];
    tangle_y(y_new);
    const auto ic_new = get_ic(x_new, y_new);
    if (cell_[ic_new] < max_capacity_) {
      const auto ic = get_ic(p);
      --cell_[ic];
      p.x = x_new;
      p.y = y_new;
      ++cell_[ic_new];
    }
  }
}

  /**
 * @brief Run simulation with periodic boundary condition
 * 
 * @tparam TPar       Template of particle
 * @tparam TRan       Template of random number generator
 * @param p_arr       Array of particles
 * @param myran       Random number generator
 * @param n_step      Total steps to run
 * @param seq_mode    Updating sequences for each time step
 */
template <typename TPar, typename TRan>
void UniDomain_AB_2::run(std::vector<TPar>& p_arr, TRan& myran,
                         int n_step, int seq_mode) {
  auto trival_move = [this](TPar &p, TRan &ran) {
    const auto r1 = ran.doub() * 2;
    const auto r2 = ran.doub();
    if (r1 > 1) {
      rot_diffuse(p, r1 * 0.5);
      trans_diffuse(p, r2);
    } else {
      trans_diffuse(p, r2);
      rot_diffuse(p, r1);
    }
    //trans_diffuse(p, ran.doub());
    //rot_diffuse(p, ran.doub());
  };
  base_run(p_arr, myran, trival_move, n_step, seq_mode);
}

/**
 * @brief Run simulation with two walls at xmin and xmax
 * 
 * @tparam TPar       Template of particle
 * @tparam TRan       Template of random number generator
 * @param p_arr       Array of particles
 * @param myran       Random number generator
 * @param n_step      Total steps to run
 * @param seq_mode    Updating sequences for each time step
 */
template <typename TPar, typename TRan>
void UniDomain_AB_2::run_wx(std::vector<TPar>& p_arr, TRan& myran,
                            int n_step, int seq_mode) {
  auto trival_move = [this](TPar &p, TRan &ran) {
    const auto r1 = ran.doub() * 2;
    const auto r2 = ran.doub();
    if (r1 > 1) {
      rot_diffuse(p, r1 * 0.5);
      trans_diffuse_wx(p, r2);
    } else {
      trans_diffuse_wx(p, r2);
      rot_diffuse(p, r1);
    }
    //trans_diffuse(p, ran.doub());
    //rot_diffuse(p, ran.doub());
  };
  base_run(p_arr, myran, trival_move, n_step, seq_mode);
}

} // end of namespace lattice.
#endif
