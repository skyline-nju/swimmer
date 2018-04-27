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
#include "rand.h"
#include "comn.h"
#include "cmdline.h"

namespace lattice {

template <typename TPar, typename UniFunc>
void update_shuffle_obj(std::vector<TPar> &p_arr, Ran &myran, UniFunc f_update) {
  shuffle(p_arr, myran);
  for (auto &i: p_arr) {
    f_update(i);
  }
}

template<typename TPar, typename UniFunc>
void update_shuffle_ptr(std::vector<TPar*> &ptr_arr, Ran &myran, UniFunc f_update) {
  shuffle(ptr_arr, myran);
  for (auto &i: ptr_arr) {
    f_update(*i);
  }
}

class UniDomain_2 {
public:
  template<typename TPar>
  explicit UniDomain_2(const cmdline::parser &cmd, std::vector<TPar> &p_arr);

  ~UniDomain_2() {delete log_; delete snap_; delete pf_;}

  template <typename TPar>
  void run(const cmdline::parser &cmd, std::vector<TPar>& p_arr, int shuffle_mode=0);

  /**
  * \brief Run and tumble model 1
  *
  * At each discrete time step, the lattice particle advances one step along
  * it's direction if not blocked by the wall at x=0, l_.x and the particles
  * in the target cell is less than max_capacity.
  * \tparam TPar         Template for the particle
  * \param p             Particle
  * \param alpha         Tumbling rate
  * \param max_capacity  Max number of particles that one cell can hold.
  */
  template <typename TPar>
  void run_and_tumble_1(TPar &p, double alpha, int max_capacity);

private:
  template <typename TPar>
  size_t get_ic(const TPar &p) const { return p.x + p.y * l_.x; }

  size_t get_ic(int ix, int iy) const { return ix + iy * l_.x; }

  Ran myran_;
  Vec_2<int> l_;
  std::vector<uint8_t> cell_;

  LogExporter *log_;
  SnapExporter_2 * snap_;
  ProfileExporter *pf_;
};

template <typename TPar>
UniDomain_2::UniDomain_2(const cmdline::parser& cmd, std::vector<TPar>& p_arr)
  : myran_(cmd.get<unsigned long long>("seed")), log_(nullptr), snap_(nullptr), pf_(nullptr) {
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
    TPar p(myran_, l_);
    const auto ic = get_ic(p);
    if (cell_[ic] == 0) {
      p_arr.push_back(p);
      cell_[ic] += 1;
    }
  }
  set_output_2(cmd, &log_, &snap_, &pf_);
}

template <typename TPar>
void UniDomain_2::run(const cmdline::parser& cmd, std::vector<TPar>& p_arr,
                      int shuffle_mode) {
  const auto n_step = cmd.get<unsigned int>("n_step");
  const auto alpha = cmd.get<double>("alpha");
  const auto max_capacity = cmd.get<int>("n_max");

  auto update_one = [this, alpha, max_capacity](TPar &p) {
    run_and_tumble_1(p, alpha, max_capacity);
  };

  std::function<void()> update_all;
  std::vector<TPar*> ptr_arr;
  if (shuffle_mode == 0) {
    update_all = [this, update_one, &p_arr]() {
      update_shuffle_obj(p_arr, myran_, update_one);
    };
  } else {
    ptr_arr.reserve(p_arr.size());
    for (auto &i : p_arr) {
      ptr_arr.push_back(&i);
    }
    update_all = [this, update_one, &ptr_arr]() {
      update_shuffle_ptr(ptr_arr, myran_, update_one);
    };
  }

  const auto t1 = std::chrono::system_clock::now();
  for (unsigned i = 1; i <= n_step; i++) {
    update_all();
    output_2(i, p_arr, cell_, log_, snap_, pf_);
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}

/**
* \brief Run and tumble model 1
*
* At each discrete time step, the lattice particle advances one step along
* it's direction if not blocked by the wall at x=0, l_.x and the particles
* in the target cell is less than max_capacity.
* \tparam TPar         Template for the particle
* \param p             Particle
* \param alpha         Tumbling rate
* \param max_capacity  Max number of particles that one cell can hold.
*/
template <typename TPar>
void UniDomain_2::run_and_tumble_1(TPar& p, double alpha, int max_capacity) {
  if (p.uy) {
    int y_new = p.y + p.uy; //!< y_new should be a signed int 
    tangle_1(y_new, 0, l_.y, l_.y);
    const auto ic_new = get_ic(p.x, y_new);
    if (cell_[ic_new] < max_capacity) {
      --cell_[get_ic(p)];
      p.y = y_new;
      ++cell_[ic_new];
    }
  } else if ((p.ux < 0 && p.x > 0) || (p.ux > 0 && p.x < l_.x - 1)) {
    const auto x_new = p.x + p.ux;
    const auto ic_new = get_ic(x_new, p.y);
    if (cell_[ic_new] < max_capacity) {
      --cell_[get_ic(p)];
      p.x = x_new;
      ++cell_[ic_new];
    }
  }
  if (myran_.doub() < alpha)
    p.tumble(myran_);
}

} // end of namespace lattice.
#endif
