#ifndef LATTICLEDOMAIN_H
#define LATTICLEDOMAIN_H

#include <cstdint>
#include <vector>
#include "vect.h"
#include "rand.h"
#include "comn.h"
#include "cmdline.h"

namespace lattice {
class UniDomain_2 {
public:
  template<typename TPar>
  explicit UniDomain_2(const cmdline::parser &cmd,
                          std::vector<TPar> &p_arr);
  template <typename TPar>
  size_t get_ic(const TPar &p) const { return p.x + p.y * l_.x; }

  size_t get_ic(int ix, int iy) const { return ix + iy * l_.x; }

  template <typename TPar>
  void run(const cmdline::parser &cmd, std::vector<TPar>& p_arr);

  template <typename TPar>
  void run(const cmdline::parser &cmd, std::vector<TPar*> &ptr_arr);

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
  Ran myran_;
  Vec_2<int> l_;
  std::vector<uint8_t> cell_;
  unsigned int n_par_;

};

template <typename TPar>
UniDomain_2::UniDomain_2(const cmdline::parser& cmd,  // NOLINT
                               std::vector<TPar>& p_arr)
  : myran_(cmd.get<uint64_t>("seed")) {
  l_.x = cmd.get<int>("Lx");
  l_.y = cmd.exist("Ly") ? cmd.get<int>("Ly") : l_.x;
  const auto phi = cmd.get<double>("pack_frac");
  const size_t ncells = l_.x * l_.y;
  n_par_ = unsigned int(ncells * phi);
  cell_.reserve(ncells);
  p_arr.reserve(n_par_);
  for (size_t i = 0; i < ncells; i++) {
    cell_.push_back(0);
  }
  while (p_arr.size() < n_par_) {
    TPar p(myran_, l_);
    const auto ic = get_ic(p);
    if (cell_[ic] == 0) {
      p_arr.push_back(p);
      cell_[ic] += 1;
    }
  }
}

template <typename TPar>
void UniDomain_2::run(const cmdline::parser& cmd, std::vector<TPar>& p_arr) {
  const auto n_step = cmd.get<unsigned int>("n_step");
  const auto alpha = cmd.get<double>("alpha");
  const auto max_capacity = cmd.get<int>("n_max");

  auto lambda = [this, alpha, max_capacity](TPar &p) {
    run_and_tumble_1(p, alpha, max_capacity);
  };
  const auto t1 = std::chrono::system_clock::now();
  for (unsigned i = 1; i <= n_step; i++) {
    shuffle(p_arr, myran_);
    auto end = p_arr.end();
    for (auto it = p_arr.begin(); it != end; ++it) {
      lambda(*it);
    }
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;
}

template <typename TPar>
void UniDomain_2::run(const cmdline::parser& cmd, std::vector<TPar*>& ptr_arr) {
  const auto n_step = cmd.get<unsigned int>("n_step");
  const auto alpha = cmd.get<double>("alpha");
  const auto max_capacity = cmd.get<int>("n_max");

  auto lambda = [this, alpha, max_capacity](TPar &p) {
    run_and_tumble_1(p, alpha, max_capacity);
  };
  const auto t1 = std::chrono::system_clock::now();
  for (unsigned i = 1; i <= n_step; i++) {
    shuffle(ptr_arr, myran_);
    auto end = ptr_arr.end();
    for (auto it = ptr_arr.begin(); it != end; ++it) {
      lambda(**it);
    }
  }
  const auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = t2 - t1;
  std::cout << "elapsed time: " << elapsed_time.count() << std::endl;

}

template <typename TPar>
void UniDomain_2::run_and_tumble_1(TPar& p, double alpha, int max_capacity) {
  if (p.uy) {
    auto y_new = p.y + p.uy;
    tangle_1(y_new, 0, l_.y, l_.y);
    const auto ic_new = get_ic(p.x, y_new);
    if (cell_[ic_new] < max_capacity) {
      p.y = y_new;
      ++cell_[ic_new];
    }
  } else if ((p.ux < 0 && p.x > 0) || (p.ux > 0 && p.x < l_.x - 1)) {
    const auto x_new = p.x + p.ux;
    const auto ic_new = get_ic(x_new, p.y);
    if (cell_[ic_new] < max_capacity) {
      p.x = x_new;
      ++cell_[ic_new];
    }
  }
  if (myran_.doub() < alpha)
    p.tumble(myran_);
}

}


#endif
