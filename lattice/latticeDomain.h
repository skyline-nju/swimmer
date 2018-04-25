#ifndef LATTICLEDOMAIN_H
#define LATTICLEDOMAIN_H

#include <cstdint>
#include <vector>
#include "vect.h"
#include "rand.h"
#include "comn.h"
#include "cmdline.h"


class LatUniDomain_2 {
public:
  template<typename TPar>
  explicit LatUniDomain_2(const cmdline::parser &cmd,
                          std::vector<TPar> &p_arr);
  template <typename TPar>
  size_t get_ic(const TPar &p) const { return p.x + p.y * l_.x; }

  template <typename TPar>
  void run(const cmdline::parser &cmd, std::vector<TPar>& p_arr);

private:
  Ran myran_;
  Vec_2<unsigned int> l_;
  std::vector<uint8_t> cell_;
  unsigned int n_par_;

};

template <typename TPar>
LatUniDomain_2::LatUniDomain_2(const cmdline::parser& cmd,
                               std::vector<TPar>& p_arr)
  : myran_(cmd.get<uint64_t>("seed")) {
  l_.x = cmd.get<unsigned int>("Lx");
  l_.y = cmd.exist("Ly") ? cmd.get<unsigned int>("Ly") : l_.x;
  const auto phi = cmd.get<double>("pack_frac");
  const size_t ncells = l_.x * l_.y;
  n_par_ = unsigned int(ncells * phi);
  cell_.reserve(ncells);
  p_arr.reserve(n_par_);
  for (size_t i = 0; i < ncells; i++) {
    cell_.push_back(0);
  }
  while(p_arr.size() < n_par_) {
    TPar p(myran_, l_);
    const auto ic = get_ic(p);
    if (cell_[ic] == 0) {
      p_arr.push_back(p);
      cell_[ic] = 1;
    }
  }
}

template <typename TPar>
void LatUniDomain_2::run(const cmdline::parser& cmd, std::vector<TPar>& p_arr) {
  const auto n_step = cmd.get<int>("n_step");
  const auto alpha = cmd.get<double>("alpha");


  auto run_and_tumble = [this, alpha](TPar &p) {
    if (p.uy) {
      int y = p.y + p.uy;
      tangle_1(y, 0, l_.y, l_.y);
    }
    p.tumble(myran_);
  };


}


#endif
