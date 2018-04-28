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

namespace lattice {

class UniDomain_2 {
public:
  template<typename TPar, typename TRan>
  explicit UniDomain_2(const cmdline::parser &cmd, std::vector<TPar> &p_arr, TRan &myran);

  ~UniDomain_2() {delete log_; delete snap_; delete pf_;}

  template <typename TPar, typename TRan>
  void run(const cmdline::parser &cmd, std::vector<TPar>& p_arr,
           TRan &myran, int updating_order);
  /**
  * \brief Run and tumble model 1
  *
  * At each discrete time step, the lattice particle advances one step along
  * it's direction if not blocked by the wall at x=0, l_.x and the particles
  * in the target cell is less than max_capacity.
  *
  * \tparam TPar         Template for the particle
  * \tparam TRan         Template for random number generator
  * \param p             Particle
  * \param alpha         Tumbling rate
  * \param myran         Random number generator
  */
  template <typename TPar, typename TRan>
  void run_and_tumble_1(TPar &p, TRan &myran, double alpha);

  /**
   * @brief Translational motion of on-lattive active particles 
   * 
   * @tparam TPar           Template of particles
   * @param p               One particle
   * @param p_sep           Separation of probabilities for each event occuring
   * @param rand_value      Random number range from 0 to 1
   */
  template <typename TPar>
  void trans_diffuse(TPar &p, const double *p_sep, double rand_value);

  /**
   * @brief Rotation of on-lattice particles
   * 
   * @tparam TPar       Template of particles
   * @param p           One particle
   * @param D_rot       Rotational diffusion rate
   * @param rand_value  Random number range from 0 to 1
   */
  template <typename TPar>
  void rot_diffuse(TPar &p, double D_rot, double rand_value) const;

private:
  template <typename TPar>
  size_t get_ic(const TPar &p) const { return p.x + p.y * l_.x; }

  size_t get_ic(int ix, int iy) const { return ix + iy * l_.x; }

  template <typename T>
  void tangle_x(T &x) const { tangle_1(x, 0, l_.x, l_.x); }

  template <typename T>
  void tangle_y(T &y) const { tangle_1(y, 0, l_.y, l_.y); }

  template <typename TPar>
  void jump_x(TPar &p, int x_new);

  template <typename TPar>
  void jump_y(TPar &p, int y_new);


  Vec_2<int> l_;
  std::vector<uint8_t> cell_;
  const int max_capacity_;

  LogExporter *log_;
  SnapExporter_2 * snap_;
  ProfileExporter *pf_;
};

template <typename TPar, typename TRan>
UniDomain_2::UniDomain_2(const cmdline::parser& cmd, std::vector<TPar>& p_arr, TRan &myran)
  : max_capacity_(cmd.get<int>("n_max")), log_(nullptr), snap_(nullptr), pf_(nullptr) {
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
  set_output_2(cmd, &log_, &snap_, &pf_);
}

template <typename TPar, typename TRan>
void UniDomain_2::run(const cmdline::parser& cmd, std::vector<TPar>& p_arr,
                      TRan &myran, int updating_order) {
  const auto n_step = cmd.get<unsigned int>("n_step");
  const auto alpha = cmd.get<double>("alpha");
  const auto D_rot = cmd.get<double>("D_rot");
  double p_sep[3];

  std::function<void(TPar &, TRan &)> update_one;
  if (cmd.exist("alpha")) {
    std::cout << "run and tumble\n";
    update_one = [this, alpha](TPar &p, TRan &ran) {
      run_and_tumble_1(p, ran, alpha);
    };
  } else if (cmd.exist("nu_f")) {
    std::cout << "active brownain motion\n";
    const auto nu_f = cmd.get<double>("nu_f");
    const auto nu_b = cmd.get<double>("nu_b");
    const auto nu_t = cmd.get<double>("nu_t");

    const auto rate_tot = nu_f + nu_b + 2 * nu_t;
    p_sep[0] = nu_f / rate_tot;
    p_sep[1] = p_sep[0] + nu_b / rate_tot;
    p_sep[2] = p_sep[1] + nu_t / rate_tot;
    update_one = [this, p_sep, D_rot](TPar &p, TRan &ran) {
      const auto rand_value = ran.doub() * 2;
      const auto rand_value2 = ran.doub();
      if (rand_value > 1) {
        rot_diffuse(p, D_rot, rand_value * 0.5);
        trans_diffuse(p, p_sep, rand_value2);
      } else {
        trans_diffuse(p, p_sep, rand_value2);
        rot_diffuse(p, D_rot, rand_value);
      }
      //trans_diffuse(p, ran, p_sep);
      //rot_diffuse(p, D_rot, ran.doub());
    };
  }

  std::function<void()> update_all;
  std::vector<TPar*> ptr_arr;
  if (updating_order == 0) {
    std::cout << "random sampling\n";
    update_all = [update_one, &myran, &p_arr]() {
      const unsigned int n = p_arr.size();
      for (unsigned int count = 0; count < n; count++) {
        auto i = int(myran.doub() * n);
        update_one(p_arr[i], myran);
      }
    };
  } else if (updating_order == 1) {
    std::cout << "random sampling 1\n";
    update_all = [update_one, &myran, &p_arr]() {
      const unsigned int n = p_arr.size();
      auto idx = new int[n];
      for (unsigned int i = 0; i < n; i++)
        idx[i] = int(myran.doub() * n);
      for (unsigned int i = 0; i < n; i++)
        update_one(p_arr[idx[i]], myran);
      delete[] idx;
    };
  } else if (updating_order == 2) {
    std::cout << "shuffle the particles\n";
    update_all = [update_one, &myran, &p_arr]() {
      shuffle(p_arr, myran);
      for (auto &i : p_arr) {
        update_one(i, myran);
      }
    };
  } else if (updating_order == 3) {
    std::cout << "shuffle the pointers\n";
    ptr_arr.reserve(p_arr.size());
    for (auto &i : p_arr) {
      ptr_arr.push_back(&i);
    }
    update_all = [update_one, &myran,  &ptr_arr]() {
      shuffle(ptr_arr, myran);
      for (auto &i : ptr_arr) {
        update_one(*i, myran);
      }
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


template<typename TPar>
void UniDomain_2::jump_x(TPar & p, int x_new) {
  const auto ic_new = get_ic(x_new, p.y);
  if (cell_[ic_new] < max_capacity_) {
    --cell_[get_ic(p)];
    p.x = x_new;
    ++cell_[ic_new];
  }
}

template <typename TPar>
void UniDomain_2::jump_y(TPar& p, int y_new) {
  const auto ic_new = get_ic(p.x, y_new);
  if (cell_[ic_new] < max_capacity_) {
    --cell_[get_ic(p)];
    p.y = y_new;
    ++cell_[ic_new];
  }
}

/**
* \brief Run and tumble model 1
*
* At each discrete time step, the lattice particle advances one step along
* it's direction if not blocked by the wall at x=0, l_.x and the particles
* in the target cell is less than max_capacity.
* 
* \tparam TPar         Template for the particle
* \tparam TRan         Template for random number generator
* \param p             Particle
* \param alpha         Tumbling rate
* \param myran         Random number generator
*/
template <typename TPar, typename TRan>
void UniDomain_2::run_and_tumble_1(TPar& p, TRan &myran, double alpha) {
  if (p.uy) {
    int y_new = p.y + p.uy; //!< y_new should be a signed int
    tangle_y(y_new);
    jump_y(p, y_new);
  } else if ((p.ux < 0 && p.x > 0) || (p.ux > 0 && p.x < l_.x - 1)) {
    const auto x_new = p.x + p.ux;
    jump_x(p, x_new);
  }
  if (myran.doub() < alpha)
    p.tumble(myran);
}

  /**
   * @brief Translational motion of on-lattive active particles 
   * 
   * @tparam TPar           Template of particles
   * @param p               One particle
   * @param p_sep           Separation of probabilities for each event occuring
   * @param rand_value      Random number range from 0 to 1
   */
template <typename TPar>
void UniDomain_2::trans_diffuse(TPar& p, const double *p_sep, double rand_value) {
  if (rand_value < p_sep[0]) { //! jump forward
    if (p.uy) {
      int y_new = p.y + p.uy;
      tangle_y(y_new);
      jump_y(p, y_new);
    } else {
      int x_new = p.x + p.ux;
      tangle_x(x_new);
      jump_x(p, x_new);
    }
  } else if (rand_value < p_sep[1]) { //! jump backward
    if (p.uy) {
      int y_new = p.y - p.uy;
      tangle_y(y_new);
      jump_y(p, y_new);
    } else {
      int x_new = p.x - p.ux;
      tangle_x(x_new);
      jump_x(p, x_new);
    }
  } else if (rand_value < p_sep[2]) { //! jump left, (ux, uy) --> (-uy, ux)
    if (p.uy) {
      int x_new = p.x - p.uy;
      tangle_x(x_new);
      jump_x(p, x_new);
    } else {
      int y_new = p.y + p.ux;
      tangle_y(y_new);
      jump_y(p, y_new);
    }
  } else { //! jump right, (ux, uy) --> (uy, -ux)
    if (p.uy) {
      int x_new = p.x + p.uy;
      tangle_x(x_new);
      jump_x(p, x_new);
    } else {
      int y_new = p.y - p.ux;
      tangle_y(y_new);
      jump_y(p, y_new);
    }
  }
}
/**
 * @brief Rotation of on-lattice particles
 * 
 * @tparam TPar       Template of particles
 * @param p           One particle
 * @param D_rot       Rotational diffusion rate
 * @param rand_value  Random number range from 0 to 1
 */
template <typename TPar>
void UniDomain_2::rot_diffuse(TPar& p, double D_rot, double rand_value) const {
  auto tmp = p.ux;
  if (rand_value < D_rot) {
    p.ux = -p.uy;
    p.uy = tmp;
  } else if (rand_value > 1 - D_rot) {
    p.ux = p.uy;
    p.uy = -tmp;
  }
}

} // end of namespace lattice.
#endif
