#include <cstdint>
#include "latticeDomain.h"
#include "latticeParticle.h"
#include "rand.h"

int main(int argc, char* argv[]) {
  cmdline::parser cmd;
  /* parameters */
  cmd.add<int>("Lx", 'L', "Domain size int the x direction", true);
  cmd.add<int>("Ly", '\0', "Domain size int the y dirction", false);
  cmd.add<double>("pack_frac", '\0', "packing fraction", true);
  cmd.add<unsigned int>("n_step", 'n', "total steps to run", true);
  cmd.add<unsigned long long>("seed", 's', "random number seed", false, 1);
  cmd.add<int>("n_max", '\0', "maximum occupancy per site", false, 1);

  cmd.add<std::string>("output", 'o', "folder to output", false, "");
  cmd.add<int>("log_dt", '\0', "Time interval to record log", false, 0);
  cmd.add<int>("snap_dt", '\0', "Time interval to output snapshort", false, 0);
  cmd.add<int>("profile_dt", '\0', "Time inteval to output wetting profile",
               false, 0);

  /* define motion type */
  //! run and tumble
  cmd.add<double>("alpha", '\0', "tumbling_rate", false, 0);

  //! active brownian particle
  cmd.add<double>("nu_f", '\0', "rate of moving forward", false, 0);
  cmd.add<double>("nu_b", '\0', "rate of moving backward", false, 1.);
  cmd.add<double>("nu_t", '\0', "rate of moving transversely", false, 1.);
  cmd.add<double>("D_rot", '\0', "rate of rotational diffusion", false, 0.1);
  cmd.parse_check(argc, argv);


  const int n_step = cmd.get<unsigned int>("n_step");

  {
    Ranq2 myran(cmd.get<unsigned long long>("seed"));
    std::vector<lattice::Par_2<uint16_t, int8_t>> p_arr;
    lattice::UniDomain_AB_2 domain(cmd, p_arr, myran);
    domain.run(p_arr, myran, n_step, 0);
  }

}
