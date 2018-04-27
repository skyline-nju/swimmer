#include <cstdint>
#include "latticeDomain.h"
#include "latticeParticle.h"

int main(int argc, char* argv[]) {
  cmdline::parser cmd;
  /* parameters */
  cmd.add<int>("Lx", 'L', "Domain size int the x direction", true);
  cmd.add<int>("Ly", '\0', "Domain size int the y dirction", false);
  cmd.add<double>("pack_frac", '\0', "packing fraction", true);
  cmd.add<unsigned int>("n_step", 'n', "total steps to run", true);
  cmd.add<unsigned long long>("seed", 's', "random number seed", false, 1);
  cmd.add<double>("alpha", '\0', "tumbling_rate", false, 0.001);
  cmd.add<int>("n_max", '\0', "maximum occupancy per site", false, 1);

  cmd.add<std::string>("output", 'o', "folder to output", false, "");
  cmd.add<int>("log_dt", '\0', "Time interval to record log", false, 0);
  cmd.add<int>("snap_dt", '\0', "Time interval to output snapshort", false, 0);
  cmd.add<int>("profile_dt", '\0', "Time inteval to output wetting profile",
               false, 0);

  cmd.parse_check(argc, argv);

  std::vector<lattice::Par_2<uint16_t, int8_t>> p_arr;
  lattice::UniDomain_2 domain(cmd, p_arr);
  domain.run(cmd, p_arr, 0);

}
