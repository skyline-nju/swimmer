#include "latticeDomain.h"
#include "latticeParticle.h"

int main(int argc, char* argv[]) {
  cmdline::parser cmd;
  /* parameters */
  cmd.add<unsigned int>("Lx", 'L', "Domain size int the x direction", true);
  cmd.add<unsigned int>("Ly", '\0', "Domain size int the y dirction", false);
  cmd.add<double>("pack_frac", '\0', "packing fraction", true);
  cmd.add<unsigned int>("n_step", 'n', "total steps to run", true);
  cmd.add<unsigned long long>("seed", 's', "random number seed", false, 1);
  cmd.add<double>("alpha", '\0', "tumbling_rate", false, 0.001);
  cmd.add<unsigned int>("n_max", '\0', "maximum occupancy per site", false, 1);

  cmd.parse_check(argc, argv);

  int i = -1;
  if (i)
    std::cout << i << std::endl;
}
