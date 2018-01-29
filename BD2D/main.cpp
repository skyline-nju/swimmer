#include <iostream>
#include "cmdline.h"
#include "brownian.h"
#include "neighborList.h"

int main(int argc, char* argv[]) {
  cmdline::parser cmd;
  cmd.add<double>("Lx", 'L', "System size in the x direction", true);
  cmd.add<double>("Ly", '\0', "System size in the y direction", false);
  cmd.add<double>("phi", '\0', "Packing fraction", false, 0);
  cmd.add<int>("nPar", 'N', "particle number", false, 0);
  cmd.add<double>("Pe", 'P', "the Peclet number", false, 0);
  cmd.add<double>("tau", 'T', "external torque", false, 0);
  cmd.add<double>("h", 'h', "time step to integrate", false, 2e-5);
  cmd.add<double>("eps", 'E', "epsilon_0", false, 1);
  cmd.add<double>("sigma", '\0', "particle size", false, 1);
  cmd.add<int>("nstep", '\0', "total time steps to run", false, 1000);
  cmd.add<int>("log_dt", '\0', "interval to record log", false, 10000);
  cmd.add<int>("XY_dt", '\0', "interval to record xy information", false, 5000);
  cmd.add<unsigned long long>("seed", 's', "seed for random number", false, 1);
  cmd.add("no_output", '\0', "no output");
  cmd.add("spatial_sort", '\0', "use spartial sorting");
  cmd.parse_check(argc, argv);

  DynamicBase_2 *bd;
  //bd = new BD_2<Vec_2<double>, NeighborList_2>(cmd);
  bd = new ABD_2<OrientPar_2, NeighborList_2>(cmd);
  bd->run();
  delete[] bd;
}