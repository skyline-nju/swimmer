#include <iostream>
#include "dynamic.h"
#include "cmdline.h"

int main(int argc, char* argv[]) {
  cmdline::parser cmd;
  cmd.add<double>("Lx", 'L', "System size in the x direction", true);
  cmd.add<double>("Ly", '\0', "System size in the y direction", false);
  cmd.add<double>("phi", '\0', "Packing fraction", false, 0);
  cmd.add<int>("nPar", 'N', "particle number", false, 0);
  cmd.add<double>("Pe", 'P', "the Peclet number", false, 0);
  cmd.add<double>("tau", 'T', "external torque", false, 0);
  cmd.add<double>("h", 'h', "time step to integrate", false, 1e-4);
  cmd.add<double>("eps", 'E', "epsilon_0", false, 1);
  cmd.add<double>("sigma", '\0', "particle size", false, 1);
  cmd.add<int>("nstep", '\0', "total time steps to run", false, 1000);
  cmd.add<int>("log_dt", '\0', "interval to record log", false, 100);
  cmd.add<int>("XY_dt", '\0', "interval to record xy information", false, 5000);
  cmd.add<unsigned long long>("seed", 's', "seed for random number", false, 1);
  cmd.parse_check(argc, argv);

  BaseDynamic_2 *simulator;
  //simulator = new BrownianDynamic<BP_2, CellList_w_list>(cmd);
  //simulator = new BrownianDynamic<Par_w_Pre_Pos<BP_2>, NeighborList_2>(cmd);
  simulator = new BrownianDynamic<BP_2, NeighborList_2>(cmd);

  simulator->run();
  delete simulator;

}