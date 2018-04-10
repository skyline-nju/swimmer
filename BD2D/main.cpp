#include <iostream>
#include "cmdline.h"
#include "brownian.h"
#include "neighborList.h"
#include "mixture.h"

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
  cmd.add<double>("dipole_strength", '\0', "strength of dipole", false, 0);
  cmd.add<double>("dipole_ratio", '\0', "ratio of two charges", false, 1);
  cmd.add("no_output", '\0', "no output");
  cmd.add("spatial_sort", '\0', "use spartial sorting");
  cmd.add<double>("phi_A", '\0', "Packing fraction of A particles", false, 0);
  cmd.parse_check(argc, argv);

  DynamicBase_2 *bd;
  //bd = new BD_2<Vec_2<double>, NeighborList_2>(cmd);
  //bd = new ABD_2<Par_w_theta_2, NeighborList_2>(cmd);

  //bd = new BD_dipole_2<Par_w_theta_2, NeighborList_2>(cmd);
  //bd = new CW_CCW_AB_2<Par_w_theta_2, NeighborList_2>(cmd);
  bd = new CABD_dipole_2<Par_w_u_2, NeighborList_2>(cmd);
  //bd = new CW_CCW_AB_Dipole_2<Par_w_u_2, NeighborList_2>(cmd);
  //bd = new Chrial_passive_2<Par_w_theta_2, NeighborList_2>(cmd);
  //bd = new Chrial_dipole_passive_2<Par_w_u_2, Vec_2<double>, CellList_list_2>(cmd);
  bd->run();
  delete bd;

}