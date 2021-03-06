#include <iostream>
#include <chrono>
#include <cmdline.h>
#include "node.h"
#include "particle.h"
#include "uniDomain.h"
//#define USE_MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

int main(int argc, char* argv[]) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int my_rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

#ifndef USE_MPI
  cmdline::parser cmd;
  /* parameters */
  cmd.add<double>("Lx", 'L', "Domain size in the x direction", true);
  cmd.add<double>("Ly", '\0', "Domain size in the y direction", false);
  cmd.add<double>("phi", '\0', "packing fraction", true);
  cmd.add<int>("n_step", '\0', "total steps to run", true);
  cmd.add<double>("h", 'h', "time step for integration", true);
  cmd.add<double>("v0", 'v', "speed", false, 1);
  cmd.add<double>("alpha", '\0', "tumbling rate", true);
  cmd.add<unsigned long long>("seed", '\0', "random number seed", false, 1);
  cmd.add<double>("sigma", '\0', "particle size", false, 1);
  cmd.add<double>("spring_const", 'k', "spring const", false, 100);
  cmd.add<int>("int_mode", '\0', "integration mode", false, 0);
  cmd.add<double>("k_wall", '\0', "hardness of the wall", false, 100);

  /* boundary condition */
  cmd.add<int>("bc_x", '\0', "boundary condition in the x direction", false, 0);
  cmd.add<int>("bc_y", '\0', "boundary condition in the y direction", false, 0);

  /* output settings */
  cmd.add<string>("output", 'o', "path for outputting", false);
  cmd.add<int>("log_dt", '\0', "interval to update log", false, 0);
  cmd.add<int>("snap_dt", '\0', "interval to save snapshot", false, 0);
  cmd.add<int>("profile_dt", '\0', "interval to calculate wetting profile", false, 0);

  /* setting for calculating wetting profile */
  cmd.add<double>("eps", '\0', "min spacing between two neighbors", false, 1.1);
  cmd.add<int>("min_pts", '\0', "min number of neighbors of a core point", false, 2);
  cmd.add<double>("height_min", '\0', "min distance for a wetting particle", false, 0.55);
  cmd.parse_check(argc, argv);

  std::vector<BiNode<ActiveBrownPar_2>> p_arr;
  CellListNode_2<BiNode<ActiveBrownPar_2>> *cl;
  UniDomain_2 domain(cmd);
  domain.ini_rand(p_arr, &cl);
  domain.run(cmd, p_arr, cl);
  delete cl;
#else
  MPI_Finalize();
#endif
}