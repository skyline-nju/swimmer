#include <iostream>
#include <chrono>
#include <cmdline.h>
#include "node.h"
#include "particle.h"
#include "singleDomain2D.h"
#include "boundary.h"
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
  cmd.add<double>("Lx", 'L', "Domain size in the x direction", true);
  cmd.add<double>("Ly", '\0', "Domain size in the y direction", false);
  cmd.add<double>("phi", '\0', "packing fraction", true);
  cmd.add<int>("n_step", '\0', "total steps to run", true);
  cmd.add<double>("h", 'h', "time step for integration", true);
  cmd.add<double>("v0", 'v', "speed", false, 1);
  cmd.add<double>("alpha", '\0', "tumbling rate", true);
  cmd.add<unsigned long long>("seed", '\0', "random number seed", false, 1);
  cmd.add<double>("sigma", '\0', "particle size", false, 1);
  cmd.add<int>("log_dt", '\0', "interval to update log", false, 1000);
  cmd.add<int>("snap_dt", '\0', "interval to save snapshot", false, 1000);
  cmd.add<double>("spring_const", 'k', "spring const", false, 100);
  cmd.add<int>("int_mode", '\0', "integration mode", false, 0);
  cmd.add<double>("k_wall", '\0', "hardness of the wall", false, 100);
  cmd.add<string>("output", 'o', "path for outputting", false);
  cmd.parse_check(argc, argv);

  {
    Single_domain_2<BiNode<ActiveBrownPar_2>, Wall_x_PBC_y_2> s(cmd);
    s.ini_rand();
    s.run(cmd);
  }
#else
  MPI_Finalize();
#endif
}