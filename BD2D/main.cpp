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
  cmd.parse_check(argc, argv);

  DynamicBase_2 *bd;
  //bd = new BD_2<Vec_2<double>, NeighborList_2>(cmd);
  //bd = new ABD_2<OrientPar_2, NeighborList_2>(cmd);

  //bd = new BD_dipole_2<OrientPar_2, NeighborList_2>(cmd);
  bd = new CW_CCW_AB_2<OrientPar_2, NeighborList_2>(cmd);
  bd->run();
  delete[] bd;


  //ExtDipoleForce fed(80, 2, -0.5, 10, 4, 3.0 / 16);
  //Vec_2<double> x1(0, 0);
  //Vec_2<double> x2(3, 0);
  //Vec_2<double> dR(x1 - x2);
  //Vec_3<double> f1;
  //Vec_3<double> f2;
  //double theta1 = PI/2;
  //double theta2 = PI/2;
  //fed.eval(f1, f2, dR, theta1, theta2);
  //std::cout << f1.x << "\t" << f1.y << "\t" << f1.z << std::endl;
  //std::cout << f2.x << "\t" << f2.y << "\t" << f2.z << std::endl;

  //Vec_3<double> f3;
  //Vec_3<double> f4;
  //fed.eval2(f3, f4, dR, theta1, theta2);
  //std::cout << f3.x << "\t" << f3.y << "\t" << f3.z << std::endl;
  //std::cout << f4.x << "\t" << f4.y << "\t" << f4.z << std::endl;

}