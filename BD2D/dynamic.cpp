#include "dynamic.h"
#include <functional>

BD_2::BD_2(const cmdline::parser & cmd) {
  Lx = cmd.get<double>("Lx");
  Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Lx;
  std::cout << "Lx = " << Lx << "\tLy = " << Ly << "\n";

  myran = new Ran(cmd.get<unsigned long long>("seed"));
  nstep = cmd.get<int>("nstep");

  euler = new Int_EM_2(cmd.get<double>("h"), Lx, Ly);
  fwca = new F_WCA_2(Lx, Ly, cmd.get<double>("eps"));
  set_output(cmd);
  create_particles(cmd);
  std::cout << "Finish initilization." << std::endl;
}

BD_2::~BD_2() {
  delete[] par;
  delete euler;
  delete fwca;
  delete xy_out;
  delete log_out;
  delete myran;
}

void BD_2::set_output(const cmdline::parser &cmd) {
  BaseWriter::set_para(cmd);
  fout.emplace_back();
  fout.emplace_back();
  xy_out = new XY_Writer(cmd, fout[0]);
  log_out = new LogWriter(cmd, fout[1]);
}

void BD_2::create_particles(const cmdline::parser & cmd) {
  double phi;  // packing fraction
  double sigma = cmd.get<double>("sigma");
  if (cmd.exist("nPar")) {
    nPar = cmd.get<int>("nPar");
    phi = cal_packing_fraction_2(nPar, Lx, Ly, sigma);
  } else if (cmd.exist("phi")) {
    phi = cmd.get<double>("phi");
    nPar = cal_particle_number_2(phi, Lx, Ly, sigma);
  } else {
      std::cout << "Error, need specify nPar or phi" << std::endl;
      exit(1);
  }
  std::cout << "Number of Particles: " << nPar << std::endl;

  par = new BP_2[nPar];
  if (phi < 0.4) {
    create_rand_2(par, nPar, sigma, Lx, Ly, myran);
  } else {
     // .......
  }
  (*xy_out)(0, par);
}

void BD_2::run() {
  auto force = [=](auto &p1, auto &p2) {
    fwca->pair(p1, p2);
  };
  auto integ = [=](auto *p1) {
    euler->int_T(p1, myran);
  };
  for (int i = 1; i <= nstep; i++) {
    for_each_pair(par, nPar, force);
    for_each_particle(par, nPar, integ);
    (*xy_out)(i, par);
    (*log_out)(i);
  }

}
