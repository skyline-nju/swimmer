#include "dynamic.h"

BaseDynamic_2::BaseDynamic_2(const cmdline::parser & cmd) {
  Lx = cmd.get<double>("Lx");
  Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Lx;
  std::cout << "Lx = " << Lx << "\tLy = " << Ly << "\n";
  myran = new Ran(cmd.get<unsigned long long>("seed"));
  nstep = cmd.get<int>("nstep");
  double sigma = cmd.get<double>("sigma");
  if (cmd.exist("nPar")) {
    nPar = cmd.get<int>("nPar");
    packing_fraction = cal_packing_fraction_2(nPar, Lx, Ly, sigma);
  } else if (cmd.exist("phi")) {
    packing_fraction = cmd.get<double>("phi");
    nPar = cal_particle_number_2(packing_fraction, Lx, Ly, sigma);
  } else {
    std::cout << "Error, need specify nPar or phi" << std::endl;
    exit(1);
  }

  euler = new Int_EM_2(cmd.get<double>("h"), Lx, Ly);
  fwca = new F_WCA_2(Lx, Ly, cmd.get<double>("eps"));

  out_on = false;
  //out_on = true;
  if (out_on) {
    BaseWriter::set_para(cmd);
    fout.emplace_back();
    fout.emplace_back();
    xy_out = new XY_Writer(cmd, fout[0]);
    log_out = new LogWriter(cmd, fout[1]);
  }
}

BaseDynamic_2::~BaseDynamic_2() {
  delete euler;
  delete fwca;
  if (out_on) {
    delete xy_out;
    delete log_out;
  }
  delete myran;
}

void BaseDynamic_2::run() {
  auto t_beg = std::chrono::system_clock::now();
  run(nstep);
  auto t_end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_sec = t_end - t_beg;
  double dt = elapsed_sec.count();
  std::cout << nPar << "particles, " << nstep << " steps, " << dt << "s\n";
  double speed = nPar * nstep / dt;
  std::cout << std::scientific << "speed = " << speed 
            << " particle * step / s\n";
  std::cout << std::endl;
}
