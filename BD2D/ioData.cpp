#include "ioData.h"

double BaseWriter::Lx;
double BaseWriter::Ly;
int BaseWriter::nstep;
double BaseWriter::phi;
int BaseWriter::nPar;
unsigned long long BaseWriter::seed;
double BaseWriter::h;
double BaseWriter::Pe;
double BaseWriter::tau;
bool BaseWriter::is_active;
bool BaseWriter::is_chiral;

void BaseWriter::set_para(const cmdline::parser & cmd) {
  Lx = cmd.get<double>("Lx");
  Ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : Lx;
  nstep = cmd.get<int>("nstep");
  double sigma = cmd.get<double>("sigma");
  if (cmd.exist("nPar")) {
    nPar = cmd.get<int>("nPar");
    phi = cal_packing_fraction_2(nPar, Lx, Ly, sigma);
  } else if (cmd.exist("phi")) {
    phi = cmd.get<double>("phi");
    nPar = cal_particle_number_2(phi, Lx, Ly, sigma);
  }
  seed = cmd.get<unsigned long long>("seed");
  h = cmd.get<double>("h");
  if (cmd.exist("Pe")) {
    is_active = true;
    Pe = cmd.get<double>("Pe");
  } else {
    is_active = false;
    Pe = 0;
  }
  if (cmd.exist("tau")) {
    is_chiral = true;
    tau = cmd.get<double>("tau");
  } else {
    is_chiral = false;
  }
}

LogWriter::LogWriter(const cmdline::parser & cmd, std::ofstream & fout):
                     BaseWriter(cmd) {
  mkdir("log");
  char filename[100];
  snprintf(filename, 100, "log%s%g_%g_%g_%g_%llu.dat",
    delimiter.c_str(), phi, Pe, tau, Lx, seed);
  fout.open(filename);
  set_frames(cmd);

  t_start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(t_start);

  fout << "Started simulation at " << std::ctime(&start_time) << "\n";
  fout << "-------- Parameters --------\n";
  fout << "Particle number: " << nPar << "\n";
  fout << "Packing fraction: " << phi << "\n";
  fout << "Pe: " << Pe << "\n";
  fout << "torque: " << tau << "\n";
  fout << "Lx: " << Lx << "\n";
  fout << "Ly: " << Ly << "\n";
  fout << "seed: " << seed << "\n";
  fout << "Total time steps: " << nstep << "\n";
  fout << "Integration time step: " << h << "\n";
  fout << "\n";
  fout << "-------- RUN --------\n";
  fout << "time step\telapsed time" << std::endl;

  ptr_fout = &fout;
}

void LogWriter::set_frames(const cmdline::parser & cmd) {
  int dn = cmd.get<int>("log_dt");
  int i = dn;
  while (i <= nstep) {
    frames.push_back(i);
    i += dn;
  }
}

void LogWriter::operator() (int i) {
  if (!frames.empty() && i == frames[idx_frame]) {
    idx_frame++;
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start;
    double dt = elapsed_seconds.count();
    int hour = int(dt / 3600);
    int min = int((dt - hour * 3600) / 60);
    int sec = dt - hour * 3600 - min * 60;
    *ptr_fout << i << "\t" << hour << ":" << min << ":" << sec << std::endl;
    if (i == nstep) {
      std::time_t end_time = std::chrono::system_clock::to_time_t(t_now);
      *ptr_fout << "Finished simulation at " << std::ctime(&end_time) << "\n";
    }
  }

}

XY_Writer::XY_Writer(const cmdline::parser &cmd, std::ofstream &fout):
    BaseWriter(cmd) {
  mkdir("XY");
  char filename[100];
  snprintf(filename, 100, "XY%s%g_%g_%g_%g_%llu.extxyz",
    delimiter.c_str(), phi, Pe, tau, Lx, seed);
  fout.open(filename);
  set_frames(cmd);
  ptr_fout = &fout;
}

void XY_Writer::set_frames(const cmdline::parser & cmd) {
  int dn = cmd.get<int>("XY_dt");
  int i = dn;
  while (i <= nstep) {
    frames.push_back(i);
    i += dn;
  }
}
