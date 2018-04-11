#include "outData.h"

BaseWriter::BaseWriter(const cmdline::parser &cmd): idx_frame_(0) {  // NOLINT
  set_para(cmd);
  snprintf(basename_, 100, "%g_%g_%g_%g_%llu", lx_, ly_, packing_frac_,
    tumbling_rate_, seed_);
  ptr_fout_ = nullptr;
}

void BaseWriter::set_para(const cmdline::parser &cmd) {
  lx_ = cmd.get<double>("Lx");
  ly_ = cmd.exist("Ly") ? cmd.get<double>("Ly") : lx_;
  n_step_ = cmd.get<int>("n_step");
  packing_frac_ = cmd.get<double>("phi");
  n_par_ = cal_particle_number_2(packing_frac_, lx_, ly_, 1);
  seed_ = cmd.get<unsigned long long>("seed");
  h_ = cmd.get<double>("h");
  v0_ = cmd.get<double>("v0");
  tumbling_rate_ = cmd.get<double>("alpha");
}

LogWriter::LogWriter(const cmdline::parser& cmd, std::ofstream &fout):
    BaseWriter(cmd) {
  mkdir("log");
  char filename[100];
  snprintf(filename, 100, "log%s%s.dat", delimiter.c_str(), basename_);
  fout.open(filename);
  t_start_ = std::chrono::system_clock::now();
  auto start_time = std::chrono::system_clock::to_time_t(t_start_);

#ifdef _MSC_VER
  char str[26];
  ctime_s(str, sizeof(str), &start_time);
  fout << "Started simulation at " << str << "\n";
 #else
  fout << "Started simulation at " << ctime(&start_time) << "\n";
#endif
  fout << "-------- Parameters --------\n";
  fout << "Particle number: " << n_par_ << "\n";
  fout << "Packing fraction: " << packing_frac_ << "\n";
  fout << "Lx: " << lx_ << "\n";
  fout << "Ly: " << ly_ << "\n";
  fout << "tumbling_rate" << tumbling_rate_ << "\n";
  fout << "h: " << h_ << "\n";
  fout << "n_step: " << n_step_ << "\n";
  fout << "seed: " << seed_ << "\n";
  fout << "\n";
  fout << "-------- RUN --------\n";
  fout << "time step\telapsed time" << std::endl;

  ptr_fout_ = &fout;
}

void LogWriter::set_frames(const cmdline::parser& cmd) {
  int dn;
  if (cmd.exist("log_dt")) {
    dn = cmd.get<int>("log_dt");
  } else {
    dn = int(1 / h_) * 10;
  }
  int i = dn;
  while (i <= n_step_) {
    frames_.push_back(i);
    i += dn;
  }
}

void LogWriter::record(int i) {
  if (!frames_.empty() && i == frames_[idx_frame_]) {
    idx_frame_++;
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    double dt = elapsed_seconds.count();
    int hour = int(dt / 3600);
    int min = int((dt - hour * 3600) / 60);
    int sec = dt - hour * 3600 - min * 60;
    *ptr_fout_ << i << "\t" << hour << ":" << min << ":" << sec << std::endl;
    if (i == n_step_) {
      auto end_time = std::chrono::system_clock::to_time_t(t_now);
#ifdef _MSC_VER
      char str[26];
      ctime_s(str, 26, &end_time);
      *ptr_fout_ << "Finished simulation at " << str << "\n";
#else
      *ptr_fout_ << "Finished simulation at " << ctime(&end_time) << "\n";
#endif
    }
  }
}

XyWriter::XyWriter(const cmdline::parser & cmd, std::ofstream & fout):
    BaseWriter(cmd) {
  mkdir("XY");
  char filename[100];
  snprintf(filename, 100, "XY%s%s.extxyz", delimiter.c_str(), basename_);
  
  fout.open(filename);
  set_frames(cmd);
  ptr_fout_ = &fout;
}

void XyWriter::set_frames(const cmdline::parser & cmd) {
  int dn;
  if (cmd.exist("XY_dt")) {
    dn = cmd.get<int>("XY_dt");
  } else {
    dn = int(1 / h_ * 0.2);
  }
  int i = dn;
  while (i <= n_step_) {
    frames_.push_back(i);
    i += dn;
  }
}
