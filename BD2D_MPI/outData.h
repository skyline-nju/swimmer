#ifndef OUTDATA_H
#define OUTDATA_H
#include <fstream>
#include <chrono>
#include <ctime>
#include "comn.h"
#include "cmdline.h"

class BaseWriter {
public:
  typedef unsigned long long ull;

  // ReSharper disable once CppNonExplicitConvertingConstructor
  BaseWriter(const cmdline::parser &cmd);

  void set_para(const cmdline::parser &cmd);

protected:
  int idx_frame_;
  std::vector<int> frames_;

  double lx_;
  double ly_;
  int n_step_;
  double packing_frac_;
  int n_par_;
  ull seed_;
  double h_;
  double v0_;
  double tumbling_rate_;
  char basename_[100];
  std::ofstream *ptr_fout_;
};

class LogWriter : public BaseWriter {
public:
  LogWriter(const cmdline::parser &cmd, std::ofstream &fout);
  void set_frames(const cmdline::parser &cmd);
  void record(int i);

protected:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
};


class XyWriter: public BaseWriter {
public:
  XyWriter(const cmdline::parser &cmd, std::ofstream &fout);

  void set_frames(const cmdline::parser &cmd);

  template<typename TPar>
  void write(int i_step, std::vector<TPar> &p_arr);
  
};

template<typename TPar>
void XyWriter::write(int i_step, std::vector<TPar>& p_arr) {
  if (i_step == 0 || (!frames_.empty() && i_step == frames_[idx_frame_])) {
    if (i_step > 0)
      idx_frame_++;
    *ptr_fout_ << p_arr.size() << "\n";
    // comment line
    *ptr_fout_ << "Lattice=\"" << lx_ << " 0 0 0 " << ly_ << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step * h_;
    auto end = p_arr.cend();
    for (auto it = p_arr.cbegin(); it != end; ++it) {
      *ptr_fout_ << "\n" << std::fixed << std::setprecision(3) << "N\t"
        << (*it).x << "\t" << (*it).y;
    }
    *ptr_fout_ << std::endl;
  }
}
#endif
