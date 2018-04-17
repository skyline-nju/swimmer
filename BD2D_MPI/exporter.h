#ifndef EXPORT_H
#define EXPORT_H
#include <fstream>
#include <chrono>
#include <iomanip>

#include <cmdline.h>
#include "vect.h"

class BaseExporter_2 {
public:
  typedef unsigned long long ull;

  explicit BaseExporter_2(const cmdline::parser &cmd, int mode=0);

  void set_lin_frame(int sep);

  bool need_export(int i_step);

  int get_n_frames() const {return frames_arr_.size();}

protected:
  Vec_2<double> domain_length_;
  double pack_frac_;
  double sigma_;
  size_t n_par_;
  int n_step_;
  double h_;
  double v0_;
  ull seed_;
  double tumbling_rate_;
  std::string folder_;
  std::string path_;

private:
  int iframe_;
  std::vector<int> frames_arr_;
};


class LogExporter_2: public BaseExporter_2 {
public:
  LogExporter_2(const cmdline::parser &cmd);
  void record(int i_step);

private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  std::ofstream fout_;
};

class XyExporter: public BaseExporter_2 {
public:
  XyExporter(const cmdline::parser &cmd);
  template<typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr);

private:
  void write_lattice();
  std::ofstream fout_;
};

template <typename TPar>
void XyExporter::write_frame(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    fout_ << p_arr.size() << "\n";
    // comment line
    write_lattice();
    fout_ << "Properties=species:S:1:pos:R:2 "
            << "Time=" << i_step * h_;
    auto end = p_arr.cend();
    for (auto it = p_arr.cbegin(); it != end; ++it) {
      fout_ << "\n" << std::fixed << std::setprecision(3) << "N\t"
        << (*it).x << "\t" << (*it).y;
    }
    fout_ << std::endl;
  }
}
#endif
