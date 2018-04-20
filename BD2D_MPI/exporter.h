#ifndef EXPORT_H
#define EXPORT_H
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cmdline.h>
#include "vect.h"

/**
 * \brief Basic class for exporter
 */
class BaseExporter_2 {
public:
  typedef unsigned long long ull;

  explicit BaseExporter_2(const cmdline::parser &cmd, int mode=0);

  void set_lin_frame(int sep);

  bool need_export(int i_step);

  size_t get_n_frames() const {return frames_arr_.size();}

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
  double particle_hardness_;
  double wall_hardness_;
  std::string folder_;
  std::string path_;
  int frame_interval_;

  double eps_;  // two particle are treated as neighbors within eps.
  int min_pts_; // the minimum neighbors that one core point should have.
  double height_min_;  // the minimum distance from the "wetting" particle to the wall.
private:
  int iframe_;
  std::vector<int> frames_arr_;
};


/**
 * \brief Class for writting the log file
 */
class LogExporter_2: public BaseExporter_2 {
public:
  explicit LogExporter_2(const cmdline::parser &cmd);
  ~LogExporter_2();
  void record(int i_step);

private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  std::ofstream fout_;
};

/**
 * \brief Class for outputing the snapshot as XYZ format.
 */
class XyExporter: public BaseExporter_2 {
public:
  explicit XyExporter(const cmdline::parser &cmd);

  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr);

  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr,
                   const std::vector<char> &p_type_arr);
private:
  void write_head_lines(int n_par, int i_step);
  std::ofstream fout_;
};

template <typename TPar>
void XyExporter::write_frame(int i_step, const std::vector<TPar>& p_arr) {
  write_head_lines(p_arr.size(), i_step);
  auto end = p_arr.cend();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    fout_ << "\nV\t" << std::fixed << std::setprecision(3)
          << (*it).x << "\t" << (*it).y;
  }
  fout_ << std::endl;
}

template <typename TPar>
void XyExporter::write_frame(int i_step, const std::vector<TPar>& p_arr,
                             const std::vector<char>& p_type_arr) {
  write_head_lines(p_arr.size(), i_step);
  auto end = p_arr.cend();
  auto it_type = p_type_arr.cbegin();
  for (auto it = p_arr.cbegin(); it != end; ++it) {
    fout_ << "\n" << *it_type << "\t"
      << std::fixed << std::setprecision(3)
      << (*it).x << "\t" << (*it).y;
    ++it_type;
  }
  fout_ << std::endl;
}

#endif
