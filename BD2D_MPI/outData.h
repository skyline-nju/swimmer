#ifndef OUTDATA_H
#define OUTDATA_H
#include <fstream>
#include <chrono>
#include <ctime>
#include "comn.h"
#include "cmdline.h"
#include "boundary.h"
#include "cellList2D.h"
#include "wetting.h"

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
  void write(int i_step, const std::vector<TPar> &p_arr);
  
  //template<typename TPar>
  //void write(int i_step, const std::vector<TPar> &p_arr,
  //           const Cell_list_2<TPar> &cl, const Wall_x_PBC_y_2 &bc);

  template<typename TPar>
  void write_cluster(int i_step, const std::vector<TPar> &p_arr,
                     const Cell_list_2<TPar> &cl, const Wall_x_PBC_y_2 &bc);
};

template<typename TPar>
void XyWriter::write(int i_step, const std::vector<TPar>& p_arr) {
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

//template <typename TPar>
//void XyWriter::write(int i_step, const std::vector<TPar>& p_arr,
//                     const Cell_list_2<TPar>& cl, const Wall_x_PBC_y_2& bc) {
//  if (i_step == 0 || (!frames_.empty() && i_step == frames_[idx_frame_])) {
//    std::vector<int> marker(p_arr.size());
//    mark_wetting_par(1.0, 1, p_arr, cl, bc, marker);
//    if (i_step > 0)
//      idx_frame_++;
//    *ptr_fout_ << p_arr.size() << "\n";
//    // comment line
//    *ptr_fout_ << "Lattice=\"" << lx_ << " 0 0 0 " << ly_ << " 0 0 0 1\" "
//      << "Properties=species:S:1:pos:R:2 "
//      << "Time=" << i_step * h_;
//    for (size_t i = 0; i < p_arr.size(); i++) {
//      *ptr_fout_ << "\n";
//      if (marker[i] == 0)
//        *ptr_fout_ << "B\t";
//      else
//        *ptr_fout_ << "A\t";
//      *ptr_fout_ << std::fixed << std::setprecision(3) << p_arr[i].x << "\t" << p_arr[i].y;
//    }
//    *ptr_fout_ << std::endl;
//  }
//  
//}

template <typename TPar>
void XyWriter::write_cluster(int i_step, const std::vector<TPar>& p_arr, const Cell_list_2<TPar>& cl,
  const Wall_x_PBC_y_2& bc) {
  if (i_step == 0 || (!frames_.empty() && i_step == frames_[idx_frame_])) {
    if (i_step > 0)
      idx_frame_++;
    std::vector<Cluster> c_arr;
    std::vector<bool> flag_clustered(p_arr.size(), false);
    dbscan(c_arr, flag_clustered, 1.1, 2, p_arr, bc);
    *ptr_fout_ << p_arr.size() << "\n";
        // comment line
    *ptr_fout_ << "Lattice=\"" << lx_ << " 0 0 0 " << ly_ << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2:mass:M:1 "
      << "Time=" << i_step * h_;
    for (size_t i = 0; i < p_arr.size(); i++) {
      *ptr_fout_ << "\nB\t" << std::fixed << std::setprecision(3)
                 << p_arr[i].x << "\t" << p_arr[i].y << "\t"
                 << flag_clustered[i];
    }
    *ptr_fout_ << std::endl;
  }
}
#endif
