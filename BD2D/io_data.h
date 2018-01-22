#ifndef IO_DATA_H
#define IO_DATA_H
#include <fstream>
#include <ctime>
#include <chrono>
#include "comn.h"
#include "cmdline.h"

class BaseWriter {
public:
  BaseWriter(const cmdline::parser &cmd) : idx_frame(0) {}
  static void set_para(const cmdline::parser &cmd);
  virtual void set_frames(const cmdline::parser &cmd) = 0;

protected:
  std::vector<int> frames;
  int idx_frame;
  std::ofstream *ptr_fout;

  static double Lx;
  static double Ly;
  static int nstep;
  static double phi;
  static int nPar;
  static unsigned long long seed;
  static double h;
  static double Pe;
  static double tau;
  static bool is_active;
  static bool is_chiral;
};

class LogWriter: public BaseWriter {
public:
  LogWriter(const cmdline::parser &cmd, std::ofstream &fout);
  void set_frames(const cmdline::parser &cmd);
  void operator() (int idx_step);

private:
  std::chrono::time_point<std::chrono::system_clock> t_start;
};

class XY_Writer : public BaseWriter {
public:
  XY_Writer(const cmdline::parser &cmd, std::ofstream &fout);

  void set_frames(const cmdline::parser &cmd);
  
  template <class Par>
  void operator() (int i, Par *p);
};

template <class Par>
void XY_Writer::operator() (int i, Par *p) {
  if (i==0 || (!frames.empty() && i == frames[idx_frame])) {
    if (i > 0)
      idx_frame++;
    *ptr_fout << nPar << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
              << "Properties=species:S:1:pos:R:2 "
              << "Time=" << i * h;
    for (int j = 0; j < nPar; j++) {
      *ptr_fout << "\n" << std::setprecision(3) << "N\t" 
                << p[j].x << "\t" << p[j].y;
    }
    *ptr_fout << std::endl;
  }
}

#endif