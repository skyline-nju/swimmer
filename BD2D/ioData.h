#ifndef IODATA_H
#define IODATA_H
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
  void operator() (int i, std::vector<Par> &p);

  template <class Par>
  void write_double(int i, std::vector<Par> &p);

  template <class Par>
  void write_mix(int i_step, std::vector<Par> &p, int n_sep);

  template <class Par1, class Par2>
  void write_mix(int i_step, std::vector<Par1> &p1, std::vector<Par2> &p2);

  template <class Par>
  void write_theta(int i_step, std::vector<Par> &p);

  template <class Par>
  void write_theta(int i_step, std::vector<Par> &p, int n_sep);

};

template <class Par>
void XY_Writer::operator() (int i, std::vector<Par> &p) {
  if (i==0 || (!frames.empty() && i == frames[idx_frame])) {
    if (i > 0)
      idx_frame++;
    *ptr_fout << nPar << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
              << "Properties=species:S:1:pos:R:2 "
              << "Time=" << i * h;
    for (int j = 0; j < nPar; j++) {
      *ptr_fout << "\n" << std::fixed << std::setprecision(3) << "N\t" 
                << p[j].x << "\t" << p[j].y;
    }
    *ptr_fout << std::endl;
  }
}

template <class Par>
void XY_Writer::write_double(int i, std::vector<Par> &p) {
  if (i == 0 || (!frames.empty() && i == frames[idx_frame])) {
    if (i > 0)
      idx_frame++;
    *ptr_fout << nPar * 2 << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i * h;
    for (int j = 0; j < nPar; j++) {
      *ptr_fout << "\n" << std::fixed << std::setprecision(6) << "N\t"
        << p[j].x << "\t" << p[j].y;
      //double dx = 0.01 * std::cos(p[j].theta);
      //double dy = 0.01 * std::sin(p[j].theta);
      double dx = 0.01 * p[j].u.x;
      double dy = 0.01 * p[j].u.y;
      *ptr_fout << "\n" << std::fixed << std::setprecision(6) << "O\t"
        << p[j].x + dx << "\t" << p[j].y + dy;
    }
    *ptr_fout << std::endl;
  }
}

template<class Par>
void XY_Writer::write_mix(int i_step, std::vector<Par>& p, int n_sep) {
  if (i_step == 0 || (!frames.empty() && i_step == frames[idx_frame])) {
    if (i_step > 0)
      idx_frame++;
    *ptr_fout << nPar << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step * h;
    for (int j = 0; j < n_sep; j++) {
      *ptr_fout << "\n" << std::fixed << std::setprecision(6) << "N\t"
        << p[j].x << "\t" << p[j].y;
    }
    for (int j = n_sep; j < nPar; j++) {
      *ptr_fout << "\n" << std::fixed << std::setprecision(6) << "O\t"
        << p[j].x << "\t" << p[j].y;
    }
    *ptr_fout << std::endl;
  }
}

template <class Par1, class Par2>
void XY_Writer::write_mix(int i_step, std::vector<Par1> &p1, std::vector<Par2> &p2) {
  if (i_step == 0 || (!frames.empty() && i_step == frames[idx_frame])) {
    if (i_step > 0)
      idx_frame++;
    *ptr_fout << nPar << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2:mass:M:1 "
      << "Time=" << i_step * h;
    int n1 = p1.size();
    for (int j = 0; j < n1; j++) {
      double theta = std::atan2(p1[j].u.y, p1[j].u.x) / PI * 180;
      if (theta < 0)
        theta += 360;
      *ptr_fout << "\n" << std::fixed << std::setprecision(6) << "N\t"
        << p1[j].x << "\t" << p1[j].y << "\t" << theta;
    }
    int n2 = p2.size();
    for (int j = 0; j < n2; j++) {
      *ptr_fout << "\n" << std::fixed << std::setprecision(6) << "O\t"
        << p2[j].x << "\t" << p2[j].y << "\t0";
    }
    *ptr_fout << std::endl;
  }
}

template<class Par>
void XY_Writer::write_theta(int i, std::vector<Par>& p) {
  if (i == 0 || (!frames.empty() && i == frames[idx_frame])) {
    if (i > 0)
      idx_frame++;
    *ptr_fout << nPar << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2:mass:M:1 "
      << "Time=" << i * h;
    for (int j = 0; j < nPar; j++) {
      double theta = std::atan2(p[j].u.y, p[j].u.x) / PI * 180;
      if (theta < 0)
        theta += 360;
      *ptr_fout << "\n" << std::fixed << std::setprecision(3) << "N\t"
        << p[j].x << "\t" << p[j].y << "\t" << theta;
    }
    *ptr_fout << std::endl;
  }
}

template<class Par>
void XY_Writer::write_theta(int i, std::vector<Par>& p, int n_sep) {
  if (i == 0 || (!frames.empty() && i == frames[idx_frame])) {
    if (i > 0)
      idx_frame++;
    *ptr_fout << nPar << "\n";
    // comment line
    *ptr_fout << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2:mass:M:1 "
      << "Time=" << i * h;
    for (int j = 0; j < nPar; j++) {
      double theta = std::atan2(p[j].u.y, p[j].u.x) / PI * 180;
      if (theta < 0)
        theta += 360;
      if (j < n_sep) {
        *ptr_fout << "\nN\t";
      } else {
        *ptr_fout << "\nO\t";
      }
      *ptr_fout << std::fixed << std::setprecision(3)
        << p[j].x << "\t" << p[j].y << "\t" << theta;
    }
    *ptr_fout << std::endl;
  }
}


#endif