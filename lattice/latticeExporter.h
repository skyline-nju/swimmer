#ifndef LATTICEEXPORTER_H
#define LATTICEEXPORTER_H
#include <fstream>
#include <chrono>
#include <iomanip>
#include "cmdline.h"
#include "comn.h"
#include "../BD2D_MPI/exporter.h"

namespace lattice {

/*************************************************************************//**
 * \brief Log exporter
 ****************************************************************************/
class LogExporter:public BaseLogExporter {
public:
  explicit LogExporter(const cmdline::parser &cmd);
};

/**************************************************************************//**
* \brief Snapshot exporter in netCDF format
******************************************************************************/
class SnapExporter_2:public BaseExporter {
public:
  explicit SnapExporter_2(const cmdline::parser &cmd);
  ~SnapExporter_2();

  void open(const cmdline::parser &cmd);
  void set_chunk_and_deflate() const;
  void put_time_step(int i_step) const;
  void put_data(const int *coor, const int8_t *ori) const;
  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr);

protected:
  /* id for each variables */
  int ncid_;
  int spatial_id_;
  int cell_spatial_id_;
  int time_id_;
  int cell_lengths_id_;
  int coordinates_id_;
  int orientations_id_;

  /* dimension lengths */
  const size_t frame_len_;
  const size_t spatial_len_;
  const size_t atom_len_;
  const size_t cell_spatial_len_;

private:
  int deflate_level_;
  size_t time_idx_[1];
};

template <typename TPar>
void SnapExporter_2::write_frame(int i_step, const std::vector<TPar>& p_arr) {
  put_time_step(i_step);

  const size_t n = p_arr.size();
  int *coor = new int[n * 2];
  int8_t *ori = new int8_t[n * 2];
  for (size_t i = 0; i < n; i++) {
    coor[i * 2]     = p_arr[i].x;
    coor[i * 2 + 1] = p_arr[i].y;
    ori[i * 2]      = p_arr[i].ux;
    ori[i * 2 + 1]  = p_arr[i].uy;
  }
  put_data(coor, ori);
  delete[] coor;
  delete[] ori;
}

/*************************************************************************//**
 * \brief  Profile exporter
 ****************************************************************************/
class ProfileExporter: public BaseExporter {
public:
  explicit ProfileExporter(const cmdline::parser &cmd);
  ~ProfileExporter();

  void write_frame(int i_step,
                   const std::vector<unsigned short> &thickness_profile,
                   const std::vector<unsigned short> &num_profile,
                   const std::vector<double> &packing_frac);
  /**
  * \brief  Cal wetting profile and order parameters
  * \param cell                Cell lattices
  * \param thickness_profile   Thickness profile
  * \param num_profile         Particle number profile
  * \param packing_frac        Fraction of wetting sites
  */
  void cal_profile(const std::vector<uint8_t> &cell,
                   std::vector<unsigned short> &thickness_profile,
                   std::vector<unsigned short> &num_profile,
                   std::vector<double> &packing_frac) const;
  size_t get_row_len() const { return row_len_; }
private:
  void set_chunk_and_deflate() const;
  int deflate_level_;

  int ncid_;
  int time_id_;
  int wetting_frac_id_;
  int thickness_profile_id_;
  int num_profile_id_;

  size_t row_len_;
  size_t frame_len_;
  size_t time_idx_[1];

  std::ofstream fout_;
};

/***********************************************************************//**
 * \brief Setting for output
 * \param cmd             Cmdline parser
 * \param log_ex          Log exporter
 * \param snap_ex         Snapshot exporter
 * \param pf_ex           Exporter for wetting data
 ***************************************************************************/
void set_output_2(const cmdline::parser &cmd, LogExporter **log_ex,
                  SnapExporter_2 ** snap_ex, ProfileExporter **pf_ex);

/***********************************************************************//**
 * \brief Output data at each proper time step
 * \tparam TPar        Tempalte particle
 * \param i_step       Current time step
 * \param p_arr        Particle array
 * \param cell         Cell lattice
 * \param log_ex       Log exporter
 * \param snap_ex      Snapshot exporter
 * \param pf_ex        Profile exporter
 ****************************************************************************/
template <typename TPar>
void output_2(int i_step, const std::vector<TPar> &p_arr,
              const std::vector<uint8_t> &cell,
              LogExporter *log_ex,
              SnapExporter_2 *snap_ex,
              ProfileExporter *pf_ex) {
  if (log_ex)
    log_ex->record(i_step);
  if (pf_ex && pf_ex->need_export(i_step)) {

    
  }
  if (snap_ex && snap_ex->need_export(i_step)) {
    snap_ex->write_frame(i_step, p_arr);
  }
}
} // end of namespace lattice

#endif
