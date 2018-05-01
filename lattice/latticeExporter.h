/**
 * @brief Data exporter for the lattice model
 * 
 * @file latticeExporter.h
 * @author skyline-nju
 * @date 2018-04-27
 */
#ifndef LATTICEEXPORTER_H
#define LATTICEEXPORTER_H
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cstdint>
#include "cmdline.h"
#include "comn.h"

namespace lattice {

/*************************************************************************//**
 * \brief Log exporter
 * 
 ****************************************************************************/
class LogExporter:public BaseLogExporter {
public:
  explicit LogExporter(const cmdline::parser &cmd);
};

/*************************************************************************//**
* \brief trajectory exporter
* 
****************************************************************************/
class TrajExporter_2 :public BaseExporter {
public:
  explicit TrajExporter_2(const cmdline::parser &cmd);

  ~TrajExporter_2();

  void open(const cmdline::parser &cmd);
  void set_chunk_and_deflate() const;
  void put_time_step(int i_step) const;
  void put_cell_lengths() const;
  //void put_coordinates(const float *data) const;
  void put_coordinates(const short *data) const;
  void put_atom_types(const char *data) const;
  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr,
                   const std::vector<uint8_t> &cell);
  char get_par_type(int x, int y, const std::vector<uint8_t> &cell) const;

  template <typename TPar>
  char get_par_type(TPar &p, const std::vector<uint8_t> &cell) const;
protected:
  /* id for each variables */
  int ncid_;
  int spatial_id_;
  int cell_spatial_id_;
  int time_id_;
  int cell_lengths_id_;
  int cell_origin_id_;
  int coordinates_id_;
  int atom_types_id_;

  /* dimension lengths */
  const size_t frame_len_;
  const size_t spatial_len_;
  const size_t atom_len_;
  const size_t cell_spatial_len_;

private:
  int deflate_level_;
  size_t time_idx_[1];
  int cell_lengths_data_[3];
  float cell_origin_data_[3];
  int type_mode_;
};

template<typename TPar>
void TrajExporter_2::write_frame(int i_step, const std::vector<TPar>& p_arr,
                                  const std::vector<uint8_t> &cell) {
  put_time_step(i_step);
  put_cell_lengths();

  //std::vector<float> coor_data;
  std::vector<short> coor_data;
  std::vector<char> type_arr;
  coor_data.reserve(p_arr.size() * 3);
  for (auto i : p_arr) {
    coor_data.push_back(i.x);
    coor_data.push_back(i.y);
    coor_data.push_back(0);
    type_arr.push_back(get_par_type(i.x, i.y, cell));
  }

  put_coordinates(&coor_data[0]);
  put_atom_types(&type_arr[0]);
  time_idx_[0] = time_idx_[0] + 1;
}

template <typename TPar>
char TrajExporter_2::get_par_type(TPar& p,
                                  const std::vector<uint8_t>& cell) const {
  char par_type;
  if (type_mode_==1) {
    par_type = cell[p.x + p.y * cell_lengths_data_[0]];
  } else if (type_mode_ == 2) {
    int x_new = p.x + p.ux;
    int y_new = p.y + p.uy;
    tangle_1(x_new, 0, cell_lengths_data_[0], cell_lengths_data_[0]);
    tangle_1(y_new, 0, cell_lengths_data_[1], cell_lengths_data_[1]);
    int ic_new = x_new + y_new * cell_lengths_data_[0];

    par_type = cell[ic_new] > 0 ? 1 : 0;   //? par_type is always zero
  } else {
    std::cerr << "line " << __LINE__ << " of " << __FILE__
      << ": wrong type mode" << std::endl;
    exit(1);
  }
  return par_type;
}

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
    ori[i * 2]      = p_arr[i].get_ux();
    ori[i * 2 + 1]  = p_arr[i].get_uy();
  }
  put_data(coor, ori);
  ++time_idx_[0];
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
                   const std::vector<uint16_t> &thickness_profile,
                   const std::vector<uint16_t> &num_profile,
                   const std::vector<double> &packing_frac);

  void write_frame(int i_step, const std::vector<uint8_t> &cell);

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
 * \param traj            Trajectory exporter
 * \param is_rt           Is run-and-tumble case ?
 ***************************************************************************/
void set_output_2(const cmdline::parser &cmd,
                  LogExporter **log_ex,
                  SnapExporter_2 ** snap_ex,
                  ProfileExporter **pf_ex,
                  TrajExporter_2 ** traj,
                  bool is_rt);

/***********************************************************************//**
 * \brief Output data at each proper time step
 * \tparam TPar        Tempalte particle
 * \param i_step       Current time step
 * \param p_arr        Particle array
 * \param cell         Cell lattice
 * \param log_ex       Log exporter
 * \param snap_ex      Snapshot exporter
 * \param pf_ex        Profile exporter
 * \param traj_ex      Trajectory exporter
 ****************************************************************************/
template <typename TPar>
void output_2(int i_step, const std::vector<TPar> &p_arr,
              const std::vector<uint8_t> &cell,
              LogExporter *log_ex,
              SnapExporter_2 *snap_ex,
              ProfileExporter *pf_ex,
              TrajExporter_2 *traj_ex) {
  if (log_ex)
    log_ex->record(i_step);
  if (pf_ex && pf_ex->need_export(i_step)) {
    pf_ex->write_frame(i_step, cell);
  }
  if (snap_ex && snap_ex->need_export(i_step)) {
    snap_ex->write_frame(i_step, p_arr);
  }
  if (traj_ex && traj_ex->need_export(i_step)) {
    traj_ex->write_frame(i_step, p_arr, cell);
  }
}

/*************************************************************************//**
 * \brief  Cal wetting profile and order parameters
 * \param cell                Cell lattices
 * \param thickness_profile   Thickness profile
 * \param num_profile         Particle number profile
 * \param packing_frac        Fraction of wetting sites
 *****************************************************************************/
void cal_profile(const std::vector<uint8_t> &cell,
                 std::vector<uint16_t> &thickness_profile,
                 std::vector<uint16_t> &num_profile,
                 std::vector<double> &packing_frac);

} // end of namespace lattice

#endif
