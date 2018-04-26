 /**
 * @brief Exporter of data for 2D run-and-tumble swimmers
 * 
 * @file exporter.cpp
 * @author skyline-nju
 * @date 2018-04-26
 */
#ifndef EXPORT_H
#define EXPORT_H
#include <fstream>
#include <chrono>
#include "cmdline.h"
#include "vect.h"
#include "comn.h"
#include "wetting.h"

/*************************************************************************//**
 * \brief Log exporter
 ****************************************************************************/
class LogExporter: public BaseExporter {
public:
  explicit LogExporter(const cmdline::parser &cmd);

  ~LogExporter();

  void record(int i_step);

private:
  std::chrono::time_point<std::chrono::system_clock> t_start_;
  std::ofstream fout_;
};

/**
* \brief Transform the 2D particle's coordinates into a 3D array by setting z=1
* \tparam TPar Template type for inputing particles
* \tparam T    Type of outputing data.
* \param p_arr Particle array as input
* \param coord Coordinate array as output
*/
template <typename TPar, typename T>
void par2_to_coord3(const std::vector<TPar> &p_arr, std::vector<T> &coord) {
  auto n = p_arr.size();
  coord.reserve(n * 3);
  for (size_t i = 0; i < n; i++) {
    coord.push_back(p_arr[i].x);
    coord.push_back(p_arr[i].y);
    coord.push_back(1);
  }
}

/*************************************************************************//**
* \brief Particle snapshot exporter with netCDF format
****************************************************************************/
class NcParExporter_2 :public BaseExporter {
public:
  explicit NcParExporter_2(const cmdline::parser &cmd);

  ~NcParExporter_2();

  void open(const cmdline::parser &cmd);
  void set_chunk_and_deflate() const;
  void put_time_step(int i_step) const;
  void put_cell_lengths() const;
  void put_coordinates(const float *data) const;
  void put_atom_types(const char *data)const;
  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr);
  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr,
                   const std::vector<char> &p_type);
protected:
  /* id for each variables */
  int ncid_;
  int spatial_id_;
  int cell_spatial_id_;
  int time_id_;
  int cell_lengths_id_;
  int coordinates_id_;
  int atom_types_id_;

  /* dimension lengths */
  const size_t frame_len_;
  const size_t spatial_len_;
  const size_t atom_len_;
  const size_t cell_spatial_len_;

private:
  bool flag_nc4_;
  bool atom_types_on_;
  int deflate_level_;
  size_t time_idx_[1];
  float cell_lengths_data_[3];
};

template <typename TPar>
void NcParExporter_2::write_frame(int i_step, const std::vector<TPar>& p_arr) {
  put_time_step(i_step);
  put_cell_lengths();

  std::vector<float> coor_data;
  par2_to_coord3(p_arr, coor_data);
  put_coordinates(&coor_data[0]);
  if (atom_types_on_) {
    std::vector<char> atom_types_data(p_arr.size(), 1);
    put_atom_types(&atom_types_data[0]);
  }
  time_idx_[0] = time_idx_[0] + 1;
}

template<typename TPar>
void NcParExporter_2::write_frame(int i_step, const std::vector<TPar>& p_arr,
                                  const std::vector<char>& p_type) {
  put_time_step(i_step);
  put_cell_lengths();

  std::vector<float> coor_data;
  par2_to_coord3(p_arr, coor_data);
  put_coordinates(&coor_data[0]);
  if (atom_types_on_) {
    put_atom_types(&p_type[0]);
  }
  time_idx_[0] = time_idx_[0] + 1;
}


/*************************************************************************//**
* \brief Wetting data exporter
****************************************************************************/
class ProfileExporter : public BaseExporter  {  // NOLINT
public:
  explicit ProfileExporter(const cmdline::parser &cmd);
  ~ProfileExporter();

  void dump_frame(int i_step,
                  const std::vector<float> &thickness_profile,
                  const std::vector<unsigned short> &num_profile,
                  const std::vector<double> &packing_frac);

  template <typename TriFunc>
  void write_frame(int i_step, TriFunc cal_profile);

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

template <typename TriFunc>
void ProfileExporter::write_frame(int i_step, TriFunc cal_profile) {
  std::vector<float> thickness_profile(row_len_ * 2, 0);
  std::vector<unsigned short> num_profile(row_len_ * 2, 0);
  std::vector<double> packing_frac(2, 0);
  cal_profile(thickness_profile, num_profile, packing_frac);
  dump_frame(i_step, thickness_profile, num_profile, packing_frac);
}

/*************************************************************************//**
 * \brief Setting for output
 * \param cmd           Cmdline parser
 * \param log_exporter  Log exporter
 * \param nc_exporter   Netcdf format snapshot exporter
 * \param pf_exporter   Exporter for data about wetting
 ****************************************************************************/
void set_output_2(const cmdline::parser &cmd,
                  LogExporter **log_exporter,
                  NcParExporter_2 **nc_exporter,
                  ProfileExporter **pf_exporter);

void get_profile_para(double &eps_out, int &min_pts_out, double &h_thres_out);

/*************************************************************************//**
 * \brief Output data at each time step
 * \tparam TPar         Template of particle
 * \tparam BiFunc       Template binary function
 * \param i_step        Time step
 * \param p_arr         Array of paticles
 * \param f_dis_square  Function to cal distance between two particles
 * \param domain_len    Domain length
 * \param log_ex        Log exporter
 * \param nc_ex         Snapshot exporter
 * \param pf_ex         Profile exporter
 ****************************************************************************/
template <typename TPar, typename BiFunc>
void output_2(int i_step,
              const std::vector<TPar> &p_arr,
              BiFunc f_dis_square,
              const Vec_2<double> &domain_len,
              LogExporter *log_ex,
              NcParExporter_2 *nc_ex,
              ProfileExporter *pf_ex) {
  if (log_ex)
    log_ex->record(i_step);
  if (pf_ex) {
    std::vector<Cluster_w_xlim> c_arr;
    std::vector<bool> flag_clustered;
    std::vector<char> flag_wetting;
    auto cal_cluster = [f_dis_square, &p_arr, &domain_len](auto &cluster, auto &flag_c, auto &flag_w) {
      double eps, h_thres;
      int min_pts;
      get_profile_para(eps, min_pts, h_thres);
      dbscan_wall(cluster, flag_c, flag_w, eps, min_pts, h_thres, p_arr,
                  f_dis_square, domain_len, Vec_2<double>());
    };
    if (pf_ex->need_export(i_step)) {
      cal_cluster(c_arr, flag_clustered, flag_wetting);
      auto lambda = [&flag_wetting, &p_arr, &domain_len](auto &thickness, auto &num, auto &frac) {
        cal_wetting_profile(thickness, num, frac, p_arr, flag_wetting,
                            domain_len, Vec_2<double>());
      };
      pf_ex->write_frame(i_step, lambda);
    }
    if (nc_ex && nc_ex->need_export(i_step)) {
      if (c_arr.empty())
        cal_cluster(c_arr, flag_clustered, flag_wetting);
      nc_ex->write_frame(i_step, p_arr, flag_wetting);
    }
  } else {
    if (nc_ex && nc_ex->need_export(i_step))
      nc_ex->write_frame(i_step, p_arr);
  }
}

#endif
