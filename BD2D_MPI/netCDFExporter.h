#ifndef NETCDFEXPORTER_H
#define NETCDFEXPORTER_H
#include "exporter.h"

/* Transform the 2D particle's coordinates into a 3D array by setting z=1 */
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

class NcParExporter_2 :public BaseExporter_2 {
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
    std::vector<char> atom_types_data(n_par_, 1);
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

/* Export the data about wetting transition */
class ProfileExporter : public BaseExporter_2 {  // NOLINT
public:
  explicit ProfileExporter(const cmdline::parser &cmd);
  ~ProfileExporter();

  void dump_frame(int i_step,
                  const std::vector<float> &thickness_profile,
                  const std::vector<unsigned short> &num_profile,
                  const std::vector<double> &packing_frac);

  //template <typename TPar, typename TBc>
  //void write_frame(int i_step, const std::vector<TPar> &p_arr,
  //                 const std::vector<char> &flag_wetting, const TBc &bc);

  template <typename TriFunc>
  void write_frame(int i_step, TriFunc cal_profile);

  double get_eps() const { return eps_; }
  unsigned int get_min_pts() const { return min_pts_; }
  double get_height_min() const { return height_min_; }
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


//template<typename TPar, typename TBc>
//void ProfileExporter::write_frame(int i_step, const std::vector<TPar>& p_arr,
//                                  const std::vector<char>& flag_wetting,
//                                  const TBc & bc) {
//  std::vector<float> thickness_profile(row_len_ * 2, 0);
//  std::vector<unsigned short> num_profile(row_len_ * 2, 0);
//  std::vector<double> packing_frac(2, 0);
//  cal_wetting_profile(thickness_profile, num_profile, packing_frac,
//                      p_arr, flag_wetting, bc);
//  dump_frame(i_step, thickness_profile, num_profile, packing_frac);
//}

template <typename TriFunc>
void ProfileExporter::write_frame(int i_step, TriFunc cal_profile) {
  std::vector<float> thickness_profile(row_len_ * 2, 0);
  std::vector<unsigned short> num_profile(row_len_ * 2, 0);
  std::vector<double> packing_frac(2, 0);
  cal_profile(thickness_profile, num_profile, packing_frac);
  dump_frame(i_step, thickness_profile, num_profile, packing_frac);
}

#endif
