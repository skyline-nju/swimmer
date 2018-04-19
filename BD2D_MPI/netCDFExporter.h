#ifndef NETCDFEXPORTER_H
#define NETCDFEXPORTER_H
#include "exporter.h"

void check_err(const int stat, const int line, const char * file);

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

/* Export particle data into netCDF file*/
class NcParExporter_2 :public BaseExporter_2 {
public:
  explicit NcParExporter_2(const cmdline::parser &cmd);

  ~NcParExporter_2();

  void open(const cmdline::parser &cmd);
  void set_chunk_and_deflate();
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
  if (need_export(i_step)) {
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

#endif
