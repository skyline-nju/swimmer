#ifndef NETCDFEXPORTER_H
#define NETCDFEXPORTER_H
#include "exporter.h"

void check_err(const int stat, const int line, const char *file);

template <typename TPar, typename T>
void par2_to_coord3(const std::vector<TPar> &p_arr, T* coord) {
  auto n = p_arr.size();
  for (size_t i = 0; i < n; i++) {
    coord[i * 3] = p_arr[i].x;
    coord[i * 3 + 1] = p_arr[i].y;
    coord[i * 3 + 2] = 1;
  }    
}

class NcParExporter_2:public BaseExporter_2 {  // NOLINT
public:
  explicit NcParExporter_2(const cmdline::parser &cmd);

  ~NcParExporter_2();
  
  void open(bool flag_nc4 = true);
  void set_dims_vars();
  void set_chunk_deflate() const;
  void write_time(int i_step);
  void write_static_vara();
  void write_coordinates(float *coordinates_data);
  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr);

protected:
  /* id for each variables */
  int ncid_;
  int spatial_id_;
  int cell_spatial_id_;
  int time_id_;
  int coordinates_id_;
  int cell_lengths_id_;
  int atom_types_id_;

  /* dimension lengths */
  const size_t frame_len_;
  const size_t spatial_len_;
  const size_t atom_len_;
  const size_t cell_spatial_len_;

private:
  size_t time_idx_[1];
  float cell_lengths_data_[3];
};

template <typename TPar>
void NcParExporter_2::write_frame(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    auto *coor_data = new float[n_par_ * 3];
    par2_to_coord3(p_arr, coor_data);
    write_time(i_step);
    write_static_vara();
    write_coordinates(coor_data);
    time_idx_[0] = time_idx_[0] + 1;
    delete[] coor_data;
  }
}

#endif
