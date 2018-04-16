#ifndef NETCDFEXPORTER_H
#define NETCDFEXPORTER_H
#include "exporter.h"

void check_err(const int stat, const int line, const char *file);

template <typename TPar, typename T>
void par2d_to_coord3d(const std::vector<TPar> &p_arr, T* coord) {
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
  
  void write_step(int i_step);
  void write_coordinates(float *coordinates_data);
  template <typename TPar>
  void write_frame(int i_step, const std::vector<TPar> &p_arr);

protected:
  int ncid_;
  int spatial_id_;
  int cell_spatial_id_;
  int time_id_;
  int coordinates_id_;
  int cell_lengths_id_;
  int atom_types_id_;
private:
  size_t time_idx_[1];
  size_t coordinates_startset_[3];
  size_t coordinates_countset_[3];
  size_t cell_lengths_startset_[2];
  size_t cell_lengths_countset_[2];
  size_t atom_types_startset_[2];
  size_t atom_types_countset_[2];
  float cell_lengths_data_[3];
};

template <typename TPar>
void NcParExporter_2::write_frame(int i_step, const std::vector<TPar>& p_arr) {
  if (need_export(i_step)) {
    auto *coor_data = new float[n_par_ * 3];
    par2d_to_coord3d(p_arr, coor_data);
    write_step(i_step);
    write_coordinates(coor_data);
    time_idx_[0] = time_idx_[0] + 1;
    delete[] coor_data;
  }
}

#endif
