#include "netCDFExporter.h"
#include "comn.h"
#include <cstdio>
#include <cstdlib>
#include <netcdf.h>

void check_err(const int stat, const int line, const char * file) {
    if (stat != NC_NOERR) {
    (void)fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

NcParExporter_2::NcParExporter_2(const cmdline::parser& cmd):  // NOLINT
                                 BaseExporter_2(cmd), time_idx_{0},
                                 coordinates_startset_{0, 0, 0},
                                 coordinates_countset_{1, size_t(n_par_), 3},
                                 cell_lengths_startset_{0, 0},
                                 cell_lengths_countset_{1, 3},
                                 atom_types_startset_{0, 0},
                                 atom_types_countset_{1, size_t(n_par_)}{
  char str[100];
  snprintf(str, 100, "%s%straj.nc", folder_.c_str(), delimiter.c_str());
  int stat = nc_create(str, NC_CLOBBER, &ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* dimension ids */
  int frame_dim;
  int spatial_dim;
  int atom_dim;
  int cell_spatial_dim;
  
  /* dimension lengths */
  const size_t frame_len = NC_UNLIMITED;
  const size_t spatial_len = 3;
  const size_t atom_len = n_par_;
  const size_t cell_spatial_len = 3;

  /* rank (number of dimensions) for each variables */
  const int rank_spatial = 1;
  const int rank_cell_spatial = 1;
  const int rank_time = 1;
  const int rank_coordinates = 3;
  const int rank_cell_lengths = 2;
  const int rank_atom_types = 2;

  /* variable shapes */
  int spatial_dims[rank_spatial];
  int cell_spatial_dims[rank_cell_spatial];
  int time_dims[rank_time];
  int coordinates_dims[rank_coordinates];
  int cell_lengths_dims[rank_cell_lengths];
  int atom_types_dims[rank_atom_types];

  /* define dimensions */
  stat = nc_def_dim(ncid_, "frame", frame_len, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", spatial_len, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "atom", atom_len, &atom_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "cell_spatial", cell_spatial_len, &cell_spatial_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  spatial_dims[0] = spatial_dim;
  stat = nc_def_var(ncid_, "spatial", NC_CHAR, rank_spatial, spatial_dims, &spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  cell_spatial_dims[0] = cell_spatial_dim;
  stat = nc_def_var(ncid_, "cell_spatial", NC_CHAR, rank_cell_spatial, cell_spatial_dims, &cell_spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  time_dims[0] = frame_dim;
  stat = nc_def_var(ncid_, "time", NC_INT, rank_time, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);

  coordinates_dims[0] = frame_dim;
  coordinates_dims[1] = atom_dim;
  coordinates_dims[2] = spatial_dim;
  stat = nc_def_var(ncid_, "coordinates", NC_FLOAT, rank_coordinates, coordinates_dims, &coordinates_id_);
  check_err(stat, __LINE__, __FILE__);

  cell_lengths_dims[0] = frame_dim;
  cell_lengths_dims[1] = cell_spatial_dim;
  stat = nc_def_var(ncid_, "cell_lengths", NC_FLOAT, rank_cell_lengths, cell_lengths_dims, &cell_lengths_id_);
  check_err(stat, __LINE__, __FILE__);

  atom_types_dims[0] = frame_dim;
  atom_types_dims[1] = atom_dim;
  stat = nc_def_var(ncid_, "atom_types", NC_BYTE, rank_atom_types, atom_types_dims, &atom_types_id_);
  check_err(stat, __LINE__, __FILE__);

  /* assign global attributes */
  {
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "title", 18, "netCDF output test");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "application", 5, "AMBER");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "program", 6, "sander");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "programVersion", 3, "9.0");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", 5, "AMBER");
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "ConventionVersion", 3, "1.0");
    check_err(stat, __LINE__, __FILE__);
  }
    /* assign per-variable attributes */
  {
    char tstr[50];
    snprintf(tstr, 50, "%.5f tau", h_);
    stat = nc_put_att_text(ncid_, time_id_, "units", 11, tstr);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_text(ncid_, coordinates_id_, "units", 5, "sigma");
    check_err(stat, __LINE__, __FILE__);
  }

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* assign variable data */
  stat = nc_put_var(ncid_, spatial_id_, "xyz");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_var(ncid_, cell_spatial_id_, "abc");
  check_err(stat, __LINE__, __FILE__);
  
  set_lin_frame(cmd.get<int>("snap_dt"));
  cell_lengths_data_[0] = cmd.get<double>("Lx");
  cell_lengths_data_[1] = cmd.exist("Ly")? cmd.get<double>("Ly"): cell_lengths_data_[0];
  cell_lengths_data_[2] = 0;
}

NcParExporter_2::~NcParExporter_2() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}


void NcParExporter_2::write_step(int  i_step) {
  auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
  check_err(stat, __LINE__, __FILE__);
  cell_lengths_startset_[0] = time_idx_[0];
  stat = nc_put_vara(ncid_, cell_lengths_id_, cell_lengths_startset_,
                     cell_lengths_countset_, cell_lengths_data_);
  check_err(stat,__LINE__,__FILE__);

  atom_types_startset_[0] = time_idx_[0];
  std::vector<char> atom_types_data(n_par_, 1);
  //std::vector<int> atom_types_data(n_par_, 1);
  stat = nc_put_vara(ncid_, atom_types_id_, atom_types_startset_,
                     atom_types_countset_, &atom_types_data[0]);
  check_err(stat,__LINE__,__FILE__);
}

void NcParExporter_2::write_coordinates(float* coordinates_data) {
  coordinates_startset_[0] = time_idx_[0];
  const auto stat = nc_put_vara(ncid_, coordinates_id_, coordinates_startset_,
                                coordinates_countset_, coordinates_data);
  check_err(stat, __LINE__, __FILE__);
}
