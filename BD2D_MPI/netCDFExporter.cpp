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
                                 BaseExporter_2(cmd),
                                 frame_len_(NC_UNLIMITED),
                                 spatial_len_(3),
                                 atom_len_(n_par_),
                                 cell_spatial_len_(3),
                                 time_idx_{0} {
  set_lin_frame(cmd.get<int>("snap_dt"));
  open();
  cell_lengths_data_[0] = cmd.get<double>("Lx");
  cell_lengths_data_[1] = cmd.exist("Ly")? cmd.get<double>("Ly"): cell_lengths_data_[0];
  cell_lengths_data_[2] = 0;
}

NcParExporter_2::~NcParExporter_2() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::open(bool flag_nc4) {
  char str[100];
  snprintf(str, 100, "%s%straj.nc", folder_.c_str(), delimiter.c_str());
  int stat;
  if (flag_nc4)
    stat = nc_create(str, NC_NETCDF4, &ncid_);
  else
    stat = nc_create(str, NC_64BIT_OFFSET, &ncid_);
  check_err(stat, __LINE__, __FILE__);

  set_dims_vars();

  if (flag_nc4)
    set_chunk_deflate();

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
}

void NcParExporter_2::set_dims_vars() {
    /* dimension ids */
  int frame_dim;
  int spatial_dim;
  int atom_dim;
  int cell_spatial_dim;

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
  auto stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", spatial_len_, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "atom", atom_len_, &atom_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "cell_spatial", cell_spatial_len_, &cell_spatial_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  spatial_dims[0] = spatial_dim;
  stat = nc_def_var(ncid_, "spatial", NC_CHAR, rank_spatial, spatial_dims, &spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  cell_spatial_dims[0] = cell_spatial_dim;
  stat = nc_def_var(ncid_, "cell_spatial", NC_CHAR, rank_cell_spatial, cell_spatial_dims, &cell_spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  cell_lengths_dims[0] = frame_dim;
  cell_lengths_dims[1] = cell_spatial_dim;
  stat = nc_def_var(ncid_, "cell_lengths", NC_FLOAT, rank_cell_lengths, cell_lengths_dims, &cell_lengths_id_);
  check_err(stat, __LINE__, __FILE__);

  time_dims[0] = frame_dim;
  stat = nc_def_var(ncid_, "time", NC_INT, rank_time, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);

  coordinates_dims[0] = frame_dim;
  coordinates_dims[1] = atom_dim;
  coordinates_dims[2] = spatial_dim;
  stat = nc_def_var(ncid_, "coordinates", NC_FLOAT, rank_coordinates, coordinates_dims, &coordinates_id_);
  check_err(stat, __LINE__, __FILE__);

  atom_types_dims[0] = frame_dim;
  atom_types_dims[1] = atom_dim;
  stat = nc_def_var(ncid_, "atom_types", NC_BYTE, rank_atom_types, atom_types_dims, &atom_types_id_);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::set_chunk_deflate() const {
  int deflate_level = 6;
  /* set time_block_size for time and cell_lengths */
  size_t time_block_size = get_n_frames() / 10;
  size_t cell_lengths_chunks[2] = {time_block_size, cell_spatial_len_};
  auto stat = nc_def_var_chunking(ncid_, cell_lengths_id_, NC_CHUNKED, cell_lengths_chunks);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var_deflate(ncid_, cell_lengths_id_, 0, 1, deflate_level);
  check_err(stat, __LINE__, __FILE__);

  size_t time_chunks[1] = {time_block_size};
  stat = nc_def_var_chunking(ncid_, time_id_, NC_CHUNKED, time_chunks);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var_deflate(ncid_, time_id_, NC_SHUFFLE, 1, deflate_level);
  check_err(stat, __LINE__, __FILE__);

  /* set time_block_size for coordinates and atom_types */
  time_block_size = 48 * 4096 / 12 / n_par_;
  if (time_block_size > get_n_frames()) {
    time_block_size = get_n_frames();
  } else if (time_block_size < 1) {
    time_block_size = 1;
  }
  std::cout << "t_chunk = " << time_block_size << std::endl;
  size_t coordinates_chunks[3] = {time_block_size, atom_len_, spatial_len_};
  stat = nc_def_var_chunking(ncid_, coordinates_id_, NC_CHUNKED, coordinates_chunks);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var_deflate(ncid_, coordinates_id_, 0, 1, deflate_level);
  check_err(stat, __LINE__, __FILE__);
  
  size_t atom_types_chunks[2] = {time_block_size, atom_len_};
  stat = nc_def_var_chunking(ncid_, atom_types_id_, NC_CHUNKED, atom_types_chunks);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_var_deflate(ncid_, atom_types_id_, NC_SHUFFLE, 1, deflate_level);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::write_time(int i_step) {
  auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::write_static_vara() {
  {
    size_t atom_types_startset[2] = {time_idx_[0], 0};
    size_t atom_types_countset[2] = {1, n_par_};
    std::vector<char> atom_types_data(n_par_, 1);
    auto stat = nc_put_vara(ncid_, atom_types_id_, atom_types_startset,
                       atom_types_countset, &atom_types_data[0]);
    check_err(stat,__LINE__,__FILE__);
  } 
  {
    size_t cell_lengths_startset[2] = {time_idx_[0], 0};
    size_t cell_lengths_countset[2] = {1, cell_spatial_len_};
    auto stat = nc_put_vara(ncid_, cell_lengths_id_, cell_lengths_startset,
                            cell_lengths_countset, cell_lengths_data_);
    check_err(stat,__LINE__,__FILE__);
  }
}

void NcParExporter_2::write_coordinates(float* coordinates_data) {
  size_t coordinates_startset[3] = {time_idx_[0], 0, 0};
  size_t coordinates_countset[3] = {1, atom_len_, spatial_len_};
  const auto stat = nc_put_vara(ncid_, coordinates_id_, coordinates_startset,
                                coordinates_countset, coordinates_data);
  check_err(stat, __LINE__, __FILE__);
}
