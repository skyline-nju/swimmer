#include "wetting.h"
#include "netCDFExporter.h"
#include <netcdf.h>

WetDataExporter::WetDataExporter(const cmdline::parser & cmd):  // NOLINT
  BaseExporter_2(cmd), time_idx_{0} {
  frame_len_ = NC_UNLIMITED;
  row_len_ = int(domain_length_.y);
  deflate_level_ = 6;

  set_lin_frame(cmd.get<int>("profile_dt"));
  char str[100];
  snprintf(str, 100, "%sprofile.nc", folder_.c_str());
  auto stat = nc_create(str, NC_NETCDF4, &ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* dimension ids */
  int frame_dim;
  int row_dim;
  int left_right_dim;

  /* define dimenstions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "row", row_len_, &row_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "left_right", 2, &left_right_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int time_dims[1] = {frame_dim};
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);

  int wetting_frac_dims[2] = {frame_dim, left_right_dim};
  stat = nc_def_var(ncid_, "wetting_fraction", NC_DOUBLE, 2, wetting_frac_dims, &wetting_frac_id_);
  check_err(stat, __LINE__, __FILE__);

  int thickness_profile_dims[3] = {frame_dim, left_right_dim, row_dim};
  stat = nc_def_var(ncid_, "thickness_profile", NC_FLOAT, 3, thickness_profile_dims, &thickness_profile_id_);
  check_err(stat, __LINE__, __FILE__);

  int num_profile_dims[3] = {frame_dim, left_right_dim, row_dim};
  stat = nc_def_var(ncid_, "particle_number_profile", NC_USHORT, 3, num_profile_dims, &num_profile_id_);
  check_err(stat, __LINE__, __FILE__);
  
  set_chunk_and_deflate();

  /* write simulating parameters as global attributes */
  {
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "packing_frac", NC_DOUBLE, 1, &pack_frac_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "v0", NC_DOUBLE, 1, &v0_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "h0", NC_DOUBLE, 1, &h_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "tumbling_rate", NC_DOUBLE, 1, &tumbling_rate_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "particle_hardness", NC_DOUBLE, 1, &particle_hardness_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "wall_hardness", NC_DOUBLE, 1, &wall_hardness_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "eps", NC_DOUBLE, 1, &eps_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_int(ncid_, NC_GLOBAL, "min_pts", NC_INT, 1, &min_pts_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "height_min", NC_DOUBLE, 1, &height_min_);
    check_err(stat, __LINE__, __FILE__);
  }

  /* assign per-variable attributes */
  {
    stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval_);
    check_err(stat, __LINE__, __FILE__);
  }

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

WetDataExporter::~WetDataExporter() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void WetDataExporter::dump_frame(int i_step,
                                 const std::vector<float>& thickness_profile,
                                 const std::vector<unsigned short>& num_profile,
                                 const std::vector<double>& packing_frac) {
  /* time step */
  {
    auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
    check_err(stat, __LINE__, __FILE__);
  }
  /* wetting fraction */
  {
    size_t startset[2] = {time_idx_[0], 0};
    size_t countset[2] = {1, 2};
    auto stat = nc_put_vara(ncid_, wetting_frac_id_,
                            startset, countset, &packing_frac[0]);
    check_err(stat, __LINE__, __FILE__);
  }
  /* thickness and particle number profile */
  {
    size_t startset[3] = {time_idx_[0], 0, 0};
    size_t countset[3] = {1, 2, row_len_};
    auto stat = nc_put_vara(ncid_, thickness_profile_id_,
                            startset, countset, &thickness_profile[0]);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_vara(ncid_, num_profile_id_,
                       startset, countset, &num_profile[0]);
    check_err(stat, __LINE__, __FILE__);
  }
  time_idx_[0] = time_idx_[0] + 1;
}

void WetDataExporter::set_chunk_and_deflate() const {
  /* time step */
  {
    size_t time_block = get_n_frames() / 10;
    if (time_block >= 4096)
      time_block = 4096;
    size_t chunk[1] = {time_block};
    auto stat = nc_def_var_chunking(ncid_, time_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, time_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }
  /* wetting fraction */
  {
    size_t time_block = get_n_frames() / 10;
    if (time_block >= 4096)
      time_block = 4096;
    size_t chunk[2] = {time_block, 2};
    auto stat = nc_def_var_chunking(ncid_, wetting_frac_id_,
                                    NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, time_id_, 0, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }
  /* thickness profile and particle profile */
  {
    size_t time_block = 48 * 4096 / row_len_;
    if (time_block > get_n_frames()) {
      time_block = get_n_frames();
    } else if (time_block < 1) {
      time_block = 1;
    }
    size_t chunk[3] = {time_block, 2, row_len_};
    auto stat = nc_def_var_chunking(ncid_, thickness_profile_id_,
                                    NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, thickness_profile_id_, 0, 1,
                              deflate_level_);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_def_var_chunking(ncid_, num_profile_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, num_profile_id_, NC_SHUFFLE, 1,
                              deflate_level_);
    check_err(stat, __LINE__, __FILE__);

  }
}
