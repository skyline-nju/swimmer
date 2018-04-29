/**
 * @brief lattice domain
 * 
 * @file latticeExporter.cpp
 * @author skyline-nju
 * @date 2018-04-27
 */
#include "latticeExporter.h"
#include <netcdf.h>

int lx;
int ly;
int lz;
double pack_frac;
size_t n_par;
int n_step;
unsigned long long seed;
int max_capacity;
bool rt_on;
double tumbling_rate;
double nu_f;
double nu_b;
double nu_t;
double D_rot;
std::string folder;
std::string base_name;

/***********************************************************************//**
* \brief Setting for output
* \param cmd             Cmdline parser
* \param log_ex          Log exporter
* \param snap_ex         Snapshot exporter
* \param pf_ex           Exporter for wetting data
* \param is_rt           Is run-and-tumble case ?
***************************************************************************/
void lattice::set_output_2(const cmdline::parser& cmd, LogExporter** log_ex,
                           SnapExporter_2** snap_ex, ProfileExporter **pf_ex,
                           bool is_rt) {
  if (cmd.exist("output")) {
    lx = cmd.get<int>("Lx");
    ly = cmd.exist("Ly") ? cmd.get<int>("Ly") : lx;
    lz = 0;

    pack_frac = cmd.get<double>("pack_frac");
    n_par = lx * ly * pack_frac;
    n_step = cmd.get<unsigned int>("n_step");
    seed = cmd.get<unsigned long long>("seed");
    max_capacity = cmd.get<int>("n_max");
    if (is_rt) {
      rt_on = true;
      tumbling_rate = cmd.get<double>("alpha");
    } else {
      rt_on = false;
      nu_f = cmd.get<double>("nu_f");
      nu_b = cmd.get<double>("nu_b");
      nu_t = cmd.get<double>("nu_t");
      D_rot = cmd.get<double>("D_rot");
    }

    folder = cmd.get<std::string>("output") + delimiter;
    mkdir(folder);
    char str[100];
    if (is_rt) {
      snprintf(str, 100, "rt_%g_%g_%d", tumbling_rate, pack_frac, max_capacity);
      base_name = str;
    }
    else {
      snprintf(str, 100, "ab_%g_%g_%d", nu_f, pack_frac, max_capacity);
      base_name = str;
    }

    if (cmd.exist("log_dt"))
      *log_ex = new LogExporter(cmd);
    else
      *log_ex = nullptr;
    if (cmd.exist("snap_dt"))
      *snap_ex = new SnapExporter_2(cmd);
    else
      *snap_ex = nullptr;
    if (cmd.exist("profile_dt"))
      *pf_ex = new ProfileExporter(cmd);
    else
      *pf_ex = nullptr;
  }
}

/*************************************************************************//**
* \brief  Cal wetting profile and order parameters
* \param cell                Cell lattices
* \param thickness_profile   Thickness profile
* \param num_profile         Particle number profile
* \param packing_frac        Fraction of wetting sites
*****************************************************************************/
void lattice::cal_profile(const std::vector<uint8_t>& cell,
                          std::vector<uint16_t>& thickness_profile,
                          std::vector<uint16_t>& num_profile,
                          std::vector<double>& packing_frac) {
  for (auto row = 0; row < ly; row++) {
    const size_t row_lx = row * lx;
    const auto i_lt = row * 2;
    const auto i_rt = i_lt + 1;
    for (auto col = 0; col < lx; col++) {
      const auto ic = col + row_lx;
      if (cell[ic]) {
        ++thickness_profile[i_lt];
        num_profile[i_lt] += cell[ic];
      } else {
        if (thickness_profile[i_lt]) {
          packing_frac[0] += 1;
        }
        break;
      }
    }
    for (auto col = lx - 1; col >= 0; col--) {
      const auto ic = col + row_lx;
      if (cell[ic]) {
        ++thickness_profile[i_rt];
        num_profile[i_rt] += cell[ic];
      } else {
        if (thickness_profile[i_rt]) {
          packing_frac[1] += 1;
        }
        break;
      }
    }
  }
  packing_frac[0] /= ly;
  packing_frac[1] /= ly;
}

/*************************************************************************//**
* \brief Constructor of LogExporter
* \param cmd Cmdline parser
*****************************************************************************/
lattice::LogExporter::LogExporter(const cmdline::parser& cmd)
: BaseLogExporter(folder + base_name + ".log", n_par,
                  n_step, cmd.get<int>("log_dt")) {
  fout_ << "\n-------- Parameters --------";
  fout_ << "\nParticle number=" << n_par;
  fout_ << "\nPacking fraction=" << pack_frac;
  fout_ << "\nLx=" << lx;
  fout_ << "\nLy=" << ly;
  fout_ << "\nLz=" << lz;
  if (rt_on) {
    fout_ << "\ntumbling rate=" << tumbling_rate;
  } else {
    fout_ << "\nforward rate=" << nu_f;
    fout_ << "\nbackward rate=" << nu_b;
    fout_ << "\ntransverse rate=" << nu_t;
    fout_ << "\nrotational diffusion rate=" << D_rot;
  }
  fout_ << "\nn_step=" << n_step;;
  fout_ << "\nseed=" << seed;
  fout_ << "\nmax_capacity=" << max_capacity;
  fout_ << "\n\n-------- RUN --------";
  fout_ << "\ntime step\telapsed time" << std::endl;
}

void check_err(const int stat, const int line, const char * file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

/*************************************************************************//**
 * \brief Constructor of SnapExporter
 * \param cmd Cmdline parser
 *****************************************************************************/
lattice::SnapExporter_2::SnapExporter_2(const cmdline::parser& cmd) // NOLINT
  : frame_len_(NC_UNLIMITED), spatial_len_(2), atom_len_(n_par),
    cell_spatial_len_(2), time_idx_{0} {
  deflate_level_ = 6;
  set_lin_frame(cmd.get<int>("snap_dt"), n_step);
  open(cmd);
}

lattice::SnapExporter_2::~SnapExporter_2() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void lattice::SnapExporter_2::open(const cmdline::parser& cmd) {
  char str[100];
  snprintf(str, 100, "%straj_%s.nc", folder.c_str(), base_name.c_str());
  auto stat = nc_create(str, NC_NETCDF4, &ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* dimension ids */
  int frame_dim;
  int spatial_dim;
  int atom_dim;
  int cell_spatial_dim;

  /* define dimensions */
  stat = nc_def_dim(ncid_, "frame", frame_len_, &frame_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "spatial", spatial_len_, &spatial_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "atom", atom_len_, &atom_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid_, "cell_spatial", cell_spatial_len_, &cell_spatial_dim);
  check_err(stat, __LINE__, __FILE__);

  /* define variables */
  int spatial_dims[1] = {spatial_dim};
  stat = nc_def_var(ncid_, "spatial", NC_CHAR, 1, spatial_dims, &spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  int cell_spatial_dims[1] = {cell_spatial_dim};
  stat = nc_def_var(ncid_, "cell_spatial", NC_CHAR, 1, cell_spatial_dims, &cell_spatial_id_);
  check_err(stat, __LINE__, __FILE__);

  int cell_lengths_dims[1] = {cell_spatial_dim};
  stat = nc_def_var(ncid_, "cell_lengths", NC_INT, 1, cell_lengths_dims, &cell_lengths_id_);
  check_err(stat, __LINE__, __FILE__);

  int time_dims[1] = {frame_dim};
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);

  int coordinates_dims[3] = {frame_dim, atom_dim, spatial_dim};
  stat = nc_def_var(ncid_, "coordinates", NC_INT, 3, coordinates_dims, &coordinates_id_);
  check_err(stat, __LINE__, __FILE__);

  int orientaions_dims[3] = {frame_dim, atom_dim, spatial_dim};
  stat = nc_def_var(ncid_, "orientations", NC_BYTE, 3, orientaions_dims, &orientations_id_);
  check_err(stat, __LINE__, __FILE__);

  set_chunk_and_deflate();

  /* assign global attributes */
  stat = nc_put_att_text(ncid_, NC_GLOBAL, "title", 21, "lattice active matter");
  check_err(stat, __LINE__, __FILE__);

  stat = nc_put_att_double(ncid_, NC_GLOBAL, "packing_frac", NC_DOUBLE, 1, &pack_frac);
  check_err(stat, __LINE__, __FILE__);
  if (rt_on) {
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "tumbling_rate", NC_DOUBLE, 1, &tumbling_rate);
    check_err(stat, __LINE__, __FILE__);
  } else {
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "forward_rate", NC_DOUBLE, 1, &nu_f);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "backward_rate", NC_DOUBLE, 1, &nu_b);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "transverse_rate", NC_DOUBLE, 1, &nu_t);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "rot_rate", NC_DOUBLE, 1, &D_rot);
    check_err(stat, __LINE__, __FILE__);
  }
  stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_int(ncid_, NC_GLOBAL, "max_capacity", NC_INT, 1, &max_capacity);
  check_err(stat, __LINE__, __FILE__);

  /* assign per-variable attributes */
  stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval_);
  check_err(stat, __LINE__, __FILE__);

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);

  /* assign variable data */
  stat = nc_put_var(ncid_, spatial_id_, "xyz");
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_var(ncid_, cell_spatial_id_, "abc");
  check_err(stat, __LINE__, __FILE__);
  size_t startset[1] = {0};
  size_t countset[1] = {2};
  int lengths[2] = {lx, ly};
  stat = nc_put_vara(ncid_, cell_lengths_id_, startset, countset, lengths);
  check_err(stat, __LINE__, __FILE__);
}

void lattice::SnapExporter_2::set_chunk_and_deflate() const {
  auto time_block = get_n_frames() / 10;
  if (time_block >= 1000)
    time_block = 1000;
  else if (time_block < 1)
    time_block = 1;

  /* time steps */
  {
    size_t chunk[1] = {time_block};
    auto stat = nc_def_var_chunking(ncid_, time_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, time_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }

  time_block = 48 * 4096 / 8 / n_par;
  if (time_block > get_n_frames()) {
    time_block = get_n_frames();
  } else if (time_block < 1) {
    time_block = 1;
  }

  /* coordinates and oritations */
  {
    size_t chunk[3] = {time_block, atom_len_, spatial_len_};
    auto stat = nc_def_var_chunking(ncid_, coordinates_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, coordinates_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_chunking(ncid_, orientations_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, orientations_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }
}

void lattice::SnapExporter_2::put_time_step(int i_step) const {
  const auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
  check_err(stat, __LINE__, __FILE__);
}

void lattice::SnapExporter_2::put_data(const int* coor, const int8_t* ori) const {
  size_t startset[3] = {time_idx_[0], 0, 0};
  size_t countset[3] = {1, atom_len_, spatial_len_};
  auto stat = nc_put_vara(ncid_, coordinates_id_, startset, countset, coor);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_vara(ncid_, orientations_id_, startset, countset, ori);
  check_err(stat, __LINE__, __FILE__);
}

/*************************************************************************//**
 * @brief Construct a new lattice::Profile Exporter::Profile Exporter object
 * 
 * @param cmd  Cmdline parser
 ***************************************************************************/
lattice::ProfileExporter::ProfileExporter(const cmdline::parser& cmd) // NOLINT
  : time_idx_{0}, fout_(folder + "profile_" + base_name + ".dat") {
  frame_len_ = NC_UNLIMITED;
  row_len_ = int(ly);
  deflate_level_ = 6;

  set_lin_frame(cmd.get<int>("profile_dt"), n_step);
  char str[100];
  snprintf(str, 100, "%sprofile_%s.nc", folder.c_str(), base_name.c_str());
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
  stat = nc_def_var(ncid_, "wetting_fraction", NC_DOUBLE, 2,
                    wetting_frac_dims, &wetting_frac_id_);
  check_err(stat, __LINE__, __FILE__);

  int thickness_profile_dims[3] = {frame_dim, left_right_dim, row_dim};
  stat = nc_def_var(ncid_, "thickness_profile", NC_USHORT, 3,
                    thickness_profile_dims, &thickness_profile_id_);
  check_err(stat, __LINE__, __FILE__);

  int num_profile_dims[3] = {frame_dim, left_right_dim, row_dim};
  stat = nc_def_var(ncid_, "particle_number_profile", NC_USHORT, 3,
                    num_profile_dims, &num_profile_id_);
  check_err(stat, __LINE__, __FILE__);

  set_chunk_and_deflate();

  /* write simulating parameters as global attributes */
  stat = nc_put_att_double(ncid_, NC_GLOBAL, "packing_frac", NC_DOUBLE, 1, &pack_frac);
  check_err(stat, __LINE__, __FILE__);
  if (rt_on) {
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "tumbling_rate", NC_DOUBLE, 1, &tumbling_rate);
    check_err(stat, __LINE__, __FILE__);
  } else {
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "forward_rate", NC_DOUBLE, 1, &nu_f);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "backward_rate", NC_DOUBLE, 1, &nu_b);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "transverse_rate", NC_DOUBLE, 1, &nu_t);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "rot_rate", NC_DOUBLE, 1, &D_rot);
    check_err(stat, __LINE__, __FILE__);
  }
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_put_att_int(ncid_, NC_GLOBAL, "max_capacity", NC_INT, 1, &max_capacity);
  check_err(stat, __LINE__, __FILE__);

  /* assign per-variable attributes */
  stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval_);
  check_err(stat, __LINE__, __FILE__);

  /* leave define mode */
  stat = nc_enddef(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

lattice::ProfileExporter::~ProfileExporter() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
  fout_.close();
}

void lattice::ProfileExporter::write_frame(int i_step,
                                           const std::vector<uint16_t>& thickness_profile,
                                           const std::vector<uint16_t>& num_profile,
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
    fout_ << i_step << std::fixed << std::setprecision(8) << "\t"
          << packing_frac[0] << "\t" << packing_frac[1] << std::endl;
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
  ++time_idx_[0]; //! update time frame
}

void lattice::ProfileExporter::write_frame(int i_step,
                                           const std::vector<uint8_t>& cell) {
  std::vector<uint16_t> thickness_profile(ly * 2, 0);
  std::vector<uint16_t> num_profile(ly * 2, 0);
  std::vector<double> packing_frac(2, 0.);
  cal_profile(cell, thickness_profile, num_profile, packing_frac);
  write_frame(i_step, thickness_profile, num_profile, packing_frac);
}

void lattice::ProfileExporter::set_chunk_and_deflate() const {
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
    auto stat = nc_def_var_chunking(ncid_, wetting_frac_id_, NC_CHUNKED, chunk);
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
    auto stat = nc_def_var_chunking(ncid_, thickness_profile_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, thickness_profile_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);

    stat = nc_def_var_chunking(ncid_, num_profile_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, num_profile_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }
}
