 /**
 * @brief Exporter of data for 2D run-and-tumble swimmers
 * 
 * @file exporter.cpp
 * @author skyline-nju
 * @date 2018-04-26
 */
#include "exporter.h"
#include <cstdio>
#include <cstdlib>
#include <netcdf.h>
#include <iomanip>

double lx;
double ly;
double lz;
double pack_frac;
double sigma;
size_t n_par;
int n_step;
double h;
double v0;
unsigned long long seed;
double tumbling_rate;
double particle_hardness;
double wall_hardness;
std::string folder;
std::string base_name;
double eps;
int min_pts;
double height_thres;

/*************************************************************************//**
 * \brief Setting for output
 * \param cmd           Cmdline parser
 * \param log_exporter  Log exporter
 * \param nc_exporter   Netcdf format snapshot exporter
 * \param pf_exporter   Exporter for data about wetting
 ****************************************************************************/
void set_output_2(const cmdline::parser& cmd,
                  LogExporter **log_exporter,
                  NcParExporter_2 ** nc_exporter,
                  ProfileExporter **pf_exporter) {
  if (cmd.exist("output")) {
    lx = cmd.get<double>("Lx");
    ly = cmd.exist("Ly") ? cmd.get<double>("Ly") : lx;
    lz = 0;

    pack_frac = cmd.get<double>("phi");
    sigma = cmd.get<double>("sigma");
    n_par = cal_particle_number_2(pack_frac, lx, ly, sigma);
    n_step = cmd.get<int>("n_step");
    h = cmd.get<double>("h");
    v0 = cmd.get<double>("v0");
    seed = cmd.get<unsigned long long>("seed");
    tumbling_rate = cmd.get<double>("alpha");
    particle_hardness = cmd.get<double>("spring_const");
    wall_hardness = cmd.get<double>("k_wall");
    eps = cmd.get<double>("eps");
    min_pts = cmd.get<int>("min_pts");
    height_thres = cmd.get<double>("height_min");

    if (cmd.exist("output")) {
      folder = cmd.get<std::string>("output") + delimiter;
      mkdir(folder);
    } else {
      std::cout << "Please specify the outputting folder!\n";
      exit(1);
    }

    char str[100];
    snprintf(str, 100, "%g_%g_%g", tumbling_rate, particle_hardness, pack_frac);
    base_name = str;

    if (cmd.exist("log_dt"))
      *log_exporter = new LogExporter(cmd);
    else
      *log_exporter = nullptr;
    if (cmd.exist("snap_dt"))
      *nc_exporter = new NcParExporter_2(cmd);
    else
      *nc_exporter = nullptr;
    if (cmd.exist("profile_dt"))
      *pf_exporter = new ProfileExporter(cmd);
    else
      *pf_exporter = nullptr;
  }
}

void get_profile_para(double& eps_out, int& min_pts_out,
                      double& h_thres_out) {
  eps_out = eps;
  min_pts_out = min_pts;
  h_thres_out = height_thres;
}

/*************************************************************************//**
 * \brief Constructor of LogExporter
 * \param cmd Cmdline parser
 ***************************************************************************/
LogExporter::LogExporter(const cmdline::parser& cmd)
  : BaseExporter(n_step, cmd.get<int>("log_dt")),
    fout_(folder + base_name + ".log") {

  t_start_ = std::chrono::system_clock::now();
  auto start_time = std::chrono::system_clock::to_time_t(t_start_);
  char str[100];
  // ReSharper disable CppDeprecatedEntity
  std::strftime(str, 100, "%c", std::localtime(&start_time));
  // ReSharper restore CppDeprecatedEntity

  fout_ << "Started simulation at " << str << "\n";
  fout_ << "\n-------- Parameters --------";
  fout_ << "\nParticle number=" << n_par;;
  fout_ << "\nPacking fraction=" << pack_frac;
  fout_ << "\nLx=" << lx;
  fout_ << "\nLy=" << ly;
  fout_ << "\nLz=" << lz;
  fout_ << "\ntumbling_rate=" << tumbling_rate;
  fout_ << "\nv0=" << v0;
  fout_ << "\nh=" << h;
  fout_ << "\nn_step=" << n_step;;
  fout_ << "\nseed=" << seed;
  fout_ << "\nspring constant=" << particle_hardness;
  fout_ << "\nwall hardness=" << wall_hardness;
  fout_ << "\nintegrate mode=" << cmd.get<int>("int_mode");
  fout_ << "\n\n-------- RUN --------";
  fout_ << "\ntime step\telapsed time" << std::endl;
}

LogExporter::~LogExporter() {
  const auto t_now = std::chrono::system_clock::now();
  auto end_time = std::chrono::system_clock::to_time_t(t_now);
  char str[100];
  // ReSharper disable CppDeprecatedEntity
  std::strftime(str, 100, "%c", std::localtime(&end_time));
  // ReSharper restore CppDeprecatedEntity
  fout_ << "Finished simulation at " << str << "\n";
  std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
  fout_ << "speed=" << n_step_ * n_par / elapsed_seconds.count()
    << " particle time step per seconds\n";
  fout_.close();
}

void LogExporter::record(int i_step) {
  if (need_export(i_step)) {
    const auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    const auto dt = elapsed_seconds.count();
    const auto hour = int(dt / 3600);
    const auto min = int((dt - hour * 3600) / 60);
    const int sec = dt - hour * 3600 - min * 60;
    fout_ << i_step << "\t" << hour << ":" << min << ":" << sec << std::endl;
  }
}

void check_err(const int stat, const int line, const char * file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

/*************************************************************************//**
 * \brief Constructor for NcParExporter_2
 * \param cmd Cmdline parser
 ****************************************************************************/
NcParExporter_2::NcParExporter_2(const cmdline::parser& cmd) // NOLINT
  : frame_len_(NC_UNLIMITED),
  spatial_len_(3), atom_len_(n_par),
  cell_spatial_len_(3), time_idx_{0} {
  flag_nc4_ = true;   // If true, use netCDF-4 format, else use 64bit-offset format
  deflate_level_ = 6;
  atom_types_on_ = true;
  set_lin_frame(cmd.get<int>("snap_dt"), n_step);
  cell_lengths_data_[0] = cmd.get<double>("Lx");
  cell_lengths_data_[1] = cmd.exist("Ly") ? cmd.get<double>("Ly") : cell_lengths_data_[0];
  cell_lengths_data_[2] = 0;
  open(cmd);
}

NcParExporter_2::~NcParExporter_2() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::open(const cmdline::parser &cmd) {
  char str[100];
  snprintf(str, 100, "%straj_%s.nc", folder.c_str(), base_name.c_str());
  int stat;
  if (flag_nc4_)
    stat = nc_create(str, NC_NETCDF4, &ncid_);
  else
    stat = nc_create(str, NC_64BIT_OFFSET, &ncid_);
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

  int time_dims[1] = {frame_dim};
  stat = nc_def_var(ncid_, "time", NC_INT, 1, time_dims, &time_id_);
  check_err(stat, __LINE__, __FILE__);

  int cell_lengths_dims[2] = {frame_dim, cell_spatial_dim};
  stat = nc_def_var(ncid_, "cell_lengths", NC_FLOAT, 2, cell_lengths_dims,
                    &cell_lengths_id_);
  check_err(stat, __LINE__, __FILE__);
  int coordinates_dims[3] = {frame_dim, atom_dim, spatial_dim};
  stat = nc_def_var(ncid_, "coordinates", NC_FLOAT, 3, coordinates_dims,
                    &coordinates_id_);
  check_err(stat, __LINE__, __FILE__);

  if (atom_types_on_) {
    int atom_types_dims[2] = {frame_dim, atom_dim};
    stat = nc_def_var(ncid_, "atom_types", NC_BYTE, 2, atom_types_dims,
                      &atom_types_id_);
    check_err(stat, __LINE__, __FILE__);
  }

  if (flag_nc4_) {
    set_chunk_and_deflate();
  }

  /* assign global attributes */
  {
    stat = nc_put_att_text(ncid_, NC_GLOBAL, "title", 18, "wetting transition");
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
  /* write simulating parameters as global attributes */
  {
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "packing_frac", NC_DOUBLE, 1, &pack_frac);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "v0", NC_DOUBLE, 1, &v0);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "h0", NC_DOUBLE, 1, &h);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "tumbling_rate", NC_DOUBLE, 1, &tumbling_rate);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "particle_hardness", NC_DOUBLE, 1, &particle_hardness);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "wall_hardness", NC_DOUBLE, 1, &wall_hardness);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
    check_err(stat, __LINE__, __FILE__);
  }
  /* assign per-variable attributes */
  {
    stat = nc_put_att_int(ncid_, time_id_, "frame_inteval", NC_INT, 1, &frame_interval_);
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

void NcParExporter_2::set_chunk_and_deflate() const {
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

  /* cell lengths */
  {
    size_t chunk[2] = {time_block, cell_spatial_len_};
    auto stat = nc_def_var_chunking(ncid_, cell_lengths_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, cell_lengths_id_, 0, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }

  time_block = 48 * 4096 / 12 / n_par;
  if (time_block > get_n_frames()) {
    time_block = get_n_frames();
  } else if (time_block < 1) {
    time_block = 1;
  }

  /* coordinates */
  {
    size_t chunk[3] = {time_block, atom_len_, spatial_len_};
    auto stat = nc_def_var_chunking(ncid_, coordinates_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, coordinates_id_, 0, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }

  /* atom types */
  if (atom_types_on_) {
    size_t chunk[2] = {time_block, atom_len_};
    auto stat = nc_def_var_chunking(ncid_, atom_types_id_, NC_CHUNKED, chunk);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_var_deflate(ncid_, atom_types_id_, NC_SHUFFLE, 1, deflate_level_);
    check_err(stat, __LINE__, __FILE__);
  }
}

void NcParExporter_2::put_time_step(int i_step) const {
  auto stat = nc_put_var1(ncid_, time_id_, time_idx_, &i_step);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::put_cell_lengths() const {
  size_t startset[2] = {time_idx_[0], 0};
  size_t countset[2] = {1, cell_spatial_len_};
  auto stat = nc_put_vara(ncid_, cell_lengths_id_, startset, countset,
                          cell_lengths_data_);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::put_coordinates(const float * data) const {
  size_t startset[3] = {time_idx_[0], 0, 0};
  size_t countset[3] = {1, atom_len_, spatial_len_};
  const auto stat = nc_put_vara(ncid_, coordinates_id_, startset,
                                countset, data);
  check_err(stat, __LINE__, __FILE__);
}

void NcParExporter_2::put_atom_types(const char * data) const {
  size_t startset[2] = {time_idx_[0], 0};
  size_t countset[2] = {1, atom_len_};
  auto stat = nc_put_vara(ncid_, atom_types_id_, startset,
                          countset, data);
  check_err(stat, __LINE__, __FILE__);
}

/*************************************************************************//**
 * \brief Constructor of ProfielExporter
 * \param cmd Cmdline parser
 ****************************************************************************/
ProfileExporter::ProfileExporter(const cmdline::parser & cmd)  // NOLINT
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
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "packing_frac", NC_DOUBLE, 1, &pack_frac);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "v0", NC_DOUBLE, 1, &v0);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "h0", NC_DOUBLE, 1, &h);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "tumbling_rate", NC_DOUBLE, 1, &tumbling_rate);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "particle_hardness", NC_DOUBLE, 1, &particle_hardness);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "wall_hardness", NC_DOUBLE, 1, &wall_hardness);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_ulonglong(ncid_, NC_GLOBAL, "seed", NC_UINT64, 1, &seed);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "eps", NC_DOUBLE, 1, &eps);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_int(ncid_, NC_GLOBAL, "min_pts", NC_INT, 1, &min_pts);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_put_att_double(ncid_, NC_GLOBAL, "height_min", NC_DOUBLE, 1, &height_thres);
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

ProfileExporter::~ProfileExporter() {
  const auto stat = nc_close(ncid_);
  check_err(stat, __LINE__, __FILE__);
  fout_.close();
}

void ProfileExporter::dump_frame(int i_step,
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
  time_idx_[0] = time_idx_[0] + 1;
}

void ProfileExporter::set_chunk_and_deflate() const {
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

