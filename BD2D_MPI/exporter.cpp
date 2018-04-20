#include "exporter.h"
#include "comn.h"

BaseExporter_2::BaseExporter_2(const cmdline::parser & cmd, int mode):  // NOLINT
  frame_interval_(0), iframe_(0) {
  domain_length_.x = cmd.get<double>("Lx");
  domain_length_.y = cmd.exist("Ly")? cmd.get<double>("Ly"): domain_length_.x;
  pack_frac_ = cmd.get<double>("phi");
  sigma_ = cmd.get<double>("sigma");
  n_par_ = cal_particle_number_2(pack_frac_, domain_length_.x, domain_length_.y, sigma_);
  n_step_ = cmd.get<int>("n_step");
  h_ = cmd.get<double>("h");
  v0_ = cmd.get<double>("v0");
  seed_ = cmd.get<ull>("seed");
  tumbling_rate_ = cmd.get<double>("alpha");
  particle_hardness_ = cmd.get<double>("spring_const");
  wall_hardness_ = cmd.get<double>("k_wall");
  eps_ = cmd.get<double>("eps");
  min_pts_ = cmd.get<int>("min_pts");
  height_min_ = cmd.get<double>("height_min");

  if (cmd.exist("output")) {
    folder_ = cmd.get<std::string>("output") + delimiter;
    mkdir(folder_);
  } else {
    std::cout<< "Please specify the outputting folder!\n";
    exit(1);
  }

  if (mode == 1) {
    char t_str[100];
    const auto now = std::chrono::system_clock::now();
    auto t_now = std::chrono::system_clock::to_time_t(now);
    // ReSharper disable once CppDeprecatedEntity
    std::strftime(t_str, 100, "%F_%H-%M-%S", std::localtime(&t_now));
    path_ = folder_ + t_str + ".log";
  } else if (mode == 2) {
    path_ = folder_ + "traj.extxyz";
  }
}

void BaseExporter_2::set_lin_frame(int sep) {
  frame_interval_ = sep;
  for (auto i = sep; i <= n_step_; i+= sep) {
    frames_arr_.push_back(i);
  }
}

bool BaseExporter_2::need_export(int i_step) {
  auto flag = false;
  if (!frames_arr_.empty() && i_step == frames_arr_[iframe_]) {
    iframe_++;
    flag = true;
  }
  return flag;
}

LogExporter_2::LogExporter_2(const cmdline::parser & cmd):
                             BaseExporter_2(cmd, 1), fout_(path_) {
  t_start_ = std::chrono::system_clock::now();
  auto start_time = std::chrono::system_clock::to_time_t(t_start_);
  char str[100];
  // ReSharper disable once CppDeprecatedEntity
  std::strftime(str, 100, "%c", std::localtime(&start_time));
  {
    fout_ << "Started simulation at " << str << "\n";
    fout_ << "\n-------- Parameters --------";
    fout_ << "\nParticle number=" << n_par_;;
    fout_ << "\nPacking fraction=" << pack_frac_;
    fout_ << "\nLx=" << domain_length_.x;
    fout_ << "\nLy=" << domain_length_.y;
    fout_ << "\ntumbling_rate=" << tumbling_rate_;
    fout_ << "\nv0=" << v0_;
    fout_ << "\nh=" << h_;
    fout_ << "\nn_step=" << n_step_;;
    fout_ << "\nseed=" << seed_;
    fout_ << "\nspring constant=" << particle_hardness_;
    fout_ << "\nwall hardness=" << wall_hardness_;
    fout_ << "\nintegrate mode=" << cmd.get<int>("int_mode");
    fout_ << "\n\n-------- RUN --------";
    fout_ << "\ntime step\telapsed time" << std::endl;    
  }
  set_lin_frame(cmd.get<int>("log_dt"));

}

LogExporter_2::~LogExporter_2() {
  const auto t_now = std::chrono::system_clock::now();
  auto end_time = std::chrono::system_clock::to_time_t(t_now);
  char str[100];
  // ReSharper disable CppDeprecatedEntity
  std::strftime(str, 100, "%c", std::localtime(&end_time));
  // ReSharper restore CppDeprecatedEntity
  fout_ << "Finished simulation at " << str << "\n";
  std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
  fout_ << "speed=" << n_step_ * n_par_ / elapsed_seconds.count()
        << " particle time step per seconds\n";
}

void LogExporter_2::record(int i_step) {
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

XyExporter::XyExporter(const cmdline::parser& cmd):
  BaseExporter_2(cmd, 2), fout_(path_) {
  set_lin_frame(cmd.get<int>("snap_dt"));
}

void XyExporter::write_head_lines(int n_par, int i_step) {
  char str[100];
  snprintf(str, 100, "Lattice=\" %g 0 0 0 %g 0 0 0 %g \"",
           domain_length_.x, domain_length_.y, 1.);
  fout_ << n_par << "\n";
  fout_ << str;
  fout_ << "Properties=species:S:1:pos:R:2 " << "Time=" << i_step * h_;
}
