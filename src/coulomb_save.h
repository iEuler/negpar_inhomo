#include <fstream>

namespace coulomb {
std::string int2str(int value, int digits = 3) {
  unsigned int uvalue = value;
  if (value < 0) {
    uvalue = -value;
  }
  std::string result;
  while (digits-- > 0) {
    result += ('0' + uvalue % 10);
    uvalue /= 10;
  }
  if (value < 0) {
    result += '-';
  }
  result = string(result.rbegin(), result.rend());
  result = '_' + result;
  return result;
}

template <class T>
void save_macro(const vector<T> &macro, string filename) {
  string filename_mod, subfix = ".txt";
  size_t Nx = macro.size();

  if (FLAG_FILENAME_WITH_NUM)
    filename_mod = "result/" + filename + int2str(K_SAVE_TIME) + subfix;
  else
    filename_mod = "result/" + filename + subfix;

  // if (FLAG_FILENAME_WITH_NUM)	filename_mod = filename +
  // int2str(K_SAVE_TIME) + subfix; else filename_mod = filename + subfix;

  ofstream file0;
  file0.open(filename_mod);
  file0 << setprecision(15);
  for (size_t kx = 0; kx < Nx; kx++) file0 << macro[kx] << '\n';
  file0.close();
}

void save_complex(int Nx, complex<double> *FS, string filename) {
  string filename_r, filename_i, subfix = ".txt";

  filename_r = "result/" + filename + "_r" + subfix;
  filename_i = "result/" + filename + "_i" + subfix;

  ofstream file0, file1;

  file0.open(filename_r);
  file1.open(filename_i);
  for (int kx = 0; kx < Nx; kx++) {
    complex<double> FSkx = *(FS + kx);
    file0 << real(FSkx) << '\n';
    file1 << imag(FSkx) << '\n';
  }
  file0.close();
  file1.close();
}

void save_2d(int Nrow, int Ncol, const vector<vector<double>> &data2d,
             string filename) {
  string filename_mod, subfix = ".txt";

  if (FLAG_FILENAME_WITH_NUM)
    filename_mod = "result/" + filename + int2str(K_SAVE_TIME) + subfix;
  else
    filename_mod = "result/" + filename + subfix;

  ofstream file0;
  file0.open(filename_mod);
  for (int krow = 0; krow < Nrow; krow++) {
    for (int kcol = 0; kcol < Ncol; kcol++) file0 << data2d[krow][kcol] << ' ';
    file0 << endl;
  }
  file0.close();
}

void save_rhouT(int Nx, const vector<double> &rho, const vector<double> &u1,
                const vector<double> &u2, const vector<double> &u3,
                const vector<double> &Tprt) {
  save_macro<double>(rho, "rho");
  save_macro<double>(u1, "u1");
  save_macro<double>(u2, "u2");
  save_macro<double>(u3, "u3");
  save_macro<double>(Tprt, "Tprt");
}

void save_rhouT(const std::vector<ParticleGroup> &Sp_x,
                const NumericGridClass &grid) {
  int Nx = grid.Nx;
  double Neff = grid.Neff;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> u2(Nx);
  vector<double> u3(Nx);
  vector<double> Tprt(Nx);

  compute_rhouT(Nx, Sp_x, Neff, rho, u1, u2, u3, Tprt);

  save_rhouT(Nx, rho, u1, u2, u3, Tprt);
}

void save_rhouT(NeParticleGroup *S_x, const NumericGridClass &grid) {
  update_rhouT(S_x, grid);
  int Nx = grid.Nx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = (S_x + kx)->rho;
    u1[kx] = (S_x + kx)->u1;
    Tprt[kx] = (S_x + kx)->Tprt;
  }

  save_macro<double>(rho, "rho");
  save_macro<double>(u1, "u1");
  save_macro<double>(Tprt, "Tprt");
}

void save_E(NeParticleGroup *S_x, const NumericGridClass &grid) {
  update_rhouT(S_x, grid);
  int Nx = grid.Nx;

  vector<double> E(Nx);
  vector<double> E_F(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    E[kx] = (S_x + kx)->elecfield;
    E_F[kx] = (S_x + kx)->elecfield_F;
  }

  save_macro<double>(E, "elecfield");
  save_macro<double>(E_F, "elecfield_F");
}

void save_rhouT_F(NeParticleGroup *S_x, const NumericGridClass &grid) {
  update_rhouT_F(S_x, grid);

  int Nx = grid.Nx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = (S_x + kx)->rhoF;
    u1[kx] = (S_x + kx)->u1F;
    Tprt[kx] = (S_x + kx)->TprtF;
  }

  save_macro<double>(rho, "rhoF");
  save_macro<double>(u1, "u1F");
  save_macro<double>(Tprt, "TprtF");
}

void save_rhouT_maxwellian(NeParticleGroup *S_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = (S_x + kx)->rhoM;
    u1[kx] = (S_x + kx)->u1M;
    Tprt[kx] = (S_x + kx)->TprtM;
  }

  save_macro<double>(rho, "rhoM");
  save_macro<double>(u1, "u1M");
  save_macro<double>(Tprt, "TprtM");
}

void save_NpNn(NeParticleGroup *S_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  vector<int> Np_all(Nx);
  vector<int> Nn_all(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    Np_all[kx] = (S_x + kx)->size('p');
    Nn_all[kx] = (S_x + kx)->size('n');
  }

  save_macro<int>(Np_all, "Np");
  save_macro<int>(Nn_all, "Nn");
}

void save_d_rhouT(NeParticleGroup *S_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = (S_x + kx)->drho;
    u1[kx] = (S_x + kx)->dm1;
    Tprt[kx] = (S_x + kx)->denergy;
  }

  save_macro<double>(rho, "drho");
  save_macro<double>(u1, "dm1");
  save_macro<double>(Tprt, "denergy");
}

void save_dx_rhouT_M(NeParticleGroup *S_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  vector<double> rho(Nx);
  vector<double> u1(Nx);
  vector<double> Tprt(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    rho[kx] = (S_x + kx)->dx_rhoM;
    u1[kx] = (S_x + kx)->dx_u1M;
    Tprt[kx] = (S_x + kx)->dx_TprtM;
  }

  save_macro<double>(rho, "dxrho");
  save_macro<double>(u1, "dxu1");
  save_macro<double>(Tprt, "dxTprt");
}

void save_m012_PN(NeParticleGroup *S_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  vector<double> m0(Nx);
  vector<double> m1(Nx);
  vector<double> m2(Nx);

  for (int kx = 0; kx < Nx; kx++) {
    (S_x + kx)->computemoments();
    m0[kx] = (S_x + kx)->m0P - (S_x + kx)->m0N;
    m1[kx] = (S_x + kx)->m11P - (S_x + kx)->m11N;
    m2[kx] = (S_x + kx)->m2P - (S_x + kx)->m2N;
  }

  save_macro<double>(m0, "m0");
  save_macro<double>(m1, "m1");
  save_macro<double>(m2, "m2");
}

void save_E(ParticleGroup *Sp_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  vector<double> elec(Nx);

  for (int kx = 0; kx < Nx; kx++) elec[kx] = (Sp_x + kx)->elecfield;

  save_macro<double>(elec, "E");
}

void save_grids(const NumericGridClass &grid) {
  vector<double> x = grid.x;
  save_macro<double>(x, "x");
  vector<double> vx = grid.vx;
  save_macro<double>(vx, "v");
}

void save_dist(ParticleGroup *Sp_x, const NumericGridClass &grid) {
  int Nx = grid.Nx;

  int Nbar = grid.Nv;
  vector<int> numinbar_one(Nbar);
  vector<vector<int>> numinbar_all(Nx, numinbar_one);

  for (int kx = 0; kx < Nx; kx++) {
    int Ndist = (Sp_x + kx)->size();
    vector<double> xdist(Ndist);
    auto &Sp = (Sp_x + kx)->list();
    for (int kp = 0; kp < Ndist; kp++) xdist[kp] = Sp[kp].velocity(0);
    histinfo_fixbar(xdist, numinbar_all[kx], grid.vmin, grid.vmax);
  }

  vector<double> dist_one(Nbar);
  vector<vector<double>> dist_all(Nx, dist_one);

  for (int kx = 0; kx < Nx; kx++) {
    for (int kbar = 0; kbar < Nbar; kbar++)
      dist_all[kx][kbar] = numinbar_all[kx][kbar] * grid.Neff;
  }

  save_2d(Nx, Nbar, dist_all, "dist");
}

void save_dist(NeParticleGroup *S_x, const NumericGridClass &grid,
               char partype) {
  int Nx = grid.Nx;
  int Nbar = grid.Nv;

  vector<int> numinbar_one(Nbar);
  vector<vector<int>> numinbar_all(Nx, numinbar_one);

  for (int kx = 0; kx < Nx; kx++) {
    int Ndist = (S_x + kx)->size(partype);
    vector<double> xdist(Ndist);
    auto &Sp = (S_x + kx)->list(partype);
    for (int kp = 0; kp < Ndist; kp++) xdist[kp] = Sp[kp].velocity(0);
    histinfo_fixbar(xdist, numinbar_all[kx], grid.vmin, grid.vmax);
  }

  double Neff = grid.Neff;
  if (partype == 'f') Neff = grid.Neff_F;
  double coeff = Neff / grid.dx / grid.dv;

  vector<double> dist_one(Nbar);
  vector<vector<double>> dist_all(Nx, dist_one);

  for (int kx = 0; kx < Nx; kx++) {
    for (int kbar = 0; kbar < Nbar; kbar++)
      dist_all[kx][kbar] = numinbar_all[kx][kbar] * coeff;
  }

  string distname;

  if (partype == 'p') {
    distname = "distp";
  } else if (partype == 'n') {
    distname = "distn";
  } else {
    distname = "distf";
  }

  save_2d(Nx, Nbar, dist_all, distname);
}

void save_dist(NeParticleGroup *S_x, const NumericGridClass &grid) {
  save_dist(S_x, grid, 'p');
  save_dist(S_x, grid, 'n');
  save_dist(S_x, grid, 'f');
}

void save_particle1d1d(ParticleGroup *Sp_x, const NumericGridClass &grid) {
  string filename_mod, subfix = ".txt", filename = "particle";

  if (FLAG_FILENAME_WITH_NUM) {
    filename_mod = "result/" + filename + int2str(K_SAVE_TIME) + subfix;
  } else {
    filename_mod = "result/" + filename + subfix;
  }

  ofstream file0;
  file0.open(filename_mod);

  int Nx = grid.Nx;
  for (int kx = 0; kx < Nx; kx++) {
    auto &Sp = (Sp_x + kx)->list();
    for (int kp = 0; kp < (Sp_x + kx)->size(); kp++)
      file0 << Sp[kp].position() << ' ' << Sp[kp].velocity(0) << ' ';
  }

  file0.close();
}

void save_particle1d1d(std::vector<NeParticleGroup> &S_x,
                       const NumericGridClass &grid, char partype,
                       int whattosave) {
  // whattosave = 1, save the first coordinates
  // whattosave = 2, save the energy

  string filename_mod, subfix = ".txt", filename;

  if (partype == 'p') {
    filename = "particleP";
  } else {
    filename = "particleN";
  }

  if (FLAG_FILENAME_WITH_NUM) {
    filename_mod = "result/" + filename + int2str(K_SAVE_TIME) + subfix;
  } else {
    filename_mod = "result/" + filename + subfix;
  }

  ofstream file0;
  file0.open(filename_mod);

  int Nx = grid.Nx;
  double savequantiti;
  for (int kx = 0; kx < Nx; kx++) {
    auto &Sp = S_x[kx].list(partype);
    for (int kp = 0; kp < S_x[kx].size(partype); kp++) {
      if (whattosave == 1) {
        savequantiti = Sp[kp].velocity(0);

      } else if (whattosave == 2) {
        savequantiti = 0.;
        auto &vp = Sp[kp].velocity();
        for (int kv = 0; kv < 3; kv++) savequantiti += vp[kv] * vp[kv];
        savequantiti = sqrt(savequantiti);
      }
      // file0 << (Sp+kp)->position()<< ' ' << savequantiti << ' ';
      file0 << savequantiti << endl;
    }
  }

  file0.close();
}

void save_particle1d1d(std::vector<NeParticleGroup> &S_x,
                       const NumericGridClass &grid) {
  save_particle1d1d(S_x, grid, 'p', 1);
  save_particle1d1d(S_x, grid, 'n', 1);
}

int count_particle_number(const std::vector<NeParticleGroup> &S_x, int Nx,
                          char partype) {
  int n = 0;
  for (int kx = 0; kx < Nx; kx++) n += S_x[kx].size(partype);
  return n;
}

void save_particleenergy(std::vector<NeParticleGroup> &S_x,
                         const NumericGridClass &grid) {
  save_particle1d1d(S_x, grid, 'p', 2);
  save_particle1d1d(S_x, grid, 'n', 2);
}

void save_homo_rdist(void) {
  // only used in homogeneous problem
  // double b = 1., c = 0.1;
  int Nr = 100;
  // double rmax = b+8.0*sqrt(c);

  double rmax = 5.;

  double dr = rmax / Nr;
  vector<double> rvec(Nr);

  for (int kr = 0; kr < Nr; kr++) {
    rvec[kr] = (kr + 0.5) * dr;
  }

  save_macro<double>(rvec, "rdist");
}

void save_homo_rdist(int Nr) {
  // only used in homogeneous problem
  // double b = 1., c = 0.1;
  // double rmax = b+8.0*sqrt(c);

  double rmax = 10.;

  double dr = rmax / Nr;
  vector<double> rvec(Nr);

  for (int kr = 0; kr < Nr; kr++) {
    rvec[kr] = (kr + 0.5) * dr;
  }

  save_macro<double>(rvec, "rdist");
}

// void save_homo_dist(NeParticleGroup * S_x, const NumericGridClass & grid, int
// flag_case) {
void save_homo_dist(const NeParticleGroup &S_x, int Nr, int flag_case) {
  // only used in homogeneous problem
  // double b = 1., c = 0.1;
  // int Nr = 400;
  // double rmax = b+8.0*sqrt(c);
  double rmax = 10.;

  vector<int> numinbar(Nr);

  auto Np = S_x.size('p');
  vector<double> pdist(Np);
  auto &Sp = S_x.list('p');
  for (int kp = 0; kp < Np; kp++) {
    auto &vf = Sp[kp].velocity();
    double vfnorm = 0.;
    for (int kv = 0; kv < 3; kv++) vfnorm += vf[kv] * vf[kv];
    pdist[kp] = sqrt(vfnorm);
  }
  histinfo_fixbar(pdist, numinbar, 0., rmax);
  if (flag_case == 0) {
    save_macro<int>(numinbar, "pdist");
  } else if (flag_case == 1) {
    save_macro<int>(numinbar, "pdist_before");
  } else if (flag_case == 2) {
    save_macro<int>(numinbar, "pdist_after");
  }

  int Nn = S_x.size('n');
  vector<double> ndist(Nn);
  auto &Sn = S_x.list('n');
  for (int kn = 0; kn < Nn; kn++) {
    auto &vf = Sn[kn].velocity();
    double vfnorm = 0.;
    for (int kv = 0; kv < 3; kv++) vfnorm += vf[kv] * vf[kv];
    ndist[kn] = sqrt(vfnorm);
  }
  histinfo_fixbar(ndist, numinbar, 0., rmax);
  if (flag_case == 0) {
    save_macro<int>(numinbar, "ndist");
  } else if (flag_case == 1) {
    save_macro<int>(numinbar, "ndist_before");
  } else if (flag_case == 2) {
    save_macro<int>(numinbar, "ndist_after");
  }

  int Nf = S_x.size('f');
  vector<double> fdist(Nf);
  auto &Sf = S_x.list('f');
  for (int kf = 0; kf < Nf; kf++) {
    auto &vf = Sf[kf].velocity();
    double vfnorm = 0.;
    for (int kv = 0; kv < 3; kv++) vfnorm += vf[kv] * vf[kv];
    fdist[kf] = sqrt(vfnorm);
  }
  histinfo_fixbar(fdist, numinbar, 0., rmax);
  if (flag_case == 0) {
    save_macro<int>(numinbar, "fdist");
  } else if (flag_case == 1) {
    save_macro<int>(numinbar, "fdist_before");
  } else if (flag_case == 2) {
    save_macro<int>(numinbar, "fdist_after");
  }
}

void save_particles(const NeParticleGroup &S_x_before,
                    const NeParticleGroup &S_x_after) {
  int Np0 = S_x_before.size('p');
  int Nn0 = S_x_before.size('n');
  int Np1 = S_x_after.size('p');
  int Nn1 = S_x_after.size('n');

  if ((Np1 + Nn1) > (Np0 + Nn0)) {
    ofstream file0;

    {
      file0.open("result/ParticleP0.txt");
      auto &Sp = S_x_before.list('p');
      for (int kp = 0; kp < Np0; kp++) {
        const auto vp = Sp[kp].velocity();
        file0 << vp[0] << ' ' << vp[1] << ' ' << vp[2] << endl;
      }
      file0.close();
    }
    {
      file0.open("result/ParticleN0.txt");
      auto &Sp = S_x_before.list('n');
      for (int kp = 0; kp < Nn0; kp++) {
        auto vp = Sp[kp].velocity();
        file0 << vp[0] << ' ' << vp[1] << ' ' << vp[2] << endl;
      }
      file0.close();
    }

    {
      file0.open("result/ParticleP1.txt");
      auto &Sp = S_x_after.list('p');
      for (int kp = 0; kp < Np1; kp++) {
        const auto vp = Sp[kp].velocity();
        file0 << vp[0] << ' ' << vp[1] << ' ' << vp[2] << endl;
      }
      file0.close();
    }

    {
      file0.open("result/ParticleN1.txt");
      auto &Sp = S_x_after.list('n');
      for (int kp = 0; kp < Nn1; kp++) {
        const auto vp = Sp[kp].velocity();
        file0 << vp[0] << ' ' << vp[1] << ' ' << vp[2] << endl;
      }
      file0.close();
    }
    exit(0);
  }
}

void saveparameter(const ParaClass &para, const NumericGridClass &grid) {
  ofstream file0;

  file0.open("result/parameter.txt");
  file0 << setprecision(15);
  file0 << "coeff_binarycoll " << para.coeff_binarycoll << endl
        << "resample_ratio " << para.resample_ratio << endl
        << "Npickup_neg " << para.Npickup_neg << endl
        << "Nfreq " << para.Nfreq << endl
        << "xmax " << grid.xmax << endl
        << "xmin " << grid.xmin << endl
        << "vmax " << grid.vmax << endl
        << "vmin " << grid.vmin << endl
        << "tmax " << grid.tmax << endl
        << "Nx " << grid.Nx << endl
        << "Nt " << grid.Nt << endl
        << "Nv " << grid.Nv << endl
        << "dx " << grid.dx << endl
        << "dv " << grid.dv << endl
        << "dt " << grid.dt << endl
        << "Neff " << grid.Neff << endl
        << "Neff_F " << grid.Neff_F << endl
        << "collision_kernel " << para.flag_collision << endl
        << "lambda_Poisson " << para.lambda_Poisson << endl
        << "resample_spatial_ratio " << para.resample_spatial_ratio << endl
        << "sync_time_interval " << para.sync_time_interval << endl
        << "resample_sync_ratio " << para.resample_sync_ratio << endl;
  file0.close();

  file0.open("result/parameter2.txt");
  file0 << setprecision(15);
  file0 << "method_binarycoll " << para.method_binarycoll << endl
        << "bdry_x " << grid.bdry_x << endl
        << "bdry_v " << grid.bdry_v << endl
        << "method " << para.method << endl;
  file0.close();
}

void save_initial(IniValClass &inidata) {
  ofstream file0;

  file0.open("result/parameter.txt", ios::app);
  file0 << setprecision(15);
  if (inidata.probname == "LandauDamping") {
    file0 << "LD_alpha " << inidata.LD_alpha << endl
          << "totalmass " << inidata.totalmass << endl;
  } else if (inidata.probname == "BumpOnTail") {
    file0 << "BOT_beta " << inidata.BOT_beta << endl
          << "BOT_rho0 " << inidata.BOT_rho0 << endl
          << "BOT_Tprt " << inidata.BOT_Tprt << endl
          << "BOT_dTprt " << inidata.BOT_dTprt << endl
          << "BOT_Tx " << inidata.BOT_Tx << endl
          << "BOT_ub " << inidata.BOT_ub << endl;
  }
  file0.close();

  file0.open("result/parameter2.txt", ios::app);
  file0 << "problem_name " << inidata.probname << endl;
  file0.close();
}

void save_macro_evolution(std::vector<NeParticleGroup> &S_x,
                          const NumericGridClass &grid) {
  for (int kx = 0; kx < grid.Nx; kx++) S_x[kx].computemoments();

  update_macro(S_x, grid);

  //	save_m012_PN(ptr_S_x, grid);
  save_rhouT(&S_x[0], grid);
  save_rhouT_F(&S_x[0], grid);
  save_E(&S_x[0], grid);
  save_dist(&S_x[0], grid);
  // save_NpNn(ptr_S_x,  grid);
  save_rhouT_maxwellian(&S_x[0], grid);

  ofstream file0;

  file0.open("result/numOfDist.txt");
  file0 << K_SAVE_TIME;
  file0.close();
}
}  // namespace coulomb