#include "Classes.h"

namespace coulomb {
void save_homo_dist(NeParticleGroup *S_x, int Nr, int flag_case);

/* ======================================================== *\
        Use Fourier transform for 3D interpolation
\* ======================================================== */

// xyz_minmax = [xmin, xmax, ymin, ymax, zmin, zmax];

void interp3d_xyzminmax(NeParticleGroup *S_x, double *xyz_minmax) {
  for (int k = 0; k < 6; k++) {
    xyz_minmax[k] = 0.;
  }

  int Np = S_x->size('p');
  int Nn = S_x->size('n');
  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');

  for (int kp = 0; kp < Np; kp++) {
    auto &v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      xyz_minmax[2 * k2] = min(xyz_minmax[2 * k2], v0[k2]);
      xyz_minmax[2 * k2 + 1] = max(xyz_minmax[2 * k2 + 1], v0[k2]);
    }
  }
  for (int kn = 0; kn < Nn; kn++) {
    auto &v0 = Sn[kn].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      xyz_minmax[2 * k2] = min(xyz_minmax[2 * k2], v0[k2]);
      xyz_minmax[2 * k2 + 1] = max(xyz_minmax[2 * k2 + 1], v0[k2]);
    }
  }

  for (int k2 = 0; k2 < 3; k2++) {
    xyz_minmax[2 * k2] -= 1e-6;
    xyz_minmax[2 * k2 + 1] += 1e+6;
  }
}

void interp3d_renormalize(NeParticleGroup *S_x, NeParticleGroup *S_x_new) {
  int Np = S_x->size('p');
  int Nn = S_x->size('n');

  Particle1d3d S_one;

  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');

  double *xyz_minmax = S_x->xyz_minmax;

  // interp3d_xyzminmax(S_x, xyz_minmax);
  double Lxyz[3];
  for (int k2 = 0; k2 < 3; k2++) {
    Lxyz[k2] = xyz_minmax[2 * k2 + 1] - xyz_minmax[2 * k2];
  }
  // for (int k2 = 0; k2 < 6; k2 ++)
  // // cout << xyz_minmax[k2] << ' ';
  // // cout << endl;

  // renormalizaed value
  double v1[3];
  for (int kp = 0; kp < Np; kp++) {
    auto &v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v1[k2] = (v0[k2] - xyz_minmax[2 * k2]) * 2.0 * pi / Lxyz[k2];
    }
    S_one.set_velocity(v1);
    S_x_new->push_back(S_one, 'p');
  }

  for (int kn = 0; kn < Nn; kn++) {
    auto &v0 = Sn[kn].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v1[k2] = (v0[k2] - xyz_minmax[2 * k2]) * 2.0 * pi / Lxyz[k2];
    }
    S_one.set_velocity(v1);
    S_x_new->push_back(S_one, 'n');
  }

  // renormalize Maxwellian
  S_x_new->u1M = (S_x->u1M - xyz_minmax[0]) * 2.0 * pi / Lxyz[0];
  S_x_new->u2M = (S_x->u2M - xyz_minmax[2]) * 2.0 * pi / Lxyz[1];
  S_x_new->u3M = (S_x->u3M - xyz_minmax[4]) * 2.0 * pi / Lxyz[2];
  S_x_new->T1M = S_x->TprtM * (4.0 * pi * pi / Lxyz[0] / Lxyz[0]);
  S_x_new->T2M = S_x->TprtM * (4.0 * pi * pi / Lxyz[1] / Lxyz[1]);
  S_x_new->T3M = S_x->TprtM * (4.0 * pi * pi / Lxyz[2] / Lxyz[2]);
}

void interp3d_rescale(Particle1d3d *Sp, int Np, double *xyz_minmax) {
  // rescale to the original coordinates stored in xyz_minmax
  double Lxyz[3];
  for (int k2 = 0; k2 < 3; k2++) {
    Lxyz[k2] = xyz_minmax[2 * k2 + 1] - xyz_minmax[2 * k2];
  }

  for (int kp = 0; kp < Np; kp++) {
    auto v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v0[k2] = xyz_minmax[2 * k2] + v0[k2] * Lxyz[k2] / (2.0 * pi);
    }
    Sp[kp].set_velocity(v0);
  }
}

/**
 the frequency grids
*/
void interp_ifreq(complex<double> *ifreq, int Nfreq) {
  for (int j = 0; j < Nfreq / 2 + 1; j++) {
    ifreq[j] = complex<double>(0., (double)j);
  }
  for (int j = Nfreq / 2 + 1; j < Nfreq; j++) {
    ifreq[j] = complex<double>(0., (double)(j - Nfreq));
  }
}

int freq_sequence(int kth, int Nfreq) {
  int kfreq = kth;
  if (kth >= Nfreq / 2 + 1) kfreq = kth - Nfreq;
  return kfreq;
}

int freq_sequence_reverse(int kfreq, int Nfreq) {
  int kth = kfreq;
  if (kfreq < 0) kth = kfreq + Nfreq;
  return kth;
}

void interp_freq(vector<double> &freq, int Nfreq) {
  for (int j = 0; j < Nfreq; j++) freq[j] = (double)(freq_sequence(j, Nfreq));
}

void interp_freq_aug(vector<int> &loc, int Nfreq, int augFactor) {
  for (int j = 0; j < Nfreq; j++) {
    int kfreq = freq_sequence(j, Nfreq);
    loc[j] = freq_sequence_reverse(kfreq, augFactor * Nfreq);
  }
}

// compute the Fourier coefficients in  3D

void interp3d_fft(NeParticleGroup *S_x, complex<double> *Fouriercoeff,
                  int Nfreq1, int Nfreq2, int Nfreq3) {
  complex<double> *ifreq1 = new complex<double>[Nfreq1];  // 1i *freq
  complex<double> *ifreq2 = new complex<double>[Nfreq2];  // 1i *freq
  complex<double> *ifreq3 = new complex<double>[Nfreq3];  // 1i *freq

  int Np = S_x->size('p');
  int Nn = S_x->size('n');

  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');

  double Neff = 1.0;

  // double Neff_temp = 1./Np;

  interp_ifreq(ifreq1, Nfreq1);
  interp_ifreq(ifreq2, Nfreq2);
  interp_ifreq(ifreq3, Nfreq3);

  double cubic_2pi = 8.0 * pi * pi * pi;

  double coeff_fft = 1. / cubic_2pi * Nfreq1 * Nfreq2 * Nfreq3;
  double maxFS = 0.0;

  // the (i,j,k)-th element of the array with size (Nx,Ny,Nz), you would use the
  // expression an_array[k + Nz * (j + Ny * i)].

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);
        Fouriercoeff[kk] = complex<double>(0., 0.);
        for (int kp = 0; kp < Np; kp++) {
          auto &vp = Sp[kp].velocity();
          complex<double> expterm =
              vp[0] * ifreq1[kk1] + vp[1] * ifreq2[kk2] + vp[2] * ifreq3[kk3];
          Fouriercoeff[kk] += exp(-expterm);
        }
        for (int kn = 0; kn < Nn; kn++) {
          auto &vp = Sn[kn].velocity();
          complex<double> expterm =
              vp[0] * ifreq1[kk1] + vp[1] * ifreq2[kk2] + vp[2] * ifreq3[kk3];
          Fouriercoeff[kk] -= exp(-expterm);
        }
        // Fouriercoeff[kk] *= Neff * coeff_fft;
        Fouriercoeff[kk] *= Neff * coeff_fft / (Nfreq1 * Nfreq2 * Nfreq3);
        maxFS = max(maxFS, abs(Fouriercoeff[kk]));
      }
    }
  }

  delete[] ifreq1;
  delete[] ifreq2;
  delete[] ifreq3;
  // delete [] Sp;
  // delete [] Sn;
}

void interp3d_fft_approx_terms(complex<double> *Fouriercoeff, vector<double> &f,
                               int Nfreq1, int Nfreq2, int Nfreq3,
                               int augFactor, int orderx, int ordery,
                               int orderz) {
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;

  vector<double> freq1(Nfreq1);  // 1i *freq
  vector<double> freq2(Nfreq2);  // 1i *freq
  vector<double> freq3(Nfreq3);  // 1i *freq

  interp_freq(freq1, Nfreq1);
  interp_freq(freq2, Nfreq2);
  interp_freq(freq3, Nfreq3);

  vector<int> loc1(Nfreq1);  // 1i *freq
  vector<int> loc2(Nfreq2);  // 1i *freq
  vector<int> loc3(Nfreq3);  // 1i *freq

  interp_freq_aug(loc1, Nfreq1, augFactor);
  interp_freq_aug(loc2, Nfreq2, augFactor);
  interp_freq_aug(loc3, Nfreq3, augFactor);

  fftw_complex *FSaug, a_fftw_complex_variable, *fin;
  fin = (fftw_complex *)fftw_malloc(sizeF * sizeof(a_fftw_complex_variable));
  FSaug = (fftw_complex *)fftw_malloc(sizeF * sizeof(a_fftw_complex_variable));
  fftw_plan plan3d_ft;

  plan3d_ft = fftw_plan_dft_3d(augFactor * Nfreq1, augFactor * Nfreq2,
                               augFactor * Nfreq3, fin, FSaug, FFTW_FORWARD,
                               FFTW_ESTIMATE);

  for (int kfc = 0; kfc < sizeF; kfc++) {
    fin[kfc][0] = f[kfc];
    fin[kfc][1] = 0.;
  }

  fftw_execute(plan3d_ft);

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);

        int kk1aug = loc1[kk1];
        int kk2aug = loc2[kk2];
        int kk3aug = loc3[kk3];
        int kkaug = kk3aug +
                    augFactor * Nfreq3 * (kk2aug + augFactor * Nfreq2 * kk1aug);

        double freq = 1.0;
        for (int kx = 0; kx < orderx; kx++) freq *= freq1[kk1];
        for (int kx = 0; kx < ordery; kx++) freq *= freq2[kk2];
        for (int kx = 0; kx < orderz; kx++) freq *= freq3[kk3];

        if ((orderx + ordery + orderz) == 2) {
          if ((orderx == 2) || (ordery == 2) || (orderz == 2))
            freq *= -.5;
          else
            freq *= -1;
        }

        if ((orderx + ordery + orderz) == 1) {
          Fouriercoeff[kk] +=
              freq * complex<double>(FSaug[kkaug][1], -FSaug[kkaug][0]);
        } else {
          Fouriercoeff[kk] +=
              freq * complex<double>(FSaug[kkaug][0], FSaug[kkaug][1]);
        }
      }
    }
  }

  fftw_destroy_plan(plan3d_ft);
  fftw_free(fin);
  fftw_free(FSaug);
}

void interp3d_fft_approx(NeParticleGroup *S_x, complex<double> *Fouriercoeff,
                         int Nfreq1, int Nfreq2, int Nfreq3) {
  int augFactor = 2;

  int Np = S_x->size('p');
  int Nn = S_x->size('n');

  auto &Sp = S_x->list('p');
  auto &Sn = S_x->list('n');

  double cubic_2pi = 8.0 * pi * pi * pi;

  double coeff_fft = 1. / cubic_2pi;

  // for the (i,j,k)-th element of the array with size (Nx,Ny,Nz), use the
  // expression an_array[k + Nz * (j + Ny * i)].

  // create f, fx, fy, fz, fxx, fyy, fzz, fxy ...
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;
  double dx = 2.0 * pi / augFactor / Nfreq1;
  double dy = 2.0 * pi / augFactor / Nfreq2;
  double dz = 2.0 * pi / augFactor / Nfreq3;

  vector<double> f(sizeF);
  vector<double> fx(sizeF);
  vector<double> fy(sizeF);
  vector<double> fz(sizeF);
  vector<double> fxx(sizeF);
  vector<double> fyy(sizeF);
  vector<double> fzz(sizeF);
  vector<double> fxy(sizeF);
  vector<double> fxz(sizeF);
  vector<double> fyz(sizeF);

  // cout << "Approx 1" << endl;

  for (int kk = 0; kk < sizeF; kk++) {
    f[kk] = 0.;
    fx[kk] = 0.;
    fy[kk] = 0.;
    fz[kk] = 0.;
    fxx[kk] = 0.;
    fyy[kk] = 0.;
    fzz[kk] = 0.;
    fxy[kk] = 0.;
    fxz[kk] = 0.;
    fyz[kk] = 0.;
  }

  for (int kp = 0; kp < Np; kp++) {
    double x0 = Sp[kp].velocity(0);
    double y0 = Sp[kp].velocity(1);
    double z0 = Sp[kp].velocity(2);
    int xloc = (int)(floor(x0 / dx + 0.5));
    int yloc = (int)(floor(y0 / dy + 0.5));
    int zloc = (int)(floor(z0 / dz + 0.5));
    if (xloc >= augFactor * Nfreq1) xloc--;
    if (yloc >= augFactor * Nfreq2) yloc--;
    if (zloc >= augFactor * Nfreq3) zloc--;
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    int loc = zloc + augFactor * Nfreq3 * (yloc + augFactor * Nfreq2 * xloc);

    if ((loc >= sizeF) || (loc < 0))
      cout << x0 << ' ' << y0 << ' ' << z0 << ' ' << dx << ' ' << loc << endl;

    f[loc]++;
    fx[loc] += xdelta;
    fy[loc] += ydelta;
    fz[loc] += zdelta;
    fxx[loc] += xdelta * xdelta;
    fyy[loc] += ydelta * ydelta;
    fzz[loc] += zdelta * zdelta;
    fxy[loc] += xdelta * ydelta;
    fyz[loc] += ydelta * zdelta;
    fxz[loc] += zdelta * xdelta;
  }

  // cout << "Approx 2" << endl;

  for (int kp = 0; kp < Nn; kp++) {
    double x0 = Sn[kp].velocity(0);
    double y0 = Sn[kp].velocity(1);
    double z0 = Sn[kp].velocity(2);
    // int xloc = floor(x0/dx);
    // int yloc = floor(y0/dy);
    // int zloc = floor(z0/dz);
    int xloc = (int)(floor(x0 / dx + 0.5));
    int yloc = (int)(floor(y0 / dy + 0.5));
    int zloc = (int)(floor(z0 / dz + 0.5));
    if (xloc >= augFactor * Nfreq1) xloc--;
    if (yloc >= augFactor * Nfreq2) yloc--;
    if (zloc >= augFactor * Nfreq3) zloc--;
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    int loc = zloc + augFactor * Nfreq3 * (yloc + augFactor * Nfreq2 * xloc);

    if ((loc >= sizeF) || (loc < 0)) {
      cout << "error: in approximation. Particle moved out of range. kx = "
           << xloc << ' ' << yloc << ' ' << zloc << ' ' << loc << endl;
      exit(0);
    }

    f[loc]--;
    fx[loc] -= xdelta;
    fy[loc] -= ydelta;
    fz[loc] -= zdelta;
    fxx[loc] -= xdelta * xdelta;
    fyy[loc] -= ydelta * ydelta;
    fzz[loc] -= zdelta * zdelta;
    fxy[loc] -= xdelta * ydelta;
    fyz[loc] -= ydelta * zdelta;
    fxz[loc] -= zdelta * xdelta;
  }
  // cout << "Approx 3" << endl;

  for (int kk = 0; kk < Nfreq1 * Nfreq2 * Nfreq3; kk++)
    Fouriercoeff[kk] = complex<double>(0., 0.);

  interp3d_fft_approx_terms(Fouriercoeff, f, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            0, 0, 0);

  interp3d_fft_approx_terms(Fouriercoeff, fx, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            1, 0, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fy, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            0, 1, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fz, Nfreq1, Nfreq2, Nfreq3, augFactor,
                            0, 0, 1);
  interp3d_fft_approx_terms(Fouriercoeff, fxx, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 2, 0, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fyy, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 0, 2, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fzz, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 0, 0, 2);
  interp3d_fft_approx_terms(Fouriercoeff, fxy, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 1, 1, 0);
  interp3d_fft_approx_terms(Fouriercoeff, fxz, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 1, 0, 1);
  interp3d_fft_approx_terms(Fouriercoeff, fyz, Nfreq1, Nfreq2, Nfreq3,
                            augFactor, 0, 1, 1);

  for (int kk = 0; kk < Nfreq1 * Nfreq2 * Nfreq3; kk++)
    Fouriercoeff[kk] *= coeff_fft;

  // cout << "Approx finished." << endl;
}

void filter_Fourier(complex<double> *Fouriercoeff,
                    vector<int> &flag_Fouriercoeff, int size_FC) {
  // double thres = 10.0;
  for (int k = 0; k < size_FC; k++) {
    flag_Fouriercoeff[k] = 1;
    /*
    double abs_FC = abs(Fouriercoeff[k]);
    if (abs_FC < thres) {
      Fouriercoeff[k] *= 0.;
      flag_Fouriercoeff[k] = 0;
    }
    */
  }
}

double interp3d_fvalue(double *Sf, complex<double> *Fouriercoeff,
                       complex<double> *ifreq1, complex<double> *ifreq2,
                       complex<double> *ifreq3, vector<int> &flag_Fouriercoeff,
                       int Nfreq1, int Nfreq2, int Nfreq3) {
  complex<double> fval_c(0., 0.);

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kkk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);
        if (flag_Fouriercoeff[kkk] > 0) {
          fval_c += Fouriercoeff[kkk] *
                    exp(ifreq1[kk1] * Sf[0] + ifreq2[kk2] * Sf[1] +
                        ifreq3[kk3] * Sf[2]);
        }
      }
    }
  }

  // return real(fval_c)/(Nfreq1*Nfreq2*Nfreq3);
  return real(fval_c);
}

double interp3d_fvalue_approx(double deltax, double deltay, double deltaz,
                              double f, double fx, double fy, double fz,
                              double fxx, double fyy, double fzz, double fxy,
                              double fxz, double fyz) {
  double f0;
  f0 = f + fx * deltax + fy * deltay + fx * deltay +
       .5 * fxx * deltax * deltax + .5 * fyy * deltay * deltay +
       .5 * fzz * deltaz * deltaz + fxy * deltax * deltay +
       fxz * deltax * deltaz + fyz * deltay * deltaz;
  return f0;
}

void interp3d_acceptsampled(double *Sf, NeParticleGroup *ptr_S_x_incell,
                            double fval, double &maxf) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new

    double keeprate = maxf / (1.5 * abs(fval));

    maxf = 1.5 * abs(fval);

    int Np_remove = myfloor((1 - keeprate) * ptr_S_x_incell->size('p'));
    int Nn_remove = myfloor((1 - keeprate) * ptr_S_x_incell->size('n'));

    for (int kp = 0; kp < Np_remove; kp++) {
      int k_remove = (int)(myrand() * ptr_S_x_incell->size('p'));
      ptr_S_x_incell->erase(k_remove, 'p');
    }

    for (int kn = 0; kn < Nn_remove; kn++) {
      int k_remove = (int)(myrand() * ptr_S_x_incell->size('n'));
      ptr_S_x_incell->erase(k_remove, 'n');
    }
  }

  // accept this particle with rate abs(fval/maxf)
  if (myrand() < (abs(fval / maxf))) {
    double sum_Sf_pi_sq = 0.;
    for (int kv = 0; kv < 3; kv++)
      sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
    if (sqrt(sum_Sf_pi_sq) < pi) {
      Particle1d3d S_one({Sf[0], Sf[1], Sf[2]});
      if (fval > 0) {
        ptr_S_x_incell->push_back(S_one, 'p');
      } else {
        ptr_S_x_incell->push_back(S_one, 'n');
      }
    }
  }
}

void resampleF_acceptsampled(double *Sf, NeParticleGroup *ptr_S_x_incell,
                             double fval, double &maxf) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new

    double keeprate = maxf / (1.5 * abs(fval));

    maxf = 1.5 * abs(fval);

    int Np_remove = myfloor((1 - keeprate) * ptr_S_x_incell->size('f'));

    for (int kp = 0; kp < Np_remove; kp++) {
      int k_remove = (int)(myrand() * ptr_S_x_incell->size('f'));
      ptr_S_x_incell->erase(k_remove, 'f');
    }
  }

  // accept this particle with rate abs(fval/maxf)
  if (myrand() < (abs(fval / maxf))) {
    double sum_Sf_pi_sq = 0.;
    for (int kv = 0; kv < 3; kv++)
      sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
    if (sqrt(sum_Sf_pi_sq) < pi) {
      Particle1d3d S_one({Sf[0], Sf[1], Sf[2]});
      ptr_S_x_incell->push_back(S_one, 'f');
    }
  }
}
/******************************************************************/
/* ---------- Use Fourier transform for 3D interpolation -------- */
/******************************************************************/

/*
  Find the coarse approximation with the given Fourier coefficients
  Need to include 'fftw3.f'
*/
void interp3d_fcoarse(complex<double> *Fouriercoeff, vector<double> &fcoarse,
                      int Nfreq1, int Nfreq2, int Nfreq3) {
  // double Lcubic = (double) (Nfreq1*Nfreq2*Nfreq3);

  fftw_complex *FC, a_fftw_complex_variable, *fcoarse_c;
  FC = (fftw_complex *)fftw_malloc(Nfreq1 * Nfreq2 * Nfreq3 *
                                   sizeof(a_fftw_complex_variable));
  fcoarse_c = (fftw_complex *)fftw_malloc(Nfreq1 * Nfreq2 * Nfreq3 *
                                          sizeof(a_fftw_complex_variable));
  fftw_plan plan3d_ift;

  for (int kfc = 0; kfc < Nfreq1 * Nfreq2 * Nfreq3; kfc++) {
    FC[kfc][0] = real(Fouriercoeff[kfc]);
    FC[kfc][1] = imag(Fouriercoeff[kfc]);
  }

  // use fftw to obtain an estimation of f_p - f_n
  // // cout << " resample 1.3" << endl;
  plan3d_ift = fftw_plan_dft_3d(Nfreq1, Nfreq2, Nfreq3, FC, fcoarse_c,
                                FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(plan3d_ift);

  for (int kfc = 0; kfc < Nfreq1 * Nfreq2 * Nfreq3; kfc++) {
    // fcoarse[kfc] = fcoarse_c[kfc][0] / Lcubic;
    fcoarse[kfc] = fcoarse_c[kfc][0];
  }

  fftw_destroy_plan(plan3d_ift);
  fftw_free(FC);
  fftw_free(fcoarse_c);
}

/** The derivatives of f on grids
 */
// CONTINUE HERE

void interp3d_fxyz_terms(complex<double> *Fouriercoeff, vector<double> &f,
                         int Nfreq1, int Nfreq2, int Nfreq3, int augFactor,
                         int orderx, int ordery, int orderz) {
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;
  complex<double> *FSaug = new complex<double>[sizeF];
  for (int k = 0; k < sizeF; k++) FSaug[k] = complex<double>(0., 0.);

  vector<double> freq1(Nfreq1);  // 1i *freq
  vector<double> freq2(Nfreq2);  // 1i *freq
  vector<double> freq3(Nfreq3);  // 1i *freq

  interp_freq(freq1, Nfreq1);
  interp_freq(freq2, Nfreq2);
  interp_freq(freq3, Nfreq3);

  vector<int> loc1(Nfreq1);  // 1i *freq
  vector<int> loc2(Nfreq2);  // 1i *freq
  vector<int> loc3(Nfreq3);  // 1i *freq

  interp_freq_aug(loc1, Nfreq1, augFactor);
  interp_freq_aug(loc2, Nfreq2, augFactor);
  interp_freq_aug(loc3, Nfreq3, augFactor);

  for (int kk1 = 0; kk1 < Nfreq1; kk1++) {
    for (int kk2 = 0; kk2 < Nfreq2; kk2++) {
      for (int kk3 = 0; kk3 < Nfreq3; kk3++) {
        int kk = kk3 + Nfreq3 * (kk2 + Nfreq2 * kk1);

        int kk1aug = loc1[kk1];
        int kk2aug = loc2[kk2];
        int kk3aug = loc3[kk3];
        int kkaug = kk3aug +
                    augFactor * Nfreq3 * (kk2aug + augFactor * Nfreq2 * kk1aug);

        double freq = 1.0;
        for (int kx = 0; kx < orderx; kx++) freq *= freq1[kk1];
        for (int kx = 0; kx < ordery; kx++) freq *= freq2[kk2];
        for (int kx = 0; kx < orderz; kx++) freq *= freq3[kk3];

        FSaug[kkaug] = freq * Fouriercoeff[kk];
      }
    }
  }

  interp3d_fcoarse(FSaug, f, augFactor * Nfreq1, augFactor * Nfreq2,
                   augFactor * Nfreq3);

  delete[] FSaug;
}

void interp3d_fxyz(complex<double> *Fouriercoeff, vector<double> &f,
                   vector<double> &fx, vector<double> &fy, vector<double> &fz,
                   vector<double> &fxx, vector<double> &fyy,
                   vector<double> &fzz, vector<double> &fxy,
                   vector<double> &fxz, vector<double> &fyz, int Nfreq1,
                   int Nfreq2, int Nfreq3, int augFactor) {
  interp3d_fxyz_terms(Fouriercoeff, f, Nfreq1, Nfreq2, Nfreq3, augFactor, 0, 0,
                      0);
  interp3d_fxyz_terms(Fouriercoeff, fx, Nfreq1, Nfreq2, Nfreq3, augFactor, 1, 0,
                      0);
  interp3d_fxyz_terms(Fouriercoeff, fy, Nfreq1, Nfreq2, Nfreq3, augFactor, 0, 1,
                      0);
  interp3d_fxyz_terms(Fouriercoeff, fz, Nfreq1, Nfreq2, Nfreq3, augFactor, 0, 0,
                      1);
  interp3d_fxyz_terms(Fouriercoeff, fxx, Nfreq1, Nfreq2, Nfreq3, augFactor, 2,
                      0, 0);
  interp3d_fxyz_terms(Fouriercoeff, fyy, Nfreq1, Nfreq2, Nfreq3, augFactor, 0,
                      2, 0);
  interp3d_fxyz_terms(Fouriercoeff, fzz, Nfreq1, Nfreq2, Nfreq3, augFactor, 0,
                      0, 2);
  interp3d_fxyz_terms(Fouriercoeff, fxy, Nfreq1, Nfreq2, Nfreq3, augFactor, 1,
                      1, 0);
  interp3d_fxyz_terms(Fouriercoeff, fxz, Nfreq1, Nfreq2, Nfreq3, augFactor, 1,
                      0, 1);
  interp3d_fxyz_terms(Fouriercoeff, fyz, Nfreq1, Nfreq2, Nfreq3, augFactor, 0,
                      1, 1);
}

void func_fourierupper3d(int N, vector<double> &fc, vector<double> &f_up);
void merge_NeParticleGroup(NeParticleGroup &S_x,
                           const NeParticleGroup &S_x_new);
void mergeF_NeParticleGroup(NeParticleGroup &S_x,
                            const NeParticleGroup &S_x_new);
// void interp3d_fft_eachlevel(NeParticleGroup * S_x, MultlLevelGroup * MLsol,
// int Nlevel); void interp3d_fft_ml(complex<double> *Fouriercoeff, int *
// flag_Fouriercoeff, MultlLevelGroup * MLsol, int Nlevel);
void save_complex(int Nx, complex<double> *FS, string filename);

/*
  Sampled particles are stored in particles_Sp_sampled and particles_Sn_sampled
  with size particles%Np_sampled and particles%Nn_sampled
*/

// void samplefromfourier3d(NeParticleGroup * S_x, NeParticleGroup * S_x_new,
// MultlLevelGroup * MLsol, int Nlevel) {
NeParticleGroup samplefromfourier3d(NeParticleGroup &S_x, int Nfreq) {
  NeParticleGroup S_x_new;

  double Neff = 1.0;
  bool flag_useApproximation = true;

  /* Normalize particle velocity to [0 2*pi] */
  S_x.computexyzrange();

  NeParticleGroup S_x_renormalized;
  NeParticleGroup *ptr_S_x_renormalized = &S_x_renormalized;

  interp3d_renormalize(&S_x, ptr_S_x_renormalized);

  /* Prepare the grids in physical space and frequence space */
  // double dx = 2.0*pi/Nfreq;

  complex<double> *ifreq = new complex<double>[Nfreq];  // 1i *freq
  interp_ifreq(ifreq, Nfreq);
  vector<double> interp_x(Nfreq);
  for (int kx = 0; kx < Nfreq; kx++) interp_x[kx] = kx * 2 * pi / Nfreq;

  /* Compute the Fourier coefficient */
  complex<double> *Fouriercoeff = new complex<double>[Nfreq * Nfreq * Nfreq];
  vector<int> flag_Fouriercoeff(Nfreq * Nfreq * Nfreq);

  if (flag_useApproximation)
    interp3d_fft_approx(ptr_S_x_renormalized, Fouriercoeff, Nfreq, Nfreq,
                        Nfreq);
  else
    interp3d_fft(ptr_S_x_renormalized, Fouriercoeff, Nfreq, Nfreq, Nfreq);

  filter_Fourier(
      Fouriercoeff, flag_Fouriercoeff,
      Nfreq * Nfreq * Nfreq);  // Apply the filter on Fourier coefficients

  // cout << " F coeff computed " << endl;

  /* Compute a coarse interpolation in physical space */
  vector<double> fcoarse(Nfreq * Nfreq * Nfreq);
  interp3d_fcoarse(Fouriercoeff, fcoarse, Nfreq, Nfreq, Nfreq);

  int augFactor = 2;
  int sizeF = augFactor * augFactor * augFactor * Nfreq * Nfreq * Nfreq;

  vector<double> f(sizeF);
  vector<double> fx(sizeF);
  vector<double> fy(sizeF);
  vector<double> fz(sizeF);
  vector<double> fxx(sizeF);
  vector<double> fyy(sizeF);
  vector<double> fzz(sizeF);
  vector<double> fxy(sizeF);
  vector<double> fxz(sizeF);
  vector<double> fyz(sizeF);

  interp3d_fxyz(Fouriercoeff, f, fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz,
                Nfreq, Nfreq, Nfreq, augFactor);

  /* evaluate the upperbound of f */
  vector<double> f_up(sizeF);
  func_fourierupper3d(augFactor * Nfreq, f, f_up);

  /* refined x grid */
  double dxaug = 2.0 * pi / Nfreq / augFactor;
  vector<double> interp_xaug(Nfreq * augFactor);
  for (int kx = 0; kx < Nfreq * augFactor; kx++)
    interp_xaug[kx] = kx * 2 * pi / Nfreq / augFactor;

  /* create a NeParticleGroup to host the P and N particles in current cell */

  /* Start sampling */

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        int kk = kz + augFactor * Nfreq * (ky + augFactor * Nfreq * kx);

        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        double fcc = f_up[kk];

        if (fcc < abs(f[kk])) cout << "ERROR: small bound!" << endl;

        double maxf = 1.5 * abs(fcc);
        int N_incell = myfloor(maxf * dxaug * dxaug * dxaug / Neff);

        int k_virtual = 0;
        NeParticleGroup S_x_incell;
        NeParticleGroup *ptr_S_x_incell = &S_x_incell;

        while (k_virtual < N_incell) {
          // create a particle in the cell
          // double Sf[3] = {xc+myrand()*dx, yc+myrand()*dx, zc+myrand()*dx};
          double deltax = myrand() * dxaug - 0.5 * dxaug;
          double deltay = myrand() * dxaug - 0.5 * dxaug;
          double deltaz = myrand() * dxaug - 0.5 * dxaug;
          double Sf[3] = {xc + deltax, yc + deltay, zc + deltaz};

          // compute f at this point
          double fval = 0;
          if (flag_useApproximation)
            fval = interp3d_fvalue_approx(deltax, deltay, deltaz, f[kk], fx[kk],
                                          fy[kk], fz[kk], fxx[kk], fyy[kk],
                                          fzz[kk], fxy[kk], fxz[kk], fyz[kk]);
          else
            fval = interp3d_fvalue(Sf, Fouriercoeff, ifreq, ifreq, ifreq,
                                   flag_Fouriercoeff, Nfreq, Nfreq, Nfreq);

          // reset current cell if fval>maxf, otherwise continue sampling in
          // current cell
          interp3d_acceptsampled(Sf, ptr_S_x_incell, fval, maxf);

          // reset N_incell if maxf is changed
          N_incell = myfloor(maxf / (Neff / (dxaug * dxaug * dxaug)));
          k_virtual++;
        }

        merge_NeParticleGroup(S_x_new, S_x_incell);
      }
    }
  }

  // cout << "Resampled." << endl;

  // rescale to the original coordinates
  auto &Sp_sampled = S_x_new.list('p');
  auto &Sn_sampled = S_x_new.list('n');
  double *xyz_minmax = S_x.xyz_minmax;
  interp3d_rescale(&Sp_sampled[0], S_x_new.size('p'), xyz_minmax);
  interp3d_rescale(&Sn_sampled[0], S_x_new.size('n'), xyz_minmax);

  // cout << "Rescaled." << endl;

  // free memory
  delete[] ifreq;
  delete[] Fouriercoeff;

  return S_x_new;
}

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
void func_fourierupper3d(int N, vector<double> &fc, vector<double> &f_up) {
  // Find f_up>abs(fc)
  double f_all[8];
  int kk[8];

  for (int kx = 0; kx < N; kx++) {
    int xr = kx + 1;
    if (kx == N - 1) xr = 0;

    for (int ky = 0; ky < N; ky++) {
      int yr = ky + 1;
      if (ky == N - 1) yr = 0;

      for (int kz = 0; kz < N; kz++) {
        int zr = kz + 1;
        if (kz == N - 1) zr = 0;

        kk[0] = kz + N * (ky + N * kx);
        kk[1] = zr + N * (ky + N * kx);
        kk[2] = kz + N * (yr + N * kx);
        kk[3] = zr + N * (yr + N * kx);
        kk[4] = kz + N * (ky + N * xr);
        kk[5] = zr + N * (ky + N * xr);
        kk[6] = kz + N * (yr + N * xr);
        kk[7] = zr + N * (yr + N * xr);

        double max_f_all = 0.;
        for (int k = 0; k < 8; k++) {
          f_all[k] = abs(fc[kk[k]]);
          max_f_all = max(max_f_all, f_all[k]);
        }

        f_up[kk[0]] = max_f_all;
      }
    }
  }
}

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
void sampleF(NeParticleGroup &S_x, double Neff_F_new, double Neff_F_old) {
  int Nf_old = S_x.size('f');
  // int Nf_new = myfloor((S_x -> size('p') + S_x -> size('n') )*resample_ratio
  // );

  // double Neff_F_new = Neff_F_old*Nf_old/Nf_new;

  int Nf_new = myfloor(Neff_F_old * Nf_old / Neff_F_new);

  // // cout << "Resample now " <<  Nf_new << ' ' << Nf_old << endl;

  if (Nf_new < Nf_old) {
    auto &Sfold = S_x.list('f');
    Particle1d3d *Sf = new Particle1d3d[Nf_new];

    const auto p = myrandperm(Nf_old, Nf_new);

    for (int kf = 0; kf < Nf_new; kf++) {
      Sf[kf].set_velocity(Sfold[p[kf] - 1].velocity());
      Sf[kf].set_position(Sfold[p[kf] - 1].position());
    }

    // enforce momentum conservation
    double uf_old[3] = {0., 0., 0.};
    double uf_new[3] = {0., 0., 0.};
    double vmod[3] = {0., 0., 0.};

    for (int kf = 0; kf < Nf_old; kf++) {
      auto &vf = Sfold[kf].velocity();
      for (int kv = 0; kv < 3; kv++) uf_old[kv] += vf[kv];
    }
    for (int kv = 0; kv < 3; kv++) uf_old[kv] *= Neff_F_old;

    for (int kf = 0; kf < Nf_new; kf++) {
      auto &vf = Sf[kf].velocity();
      for (int kv = 0; kv < 3; kv++) uf_new[kv] += vf[kv];
    }
    for (int kv = 0; kv < 3; kv++) uf_new[kv] *= Neff_F_new;

    for (int kv = 0; kv < 3; kv++)
      vmod[kv] = (uf_new[kv] - uf_old[kv]) / Neff_F_new / Nf_new;

    for (int kf = 0; kf < Nf_new; kf++) {
      auto vf = Sf[kf].velocity();
      for (int kv = 0; kv < 3; kv++) vf[kv] -= vmod[kv];
      Sf[kf].set_velocity(vf);
    }

    // enforce energy conservation
    // mu2f*T  + Nf*c^2 = Ep_old

    double c[3];
    for (int kv = 0; kv < 3; kv++) c[kv] = uf_old[kv] / Neff_F_old / Nf_old;

    double Told[3] = {0., 0., 0.};
    double Tnew[3] = {0., 0., 0.};
    double sigma[3] = {0., 0., 0.};  // sigma = sqrt(Told/Tnew)

    for (int kf = 0; kf < Nf_old; kf++) {
      auto &vf = Sfold[kf].velocity();
      for (int kv = 0; kv < 3; kv++)
        Told[kv] += (vf[kv] - c[kv]) * (vf[kv] - c[kv]);
    }
    for (int kv = 0; kv < 3; kv++) Told[kv] *= Neff_F_old;

    for (int kf = 0; kf < Nf_new; kf++) {
      auto &vf = Sf[kf].velocity();
      for (int kv = 0; kv < 3; kv++)
        Tnew[kv] += (vf[kv] - c[kv]) * (vf[kv] - c[kv]);
    }
    for (int kv = 0; kv < 3; kv++) Tnew[kv] *= Neff_F_new;

    for (int kv = 0; kv < 3; kv++) sigma[kv] = sqrt(Told[kv] / Tnew[kv]);

    for (int kf = 0; kf < Nf_new; kf++) {
      auto vf = Sf[kf].velocity();
      for (int kv = 0; kv < 3; kv++)
        vf[kv] = c[kv] + sigma[kv] * (vf[kv] - c[kv]);
      Sf[kf].set_velocity(vf);
    }

    // update F list
    S_x.clear('f');
    for (int kf = 0; kf < Nf_new; kf++) {
      S_x.push_back((Sf + kf), 'f');
    }

    // para.Neff_F = Neff_F_new;

    delete[] Sf;

  } else {
    cout << "CHECK F RESMAPLING!!!" << endl;
  }
}

int count_particle_number(const std::vector<NeParticleGroup> &S_x, int Nx,
                          char partype);

void sampleF_inhomo(std::vector<NeParticleGroup> &S_x, NumericGridClass &grid,
                    ParaClass &para) {
  int flag_resampled_tot = 0;
  for (int kx = 0; kx < grid.Nx; kx++) {
    flag_resampled_tot += S_x[kx].flag_resampled;
  }

  // // cout << " Resample F " << flag_resampled_tot << endl;

  if (flag_resampled_tot == grid.Nx) {
    int Np_tot = count_particle_number(S_x, grid.Nx, 'p');
    int Nn_tot = count_particle_number(S_x, grid.Nx, 'n');
    int Nf_tot = count_particle_number(S_x, grid.Nx, 'f');

    int Nf_tot_new = myfloor((Np_tot + Nn_tot) * para.resample_ratio);

    if (Nf_tot_new < Nf_tot) {
      double Neff_F_old = grid.Neff_F;
      double Neff_F_new = Neff_F_old * Nf_tot / Nf_tot_new;

      for (int kx = 0; kx < grid.Nx; kx++) {
        sampleF(S_x[kx], Neff_F_new, grid.Neff_F);
      }

      grid.Neff_F = Neff_F_new;
    }

    for (int kx = 0; kx < grid.Nx; kx++) {
      S_x[kx].reset_flag_resampled();
    }
  }
}

void addMaxwellian_terms(double rhoM, vector<double> uM, vector<double> TM,
                         double Neff, vector<double> &f, int Nfreq,
                         int augFactor, int orderx, int ordery, int orderz) {
  double Mcc_coe = rhoM / sqrt(8.0 * pi * pi * pi * TM[0] * TM[1] * TM[2]);

  vector<double> interp_xaug(Nfreq * augFactor);
  vector<double> exp_x(Nfreq * augFactor);
  vector<double> exp_y(Nfreq * augFactor);
  vector<double> exp_z(Nfreq * augFactor);
  for (int kx = 0; kx < Nfreq * augFactor; kx++) {
    double xk = kx * 2 * pi / Nfreq / augFactor;
    interp_xaug[kx] = xk;
    exp_x[kx] = exp(-(xk - uM[0]) * (xk - uM[0]) / 2 / TM[0]);
    exp_y[kx] = exp(-(xk - uM[1]) * (xk - uM[1]) / 2 / TM[1]);
    exp_z[kx] = exp(-(xk - uM[2]) * (xk - uM[2]) / 2 / TM[2]);
  }

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        int kk = kz + augFactor * Nfreq * (ky + augFactor * Nfreq * kx);

        double Mcc = Mcc_coe * exp_x[kx] * exp_y[ky] * exp_z[kz];
        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        if (orderx == 1) Mcc *= -(xc - uM[0]) / TM[0];
        if (orderx == 2)
          Mcc *= ((xc - uM[0]) * (xc - uM[0]) - TM[0]) / TM[0] / TM[0];
        if (ordery == 1) Mcc *= -(yc - uM[1]) / TM[1];
        if (ordery == 2)
          Mcc *= ((yc - uM[1]) * (yc - uM[1]) - TM[1]) / TM[1] / TM[1];
        if (orderz == 1) Mcc *= -(zc - uM[2]) / TM[2];
        if (orderz == 2)
          Mcc *= ((zc - uM[2]) * (zc - uM[2]) - TM[2]) / TM[2] / TM[2];

        f[kk] = Neff * f[kk] + Mcc;
      }
    }
  }
}

void addMaxwellian(double rhoM, vector<double> uM, vector<double> TM,
                   double Neff, vector<double> &f, vector<double> &fx,
                   vector<double> &fy, vector<double> &fz, vector<double> &fxx,
                   vector<double> &fyy, vector<double> &fzz,
                   vector<double> &fxy, vector<double> &fxz,
                   vector<double> &fyz, int Nfreq, int augFactor) {
  addMaxwellian_terms(rhoM, uM, TM, Neff, f, Nfreq, augFactor, 0, 0, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fx, Nfreq, augFactor, 1, 0, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fy, Nfreq, augFactor, 0, 1, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fz, Nfreq, augFactor, 0, 0, 1);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fxx, Nfreq, augFactor, 2, 0, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fyy, Nfreq, augFactor, 0, 2, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fzz, Nfreq, augFactor, 0, 0, 2);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fxy, Nfreq, augFactor, 1, 1, 0);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fxz, Nfreq, augFactor, 1, 0, 1);
  addMaxwellian_terms(rhoM, uM, TM, Neff, fyz, Nfreq, augFactor, 0, 1, 1);
}

NeParticleGroup resample_F_from_MPN(NeParticleGroup &S_x, int Nfreq,
                                    double Neff, double Neff_F,
                                    double dx_space) {
  NeParticleGroup S_x_new;
  /* Normalize particle velocity to [0 2*pi] */
  S_x.computexyzrange();

  NeParticleGroup S_x_renormalized;
  NeParticleGroup *ptr_S_x_renormalized = &S_x_renormalized;

  interp3d_renormalize(&S_x, ptr_S_x_renormalized);

  /* Prepare the grids in physical space and frequence space */
  // double dx = 2.0*pi/Nfreq;

  complex<double> *ifreq = new complex<double>[Nfreq];  // 1i *freq
  interp_ifreq(ifreq, Nfreq);
  vector<double> interp_x(Nfreq);
  for (int kx = 0; kx < Nfreq; kx++) interp_x[kx] = kx * 2 * pi / Nfreq;

  /* Compute the Fourier coefficient */
  complex<double> *Fouriercoeff = new complex<double>[Nfreq * Nfreq * Nfreq];
  vector<int> flag_Fouriercoeff(Nfreq * Nfreq * Nfreq);

  interp3d_fft_approx(ptr_S_x_renormalized, Fouriercoeff, Nfreq, Nfreq, Nfreq);
  filter_Fourier(
      Fouriercoeff, flag_Fouriercoeff,
      Nfreq * Nfreq * Nfreq);  // Apply the filter on Fourier coefficients

  // cout << " F coeff computed " << endl;

  /* Compute a coarse interpolation in physical space */
  vector<double> fcoarse(Nfreq * Nfreq * Nfreq);
  interp3d_fcoarse(Fouriercoeff, fcoarse, Nfreq, Nfreq, Nfreq);

  int augFactor = 2;
  int sizeF = augFactor * augFactor * augFactor * Nfreq * Nfreq * Nfreq;

  vector<double> f(sizeF);
  vector<double> fx(sizeF);
  vector<double> fy(sizeF);
  vector<double> fz(sizeF);
  vector<double> fxx(sizeF);
  vector<double> fyy(sizeF);
  vector<double> fzz(sizeF);
  vector<double> fxy(sizeF);
  vector<double> fxz(sizeF);
  vector<double> fyz(sizeF);

  interp3d_fxyz(Fouriercoeff, f, fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz,
                Nfreq, Nfreq, Nfreq, augFactor);

  vector<double> uM(3);
  vector<double> TM(3);
  double rhoM = S_x.rhoM * dx_space;
  uM[0] = ptr_S_x_renormalized->u1M;
  uM[1] = ptr_S_x_renormalized->u2M;
  uM[2] = ptr_S_x_renormalized->u3M;
  TM[0] = ptr_S_x_renormalized->T1M;
  TM[1] = ptr_S_x_renormalized->T2M;
  TM[2] = ptr_S_x_renormalized->T3M;

  addMaxwellian(rhoM, uM, TM, Neff, f, fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz,
                Nfreq, augFactor);

  // Add Maxwellian Here

  /* evaluate the upperbound of f */
  vector<double> f_up(sizeF);
  func_fourierupper3d(augFactor * Nfreq, f, f_up);

  /* refined x grid */
  double dxaug = 2.0 * pi / Nfreq / augFactor;
  vector<double> interp_xaug(Nfreq * augFactor);
  for (int kx = 0; kx < Nfreq * augFactor; kx++)
    interp_xaug[kx] = kx * 2 * pi / Nfreq / augFactor;

  /* create a NeParticleGroup to host the P and N particles in current cell */

  /* Start sampling */

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        int kk = kz + augFactor * Nfreq * (ky + augFactor * Nfreq * kx);

        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        double fcc = f_up[kk];

        double maxf = 1.5 * abs(fcc);
        int N_incell = myfloor(maxf * dxaug * dxaug * dxaug / Neff_F);

        int k_virtual = 0;
        NeParticleGroup S_x_incell;
        NeParticleGroup *ptr_S_x_incell = &S_x_incell;

        while (k_virtual < N_incell) {
          // create a particle in the cell
          // double Sf[3] = {xc+myrand()*dx, yc+myrand()*dx, zc+myrand()*dx};
          double deltax = myrand() * dxaug - 0.5 * dxaug;
          double deltay = myrand() * dxaug - 0.5 * dxaug;
          double deltaz = myrand() * dxaug - 0.5 * dxaug;
          double Sf[3] = {xc + deltax, yc + deltay, zc + deltaz};

          // compute f at this point
          double fval = interp3d_fvalue_approx(
              deltax, deltay, deltaz, f[kk], fx[kk], fy[kk], fz[kk], fxx[kk],
              fyy[kk], fzz[kk], fxy[kk], fxz[kk], fyz[kk]);

          // reset current cell if fval>maxf, otherwise continue sampling in
          // current cell
          resampleF_acceptsampled(Sf, ptr_S_x_incell, fval, maxf);

          // reset N_incell if maxf is changed
          N_incell = myfloor(maxf / (Neff_F / (dxaug * dxaug * dxaug)));
          k_virtual++;
        }

        mergeF_NeParticleGroup(S_x_new, S_x_incell);
      }
    }
  }

  // cout << "Resampled." << endl;

  // rescale to the original coordinates
  auto &Sp_sampled = S_x_new.list('f');
  double *xyz_minmax = S_x.xyz_minmax;
  interp3d_rescale(&Sp_sampled[0], S_x_new.size('f'), xyz_minmax);

  // cout << "Rescaled." << endl;

  // free memory
  delete[] ifreq;
  delete[] Fouriercoeff;
  return S_x_new;
}

void assign_positions(NeParticleGroup &S_new, double xmin, double xmax);
void merge_NeParticleGroup(NeParticleGroup &S_x,
                           const NeParticleGroup &S_x_new);
void mergeF_NeParticleGroup(NeParticleGroup &S_x,
                            const NeParticleGroup &S_x_new);
void save_particles(NeParticleGroup *S_x_before, NeParticleGroup *S_x_after);
// void particleresample_homo(NeParticleGroup * S_x, const ParaClass & para) {

// void particleresample_homo(NeParticleGroup * S_x, const ParaClass & para,
// MultlLevelGroup * MLsol) {
bool particleresample_homo(NeParticleGroup &S_x, const ParaClass &para) {
  // // cout << " resample 0" << endl;

  NUM_RESAMPLE++;

  int Np_old = S_x.size('p'), Nn_old = S_x.size('n');

  // int Nmax = 2*max(S_x -> size('p'), S_x -> size('n'));

  // cout << " resample 1" << endl;

  // resample particles
  auto S_x_new = samplefromfourier3d(S_x, para.Nfreq);
  // samplefromfourier3d(S_x, ptr_S_x_new, MLsol, para.Nlevel);

  int Np_new = S_x_new.size('p'), Nn_new = S_x_new.size('n');

  // cout << " Resample finished." << endl;
  // cout << "After resampling N = (" << ptr_S_x_new ->size('p') << ", " <<
  // ptr_S_x_new ->size('n') << ");" << endl;
  assign_positions(S_x_new, S_x.get_xmin(), S_x.get_xmax());

  S_x.flag_resampled = 1;

  // replace old particles by new sampled particles

  // save_particles(S_x, ptr_S_x_new);

  // Replace the original particles by new sampled particles
  if ((Np_new < Np_old) && (Nn_new < Nn_old)) {
    // cout << "Replace by new sampled particles" << endl;
    S_x.clear('p');
    S_x.clear('n');
    merge_NeParticleGroup(S_x, S_x_new);
    return true;
  } else {
    cout << "New sampled particles rejected." << endl;
    return false;
  }
}

void save_homo_dist(NeParticleGroup *S_x, const NumericGridClass &grid,
                    int flag_case);
void resampleF_inhomo(std::vector<NeParticleGroup> &S_x, double Neff_F_new,
                      NumericGridClass &grid, int Nfreq);

// void particleresample_inhomo(NeParticleGroup * S_x, NumericGridClass & grid,
// ParaClass & para, MultlLevelGroup * MLsol) {
void particleresample_inhomo(std::vector<NeParticleGroup> &S_x,
                             NumericGridClass &grid, ParaClass &para) {
  bool needGlobalResample = false;

  bool flag_resample_success = true;

  for (int kx = 0; kx < grid.Nx; kx++) {
    if ((S_x[kx].size('p') + S_x[kx].size('n')) >= S_x[kx].size('f'))
      needGlobalResample = true;
  }
  if (needGlobalResample) {
    double resample_spatial_ratio = 0;
    for (int kx = 0; kx < grid.Nx; kx++) {
      resample_spatial_ratio +=
          (S_x[kx].size('p') + S_x[kx].size('n')) / (S_x[kx].size('f'));
    }
    resample_spatial_ratio /= grid.Nx;

    resample_spatial_ratio = para.resample_spatial_ratio;

    // for (int kx = 0; kx < grid.Nx; kx ++) {
    int kx = 0;
    while ((flag_resample_success) && (kx < grid.Nx)) {
      if ((S_x[kx].size('p') + S_x[kx].size('n')) >=
          resample_spatial_ratio * S_x[kx].size('f')) {
        cout << "Particles resampling: ( " << S_x[kx].size('p') << ", "
             << S_x[kx].size('n') << ", " << S_x[kx].size('f') << ") " << endl;

        flag_resample_success = particleresample_homo(S_x[kx], para);

        cout << "After resampling: ( " << S_x[kx].size('p') << ", "
             << S_x[kx].size('n') << ", " << S_x[kx].size('f') << ") " << endl;
      }
      kx++;
    }
  }

  if (!flag_resample_success) {
    resampleF_inhomo(S_x, grid.Neff_F / 2, grid, para.Nfreq);

    int Nx = grid.Nx;
    vector<double> rho(Nx), rho_F(Nx);

    for (int kx = 0; kx < Nx; kx++) S_x[kx].computemoments();

    for (int kx = 0; kx < Nx; kx++)
      rho[kx] =
          S_x[kx].rhoM + (S_x[kx].m0P - S_x[kx].m0N) * grid.Neff / grid.dx;

    for (int kx = 0; kx < Nx; kx++)
      rho_F[kx] = S_x[kx].m0F * grid.Neff_F / grid.dx;

    save_macro<double>(rho, "rho_test");
    save_macro<double>(rho_F, "rhoF_test");
  }
}

void resampleF_homo(NeParticleGroup &S_x, double Neff_F_new, double Neff,
                    int Nfreq, double dx_space) {
  // resample particles
  auto S_x_new = resample_F_from_MPN(S_x, Nfreq, Neff, Neff_F_new, dx_space);

  assign_positions(S_x_new, S_x.get_xmin(), S_x.get_xmax());

  // replace old particles by new sampled particles

  S_x.clear('f');
  mergeF_NeParticleGroup(S_x, S_x_new);
}

void resampleF_inhomo(std::vector<NeParticleGroup> &S_x, double Neff_F_new,
                      NumericGridClass &grid, int Nfreq) {
  for (int kx = 0; kx < grid.Nx; kx++) {
    resampleF_homo(S_x[kx], Neff_F_new, grid.Neff, Nfreq, grid.dx);
  }

  grid.Neff_F = Neff_F_new;
  SYNC_TIME = 0;

  cout << "F particle resampled." << endl;
}

void resampleF_keeptotalmass(std::vector<NeParticleGroup> &S_x,
                             NumericGridClass &grid, int Nf_old) {
  int Nf_new = count_particle_number(S_x, grid.Nx, 'f');
  if (Nf_new > Nf_old) {
    double Neff_F_new = grid.Neff_F;
    double totalmass = 0;

    for (int kx = 0; kx < grid.Nx; kx++) {
      double mass = S_x[kx].rhoM * grid.dx;
      totalmass += mass;
      int Nk = (int)(mass / Neff_F_new);

      int Nk_remove = S_x[kx].size('f') - Nk;

      for (int kp = 0; kp < Nk_remove; kp++) {
        int k_remove = (int)(myrand() * S_x[kx].size('f'));
        S_x[kx].erase(k_remove, 'f');
      }
    }

    grid.Neff_F = totalmass / count_particle_number(S_x, grid.Nx, 'f');
  }
}

void sync_coarse(std::vector<NeParticleGroup> &S_x, NumericGridClass &grid,
                 ParaClass &para) {
  if (para.flag_collision == 1) {
    if (SYNC_TIME > para.sync_time_interval) {
      // cout << "Start resample F" << endl;

      // cout << "First resample P and N" << endl;
      SYNC_TIME = 0;
      particleresample_inhomo(S_x, grid, para);

      // cout << "P and N resampled" << endl;

      /*
      double totalmass = 0;
      for (int kx = 0; kx < grid.Nx; kx ++)  totalmass += (S_x+kx) -> rhoM;
      totalmass *= grid.dx;

      int Nd = 0;
      for (int kx = 0; kx < grid.Nx; kx ++)  Nd += ( (S_x+kx) -> size('p') +
      (S_x+kx) -> size('n') ); int Nf = (int)(1.2*Nd);

      double Neff_F_new = totalmass/Nf;
      if (Neff_F_new < grid.Neff_F) Neff_F_new = grid.Neff_F;
      */

      double Neff_F_new = 100;
      for (int kx = 0; kx < grid.Nx; kx++) {
        int N_one = (S_x[kx].size('p') + S_x[kx].size('n'));
        double Neff_F_one = (S_x[kx].rhoM) * grid.dx / N_one / 1.1;
        if (Neff_F_new > Neff_F_one) Neff_F_new = Neff_F_one;
      }

      if (Neff_F_new < grid.Neff_F) Neff_F_new = grid.Neff_F;

      // cout << "Next resample F" << endl;

      int Nf_old = count_particle_number(S_x, grid.Nx, 'f');

      resampleF_inhomo(S_x, Neff_F_new, grid, para.Nfreq);
      // cout << "F resampled" << endl;
      resampleF_keeptotalmass(S_x, grid, Nf_old);
    }
  }
}
}  // namespace coulomb