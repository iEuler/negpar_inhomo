#include "Classes.h"

namespace coulomb {
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
  auto Sp = S_x->list('p');
  auto Sn = S_x->list('n');

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

  auto Sp = S_x->list('p');
  auto Sn = S_x->list('n');

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
  for (int kp = 0; kp < Np; kp++) {
    auto &v0 = Sp[kp].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v0[k2] = (v0[k2] - xyz_minmax[2 * k2]) * 2.0 * pi / Lxyz[k2];
    }
    S_one.set_velocity(v0);
    S_x_new->push_back(S_one, 'p');
  }
  for (int kn = 0; kn < Nn; kn++) {
    auto &v0 = Sn[kn].velocity();
    for (int k2 = 0; k2 < 3; k2++) {
      v0[k2] = (v0[k2] - xyz_minmax[2 * k2]) * 2.0 * pi / Lxyz[k2];
    }
    S_one.set_velocity(v0);
    S_x_new->push_back(S_one, 'n');
  }
}

void interp3d_rescale(Particle1d3d *Sp, int Np, double *xyz_minmax) {
  // rescale to the original coordinates stored in xyz_minmax
  double Lxyz[3];
  for (int k2 = 0; k2 < 3; k2++) {
    Lxyz[k2] = xyz_minmax[2 * k2 + 1] - xyz_minmax[2 * k2];
  }

  for (int kp = 0; kp < Np; kp++) {
    auto &v0 = Sp[kp].velocity();
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

void interp_freq_aug(int *loc, int Nfreq, int augFactor) {
  for (int j = 0; j < Nfreq; j++) {
    int kfreq = freq_sequence(j, Nfreq);
    *(loc + j) = freq_sequence_reverse(kfreq, augFactor * Nfreq);
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

  auto Sp = S_x->list('p');
  auto Sn = S_x->list('n');

  double Neff = 1.0;

  // double Neff_temp = 1./Np;

  interp_ifreq(ifreq1, Nfreq1);
  interp_ifreq(ifreq2, Nfreq2);
  interp_ifreq(ifreq3, Nfreq3);

  double cubic_2pi = 8.0 * pi * pi * pi;
  // double dx = 2.0*pi/Nfreq;
  // double Neff = Neff_temp * Lxyz[0]*Lxyz[1]*Lxyz[2]/cubic_2pi; // need to
  // multiply para.Neff, here we take para.Neff = 1

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

void interp3d_fft_approx_terms(complex<double> *Fouriercoeff, double *f,
                               int Nfreq1, int Nfreq2, int Nfreq3,
                               int augFactor, int orderx, int ordery,
                               int orderz) {
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;

  double *freq1 = new double[Nfreq1];  // 1i *freq
  double *freq2 = new double[Nfreq2];  // 1i *freq
  double *freq3 = new double[Nfreq3];  // 1i *freq

  interp_freq(freq1, Nfreq1);
  interp_freq(freq2, Nfreq2);
  interp_freq(freq3, Nfreq3);

  int *loc1 = new int[Nfreq1];  // 1i *freq
  int *loc2 = new int[Nfreq2];  // 1i *freq
  int *loc3 = new int[Nfreq3];  // 1i *freq

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

  delete[] freq1;
  delete[] freq2;
  delete[] freq3;

  delete[] loc1;
  delete[] loc2;
  delete[] loc3;

  fftw_destroy_plan(plan3d_ft);
  fftw_free(fin);
  fftw_free(FSaug);
}

void interp3d_fft_approx(NeParticleGroup *S_x, complex<double> *Fouriercoeff,
                         int Nfreq1, int Nfreq2, int Nfreq3) {
  int augFactor = 4;

  int Np = S_x->size('p');
  int Nn = S_x->size('n');

  cout << "Number " << Np << " " << Nn << endl;
  auto Sp = S_x->list('p');
  auto Sn = S_x->list('n');

  double Neff = 1.0;

  // double Neff_temp = 1./Np;

  double cubic_2pi = 8.0 * pi * pi * pi;
  // double dx = 2.0*pi/Nfreq;
  // double Neff = Neff_temp * Lxyz[0]*Lxyz[1]*Lxyz[2]/cubic_2pi; // need to
  // multiply para.Neff, here we take para.Neff = 1

  double coeff_fft = 1. / cubic_2pi;
  double maxFS = 0.0;

  // for the (i,j,k)-th element of the array with size (Nx,Ny,Nz), use the
  // expression an_array[k + Nz * (j + Ny * i)].

  // create f, fx, fy, fz, fxx, fyy, fzz, fxy ...
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;
  double dx = 2.0 * pi / augFactor / Nfreq1;
  double dy = 2.0 * pi / augFactor / Nfreq2;
  double dz = 2.0 * pi / augFactor / Nfreq3;

  double *f = new double[sizeF];
  double *fx = new double[sizeF];
  double *fy = new double[sizeF];
  double *fz = new double[sizeF];
  double *fxx = new double[sizeF];
  double *fyy = new double[sizeF];
  double *fzz = new double[sizeF];
  double *fxy = new double[sizeF];
  double *fxz = new double[sizeF];
  double *fyz = new double[sizeF];

  for (int kk = 0; kk < sizeF; kk++) {
    *(f + kk) = 0.;
    *(fx + kk) = 0.;
    *(fy + kk) = 0.;
    *(fz + kk) = 0.;
    *(fxx + kk) = 0.;
    *(fyy + kk) = 0.;
    *(fzz + kk) = 0.;
    *(fxy + kk) = 0.;
    *(fxz + kk) = 0.;
    *(fyz + kk) = 0.;
  }

  for (int kp = 0; kp < Np; kp++) {
    double x0 = Sp[kp].velocity(0);
    double y0 = Sp[kp].velocity(1);
    double z0 = Sp[kp].velocity(2);
    int xloc = floor(x0 / dx);
    int yloc = floor(y0 / dy);
    int zloc = floor(z0 / dz);
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    int loc = zloc + augFactor * Nfreq3 * (yloc + augFactor * Nfreq2 * xloc);

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

  for (int kp = 0; kp < Nn; kp++) {
    double x0 = (Sn + kp)->velocity(0);
    double y0 = (Sn + kp)->velocity(1);
    double z0 = (Sn + kp)->velocity(2);
    int xloc = floor(x0 / dx);
    int yloc = floor(y0 / dy);
    int zloc = floor(z0 / dz);
    double xdelta = x0 - xloc * dx;
    double ydelta = y0 - yloc * dy;
    double zdelta = z0 - zloc * dz;

    int loc = zloc + augFactor * Nfreq3 * (yloc + augFactor * Nfreq2 * xloc);

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

  delete[] f;
  delete[] fx;
  delete[] fy;
  delete[] fz;
  delete[] fxx;
  delete[] fyy;
  delete[] fzz;
  delete[] fxy;
  delete[] fxz;
  delete[] fyz;
}

void filter_Fourier(complex<double> *Fouriercoeff, int *flag_Fouriercoeff,
                    int size_FC) {
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
                       complex<double> *ifreq3, int *flag_Fouriercoeff,
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

void interp3d_acceptsampled(double *Sf, Particle1d3d *Sp_incell,
                            Particle1d3d *Sn_incell, double fval, double &maxf,
                            double Neff_cell, int &Nskx, int &kskx,
                            int &kp_incell, int &kn_incell) {
  if (abs(fval) > maxf) {
    // keep sampled particles with rate maxf/maxf_new

    double keeprate = maxf / (1.5 * abs(fval));

    maxf = 1.5 * abs(fval);
    Nskx = myfloor(maxf / Neff_cell);
    // Nskx = myfloor(maxf*dx*dx*dx/Neff);
    int kkkk = 0;
    while (kkkk < kp_incell) {
      if (myrand() > keeprate) {
        auto &Sf_last = (Sp_incell + kp_incell - 1)->velocity();
        (Sp_incell + kkkk)->set_velocity(Sf_last);
        kp_incell--;
      } else {
        kkkk++;
      }
    }
    kkkk = 0;
    while (kkkk < kn_incell) {
      if (myrand() > keeprate) {
        auto &Sf_last = (Sn_incell + kn_incell - 1)->velocity();
        (Sn_incell + kkkk)->set_velocity(Sf_last);
        kn_incell--;
      } else {
        kkkk++;
      }
    }

  } else {
    // accept this particle with rate abs(fval/maxf)
    if (myrand() < (abs(fval / maxf))) {
      double sum_Sf_pi_sq = 0.;
      for (int kv = 0; kv < 3; kv++)
        sum_Sf_pi_sq += (Sf[kv] - pi) * (Sf[kv] - pi);
      if (sqrt(sum_Sf_pi_sq) < pi) {
        if (fval > 0) {
          (Sp_incell + kp_incell)->set_velocity(Sf);
          kp_incell++;
        } else {
          (Sn_incell + kn_incell)->set_velocity(Sf);
          kn_incell++;
        }
      }
    }
    kskx++;
  }
}

/******************************************************************/
/* ---------- Use Fourier transform for 3D interpolation -------- */
/******************************************************************/

/*
  Find the coarse approximation with the given Fourier coefficients
  Need to include 'fftw3.f'
*/
void interp3d_fcoarse(complex<double> *Fouriercoeff, double *fcoarse,
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

void interp3d_fxyz_terms(complex<double> *Fouriercoeff, double *f, int Nfreq1,
                         int Nfreq2, int Nfreq3, int augFactor, int orderx,
                         int ordery, int orderz) {
  int sizeF = augFactor * augFactor * augFactor * Nfreq1 * Nfreq2 * Nfreq3;
  complex<double> *FSaug = new complex<double>[sizeF];
  for (int k = 0; k < sizeF; k++) FSaug[k] = complex<double>(0., 0.);

  double *freq1 = new double[Nfreq1];  // 1i *freq
  double *freq2 = new double[Nfreq2];  // 1i *freq
  double *freq3 = new double[Nfreq3];  // 1i *freq

  interp_freq(freq1, Nfreq1);
  interp_freq(freq2, Nfreq2);
  interp_freq(freq3, Nfreq3);

  int *loc1 = new int[Nfreq1];  // 1i *freq
  int *loc2 = new int[Nfreq2];  // 1i *freq
  int *loc3 = new int[Nfreq3];  // 1i *freq

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

  delete[] freq1;
  delete[] freq2;
  delete[] freq3;

  delete[] FSaug;
  delete[] loc1;
  delete[] loc2;
  delete[] loc3;
}

void interp3d_fxyz(complex<double> *Fouriercoeff, double *f, double *fx,
                   double *fy, double *fz, double *fxx, double *fyy,
                   double *fzz, double *fxy, double *fxz, double *fyz,
                   int Nfreq1, int Nfreq2, int Nfreq3, int augFactor) {
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

void func_fourierupper3d(int N, double *fc, double *f_up);
void interp3d_fft_eachlevel(NeParticleGroup *S_x, MultlLevelGroup *MLsol,
                            int Nlevel);
void interp3d_fft_ml(complex<double> *Fouriercoeff, int *flag_Fouriercoeff,
                     MultlLevelGroup *MLsol, int Nlevel);
void save_complex(int Nx, complex<double> *FS, string filename);

/*
  Sampled particles are stored in particles_Sp_sampled and particles_Sn_sampled
  with size particles%Np_sampled and particles%Nn_sampled
*/

void samplefromfourier3d(NeParticleGroup *S_x, NeParticleGroup *S_x_new,
                         MultlLevelGroup *MLsol, int Nlevel) {
  // void samplefromfourier3d(NeParticleGroup * S_x, NeParticleGroup * S_x_new,
  // int Nfreq) {

  // some constants

  int Nfreq = 1 << Nlevel;
  double Neff = 1.0;

  // double Neff_temp = 1./( S_x-> size('p'));

  /* Normalize particle velocity to [0 2*pi] */
  S_x->computexyzrange();

  NeParticleGroup S_x_renormalized(S_x->size('p'), S_x->size('n'), 0);
  NeParticleGroup *ptr_S_x_new = &S_x_renormalized;

  interp3d_renormalize(S_x, ptr_S_x_new);

  /* Prepare the grids in physical space and frequence space */

  // double cubic_2pi = 8.0*pi*pi*pi;
  // double Lcubic = (double) (Nfreq*Nfreq*Nfreq);
  double dx = 2.0 * pi / Nfreq;
  // double Neff = Neff_temp*Lxyz[0]*Lxyz[1]*Lxyz[2]/cubic_2pi; // need to
  // multiply para.Neff, here we take para.Neff = 1

  complex<double> *ifreq = new complex<double>[Nfreq];  // 1i *freq
  interp_ifreq(ifreq, Nfreq);

  double *interp_x = new double[Nfreq];
  for (int kx = 0; kx < Nfreq; kx++) interp_x[kx] = kx * 2 * pi / Nfreq;

  /* Compute the Fourier coefficient */

  complex<double> *Fouriercoeff = new complex<double>[Nfreq * Nfreq * Nfreq];
  int *flag_Fouriercoeff = new int[Nfreq * Nfreq * Nfreq];

  int flag_multilevel = 0;
  if (flag_multilevel > 0) {
    /* multi level */
    interp3d_fft_eachlevel(ptr_S_x_new, MLsol, Nlevel);
    interp3d_fft_ml(Fouriercoeff, flag_Fouriercoeff, MLsol, Nlevel);
  } else {
    /* single level */
    // interp3d_fft(ptr_S_x_new, Fouriercoeff, Nfreq, Nfreq, Nfreq);
    interp3d_fft_approx(ptr_S_x_new, Fouriercoeff, Nfreq, Nfreq, Nfreq);
    filter_Fourier(
        Fouriercoeff, flag_Fouriercoeff,
        Nfreq * Nfreq * Nfreq);  // Apply the filter on Fourier coefficients
  }

  // cout << "Nfreq = " << Nfreq << endl;
  // save_complex(Nfreq*Nfreq*Nfreq, Fouriercoeff, "FS");

  // interp3d_fft_approx(ptr_S_x_new, Fouriercoeff, Nfreq, Nfreq, Nfreq);
  // save_complex(Nfreq*Nfreq*Nfreq, Fouriercoeff, "FS2");

  cout << " F coeff computed " << endl;

  /* Compute a coarse interpolation in physical space */
  // the value of f at grid points
  double *fcoarse = new double[Nfreq * Nfreq * Nfreq];
  interp3d_fcoarse(Fouriercoeff, fcoarse, Nfreq, Nfreq, Nfreq);
  // fx

  int augFactor = 2;
  int sizeF = augFactor * augFactor * augFactor * Nfreq * Nfreq * Nfreq;

  double *f = new double[sizeF];
  double *fx = new double[sizeF];
  double *fy = new double[sizeF];
  double *fz = new double[sizeF];
  double *fxx = new double[sizeF];
  double *fyy = new double[sizeF];
  double *fzz = new double[sizeF];
  double *fxy = new double[sizeF];
  double *fxz = new double[sizeF];
  double *fyz = new double[sizeF];

  interp3d_fxyz(Fouriercoeff, f, fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz,
                Nfreq, Nfreq, Nfreq, augFactor);

  /* extract P and N lists to hold the newly sampled particles */

  auto Sp_sampled = S_x_new->list('p');
  auto Sn_sampled = S_x_new->list('n');

  /* evaluate the upperbound of f */

  // double *f_up = new double [Nfreq*Nfreq*Nfreq];
  // func_fourierupper3d(Nfreq, fcoarse, f_up);
  double *f_up = new double[sizeF];
  func_fourierupper3d(augFactor * Nfreq, f, f_up);

  double dxaug = 2.0 * pi / Nfreq / augFactor;
  double *interp_xaug = new double[Nfreq * augFactor];
  for (int kx = 0; kx < Nfreq * augFactor; kx++)
    interp_xaug[kx] = kx * 2 * pi / Nfreq / augFactor;

  /* create two lists to host the P and N particles in current cell */

  double max_f_up = 0.;
  for (int kfc = 0; kfc < sizeF; kfc++) {
    max_f_up = max(max_f_up, abs(f_up[kfc]));
  }
  int N_incell = 3 * (ceil(max_f_up * dxaug * dxaug * dxaug / Neff));

  Particle1d3d *Sp_incell = new Particle1d3d[N_incell];
  Particle1d3d *Sn_incell = new Particle1d3d[N_incell];

  /* Start sampling */
  int kp = 0, kn = 0;

  for (int kx = 0; kx < augFactor * Nfreq; kx++) {
    for (int ky = 0; ky < augFactor * Nfreq; ky++) {
      for (int kz = 0; kz < augFactor * Nfreq; kz++) {
        int kk = kz + augFactor * Nfreq * (ky + augFactor * Nfreq * kx);

        double xc = interp_xaug[kx];
        double yc = interp_xaug[ky];
        double zc = interp_xaug[kz];

        double fcc = f_up[kk];

        // if (fcc<abs(fcoarse[kk]))  cout << "ERROR: small bound!" << endl;
        if (fcc < abs(f[kk])) cout << "ERROR: small bound!" << endl;

        double maxf = 1.5 * abs(fcc);
        int Nskx = myfloor(maxf * dxaug * dxaug * dxaug / Neff);

        bool flag_conti = true;
        if (Nskx == 0) flag_conti = false;

        // sum_Nk += Nskx;

        int kskx = 0;
        int kp_incell = 0;
        int kn_incell = 0;

        while (flag_conti) {
          // create a particle in the cell
          // double Sf[3] = {xc+myrand()*dx, yc+myrand()*dx, zc+myrand()*dx};
          double deltax = myrand() * dxaug - 0.5 * dxaug;
          double deltay = myrand() * dxaug - 0.5 * dxaug;
          double deltaz = myrand() * dxaug - 0.5 * dxaug;
          double Sf[3] = {xc + deltax, yc + deltay, zc + deltaz};

          // compute f at this point
          // double fval = interp3d_fvalue(Sf, Fouriercoeff, ifreq, ifreq,
          // ifreq, flag_Fouriercoeff, Nfreq, Nfreq, Nfreq);

          double fval = interp3d_fvalue_approx(
              deltax, deltay, deltaz, f[kk], fx[kk], fy[kk], fz[kk], fxx[kk],
              fyy[kk], fzz[kk], fxy[kk], fxz[kk], fyz[kk]);

          // reset current cell if fval>maxf, otherwise continue sampling in
          // current cell

          interp3d_acceptsampled(Sf, Sp_incell, Sn_incell, fval, maxf,
                                 Neff / (dxaug * dxaug * dxaug), Nskx, kskx,
                                 kp_incell, kn_incell);

          if (kskx >= Nskx) {
            flag_conti = false;
            for (int kp_c = 0; kp_c < kp_incell; kp_c++) {
              Sp_sampled[kp + kp_c].set_velocity(
                  (Sp_incell + kp_c)->velocity());
            }
            for (int kn_c = 0; kn_c < kn_incell; kn_c++) {
              Sn_sampled[kn + kn_c].set_velocity(
                  (Sn_incell + kn_c)->velocity());
            }
            kp += kp_incell;
            kn += kn_incell;
          }
        }
      }
    }
  }

  cout << "Resampled." << endl;
  S_x_new->update_NpNn(kp, kn);

  // rescale to the original coordinates
  double *xyz_minmax = S_x->xyz_minmax;
  interp3d_rescale(Sp_sampled, kp, xyz_minmax);
  interp3d_rescale(Sn_sampled, kn, xyz_minmax);

  cout << "Rescaled." << endl;

  // free memory
  delete[] f_up;
  delete[] interp_x;
  delete[] interp_xaug;
  delete[] ifreq;
  delete[] fcoarse;
  delete[] Sp_incell;
  delete[] Sn_incell;
  delete[] Fouriercoeff;
  delete[] flag_Fouriercoeff;

  delete[] f;
  delete[] fx;
  delete[] fy;
  delete[] fz;
  delete[] fxx;
  delete[] fyy;
  delete[] fzz;
  delete[] fxy;
  delete[] fxz;
  delete[] fyz;
}

/******************************************************************/
/* ------ Find an upper bound the for interpolated function ----- */
/******************************************************************/
void func_fourierupper3d(int N, double *fc, double *f_up) {
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
void sampleF(NeParticleGroup *S_x, double Neff_F_new, const ParaClass &para) {
  int Nf_old = S_x->size('f');
  // int Nf_new = myfloor((S_x -> size('p') + S_x -> size('n') )*resample_ratio
  // );

  double Neff_F_old = para.Neff_F;
  // double Neff_F_new = Neff_F_old*Nf_old/Nf_new;

  int Nf_new = myfloor(Neff_F_old * Nf_old / Neff_F_new);

  // // cout << "Resample now " <<  Nf_new << ' ' << Nf_old << endl;

  if (Nf_new < Nf_old) {
    auto Sfold = S_x->list('f');
    Particle1d3d *Sf = new Particle1d3d[Nf_new];

    vector<int> p(Nf_new);
    myrandperm(Nf_old, Nf_new, p);

    for (int kf = 0; kf < Nf_new; kf++) {
      (Sf + kf)->set_velocity(Sfold[p[kf] - 1].velocity());
      (Sf + kf)->set_position((Sfold[p[kf] - 1].position());
    }

    // enforce momentum conservation
    double uf_old[3] = {0., 0., 0.};
    double uf_new[3] = {0., 0., 0.};
    double vmod[3] = {0., 0., 0.};

    for (int kf = 0; kf < Nf_old; kf++) {
      auto &vf = (Sfold + kf)->velocity();
      for (int kv = 0; kv < 3; kv++) uf_old[kv] += vf[kv];
    }
    for (int kv = 0; kv < 3; kv++) uf_old[kv] *= Neff_F_old;

    for (int kf = 0; kf < Nf_new; kf++) {
      auto &vf = (Sf + kf)->velocity();
      for (int kv = 0; kv < 3; kv++) uf_new[kv] += vf[kv];
    }
    for (int kv = 0; kv < 3; kv++) uf_new[kv] *= Neff_F_new;

    for (int kv = 0; kv < 3; kv++)
      vmod[kv] = (uf_new[kv] - uf_old[kv]) / Neff_F_new / Nf_new;

    for (int kf = 0; kf < Nf_new; kf++) {
      auto &vf = (Sf + kf)->velocity();
      for (int kv = 0; kv < 3; kv++) vf[kv] -= vmod[kv];
      (Sf + kf)->set_velocity(vf);
    }

    // enforce energy conservation
    // mu2f*T  + Nf*c^2 = Ep_old

    double c[3];
    for (int kv = 0; kv < 3; kv++) c[kv] = uf_old[kv] / Neff_F_old / Nf_old;

    double Told[3] = {0., 0., 0.};
    double Tnew[3] = {0., 0., 0.};
    double sigma[3] = {0., 0., 0.};  // sigma = sqrt(Told/Tnew)

    for (int kf = 0; kf < Nf_old; kf++) {
      auto &vf = (Sfold + kf)->velocity();
      for (int kv = 0; kv < 3; kv++)
        Told[kv] += (vf[kv] - c[kv]) * (vf[kv] - c[kv]);
    }
    for (int kv = 0; kv < 3; kv++) Told[kv] *= Neff_F_old;

    for (int kf = 0; kf < Nf_new; kf++) {
      auto &vf = (Sf + kf)->velocity();
      for (int kv = 0; kv < 3; kv++)
        Tnew[kv] += (vf[kv] - c[kv]) * (vf[kv] - c[kv]);
    }
    for (int kv = 0; kv < 3; kv++) Tnew[kv] *= Neff_F_new;

    for (int kv = 0; kv < 3; kv++) sigma[kv] = sqrt(Told[kv] / Tnew[kv]);

    for (int kf = 0; kf < Nf_new; kf++) {
      auto &vf = (Sf + kf)->velocity();
      for (int kv = 0; kv < 3; kv++)
        vf[kv] = c[kv] + sigma[kv] * (vf[kv] - c[kv]);
      (Sf + kf)->set_velocity(vf);
    }

    // update F list
    S_x->update_Nf(0);
    for (int kf = 0; kf < Nf_new; kf++) {
      S_x->push_back((Sf + kf), 'f');
    }

    // para.Neff_F = Neff_F_new;

    delete[] Sf;

  } else {
    cout << "CHECK F RESMAPLING!!!" << endl;
  }
}

int count_particle_number(const std::vector<NeParticleGroup> &S_x, int Nx,
                          char partype);

void sampleF_inhomo(NeParticleGroup *S_x, NumericGridClass &grid,
                    ParaClass &para) {
  int flag_resampled_tot = 0;
  for (int kx = 0; kx < grid.Nx; kx++) {
    flag_resampled_tot += (S_x + kx)->flag_resampled;
  }

  // // cout << " Resample F " << flag_resampled_tot << endl;

  if (flag_resampled_tot == grid.Nx) {
    int Np_tot = count_particle_number(S_x, grid.Nx, 'p');
    int Nn_tot = count_particle_number(S_x, grid.Nx, 'n');
    int Nf_tot = count_particle_number(S_x, grid.Nx, 'f');

    int Nf_tot_new = myfloor((Np_tot + Nn_tot) * para.resample_ratio);

    if (Nf_tot_new < Nf_tot) {
      double Neff_F_old = para.Neff_F;
      double Neff_F_new = Neff_F_old * Nf_tot / Nf_tot_new;

      for (int kx = 0; kx < grid.Nx; kx++) {
        sampleF(S_x + kx, Neff_F_new, para);
      }

      grid.Neff_F = Neff_F_new;
      para.Neff_F = Neff_F_new;
    }

    for (int kx = 0; kx < grid.Nx; kx++) {
      (S_x + kx)->reset_flag_resampled();
    }
  }
}

void assign_positions(NeParticleGroup *S_new, double xmin, double xmax);
void merge_NeParticleGroup(NeParticleGroup *S_x, NeParticleGroup *S_x_new);

// void particleresample_homo(NeParticleGroup * S_x, const ParaClass & para) {

void particleresample_homo(NeParticleGroup *S_x, const ParaClass &para,
                           MultlLevelGroup *MLsol) {
  // // cout << " resample 0" << endl;

  int Nmax = max(S_x->size('p'), S_x->size('n'));
  NeParticleGroup S_x_new(Nmax, Nmax, 0);
  NeParticleGroup *ptr_S_x_new = &S_x_new;

  // cout << " resample 1" << endl;

  // resample particles
  // samplefromfourier3d(S_x, ptr_S_x_new, para.Nfreq);
  samplefromfourier3d(S_x, ptr_S_x_new, MLsol, para.Nlevel);

  // cout << " Resample finished." << endl;

  assign_positions(ptr_S_x_new, S_x->get_xmin(), S_x->get_xmax());

  S_x->flag_resampled = 1;

  S_x->update_NpNn(0, 0);

  merge_NeParticleGroup(S_x, ptr_S_x_new);

  // cout << " Merge finished." << endl;
}

void save_homo_dist(NeParticleGroup *S_x, const NumericGridClass &grid,
                    int flag_case);

void particleresample_inhomo(NeParticleGroup *S_x, NumericGridClass &grid,
                             ParaClass &para, MultlLevelGroup *MLsol) {
  for (int kx = 0; kx < grid.Nx; kx++) {
    if (((S_x + kx)->size('p') + (S_x + kx)->size('n')) >=
        (S_x + kx)->size('f')) {
      cout << "Particles resampling: ( " << (S_x + kx)->size('p') << ", "
           << (S_x + kx)->size('n') << ", " << (S_x + kx)->size('f') << ") "
           << endl;

      // FLAG_FILENAME_WITH_NUM = false;
      save_homo_dist(S_x, grid, 1);
      particleresample_homo(S_x + kx, para, MLsol);
      save_homo_dist(S_x, grid, 2);

      cout << "After resampling: ( " << (S_x + kx)->size('p') << ", "
           << (S_x + kx)->size('n') << ", " << (S_x + kx)->size('f') << ") "
           << endl;

      // FLAG_FILENAME_WITH_NUM = true;
    }
  }
  /*
  int flag_resampled_tot = 0;
  for (int kx = 0; kx < grid.Nx; kx ++) {
    flag_resampled_tot += (S_x+kx) -> flag_resampled;
  }
  // cout << "Resampled " << flag_resampled_tot << endl;

  sampleF_inhomo(S_x, grid, para);
  */
}

void particleresample_inhomo_nocoll(NeParticleGroup *S_x,
                                    NumericGridClass &grid, ParaClass &para,
                                    MultlLevelGroup *MLsol) {
  if ((count_particle_number(S_x, grid.Nx, 'p') +
       count_particle_number(S_x, grid.Nx, 'n')) >= (10000 * grid.Nx)) {
    for (int kx = 0; kx < grid.Nx; kx++) {
      cout << "Particles resampling: ( " << (S_x + kx)->size('p') << ", "
           << (S_x + kx)->size('n') << ", " << (S_x + kx)->size('f') << ") "
           << endl;

      particleresample_homo(S_x + kx, para, MLsol);

      cout << "After resampling: ( " << (S_x + kx)->size('p') << ", "
           << (S_x + kx)->size('n') << ", " << (S_x + kx)->size('f') << ") "
           << endl;
    }
  }
}

/* Compute the Fourier coefficients on each level*/
void interp3d_fft_eachlevel(NeParticleGroup *S_x, MultlLevelGroup *MLsol,
                            int Nlevel) {
  int Num_grids = (3 * Nlevel * Nlevel - 3 * Nlevel + 2) / 2;
  for (int k_grids = 0; k_grids < Num_grids; k_grids++) {
    int Nx = (MLsol + k_grids)->Nx;
    int Ny = (MLsol + k_grids)->Ny;
    int Nz = (MLsol + k_grids)->Nz;
    complex<double> *Fouriercoeff = (MLsol + k_grids)->complexsol;
    interp3d_fft(S_x, Fouriercoeff, Nx, Ny, Nz);
  }
}

void interp3d_fft_ml(complex<double> *Fouriercoeff, int *flag_Fouriercoeff,
                     MultlLevelGroup *MLsol, int Nlevel) {
  int Nfine = 1 << Nlevel;
  for (int kk = 0; kk < Nfine * Nfine * Nfine; kk++)
    Fouriercoeff[kk] = complex<double>(0., 0.);
  for (int kk = 0; kk < Nfine * Nfine * Nfine; kk++) flag_Fouriercoeff[kk] = 0;

  int Num_grids = (3 * Nlevel * Nlevel - 3 * Nlevel + 2) / 2;
  for (int k_grids = 0; k_grids < Num_grids; k_grids++) {
    int Nx = (MLsol + k_grids)->Nx;
    int Ny = (MLsol + k_grids)->Ny;
    int Nz = (MLsol + k_grids)->Nz;

    int *position = (MLsol + k_grids)->position;
    double coe = (MLsol + k_grids)->coe_level;
    complex<double> *Fouriercoeff_kgrids = (MLsol + k_grids)->complexsol;

    for (int k1 = 0; k1 < Nx; k1++) {
      for (int k2 = 0; k2 < Ny; k2++) {
        for (int k3 = 0; k3 < Nz; k3++) {
          int kcoarse = k3 + Nz * (k2 + Ny * k1);
          int kfine = position[kcoarse];

          Fouriercoeff[kfine] += coe * Fouriercoeff_kgrids[kcoarse];
          flag_Fouriercoeff[kfine] += 1;
        }
      }
    }
  }
}

/* Compute the fcoarse on each level*/
/*
void interp3d_fcoarse_eachlevel(MultlLevelGroup * MLsol, int Nlevel) {
  int Num_grids = (3*Nlevel*Nlevel - 3*Nlevel + 2)/2;
  for (int k_grids = 0; k_grids < Num_grids; k_grids ++) {
    int Nx = (MLsol + k_grids) -> Nx;
    int Ny = (MLsol + k_grids) -> Ny;
    int Nz = (MLsol + k_grids) -> Nz;
    complex<double> *Fouriercoeff = (MLsol + k_grids) -> complexsol;
    double *fcoarse = (MLsol + k_grids) -> sol;

    interp3d_fcoarse(Fouriercoeff, fcoarse, Nx, Ny, Nz);
  }

}

*/
}  // namespace coulomb