#include "ElectricField.h"

#include <complex>

#include "Constants.h"
#include "FFT.h"
#include "Grid.h"
#include "ParticleGroup.h"

namespace coulomb {

std::vector<double>
ElectricFieldSolver::solve_poisson(const std::vector<double> &rho,
                                   double lambda) {
  const int grid_size = grid_.Nx;
  const double domain_size = grid_.xmax - grid_.xmin;
  std::vector<double> electric_field(grid_size);
  FFT1D fft_calculator(grid_size);
  const auto rho_fft = fft_calculator.fft(rho);

  std::vector<double> wave_numbers(grid_size);
  for (int index = 0; index < grid_size / 2 + 1; ++index)
    wave_numbers[index] = index * 2.0 * pi / domain_size;
  for (int index = grid_size / 2 + 1; index < grid_size; ++index)
    wave_numbers[index] = (index - grid_size) * 2.0 * pi / domain_size;

  std::vector<std::complex<double>> electric_fft(grid_size);
  electric_fft[0] = {0.0, 0.0};
  for (int index = 1; index < grid_size; ++index) {
    electric_fft[index] = {
        rho_fft[index].imag() / wave_numbers[index] / grid_size,
        -rho_fft[index].real() / wave_numbers[index] / grid_size};
  }

  const auto electric = fft_calculator.ifft(electric_fft);
  for (int index = 0; index < grid_size; ++index)
    electric_field[index] = lambda * electric[index];
  return electric_field;
}

void ElectricFieldSolver::update(std::vector<ParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].computemoments();

  std::vector<double> rho(grid.Nx);
  for (int cell = 0; cell < grid.Nx; ++cell)
    rho[cell] = groups[cell].moments.m0 * grid.Neff / grid.dx;

  const auto electric_field = solve_poisson(rho, 10.0);
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].elecfield = electric_field[cell];
}

void ElectricFieldSolver::update(std::vector<NeParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].computemoments();

  std::vector<double> rho(grid.Nx);
  for (int cell = 0; cell < grid.Nx; ++cell)
    rho[cell] = groups[cell].rhoM + (groups[cell].positive_moments.m0 -
                                     groups[cell].negative_moments.m0) *
                                        grid.Neff / grid.dx;

  auto electric_field = solve_poisson(rho, grid.lambda_Poisson);
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].elecfield = electric_field[cell];

  for (int cell = 0; cell < grid.Nx; ++cell)
    rho[cell] = groups[cell].full_moments.m0 * grid.Neff_F / grid.dx;
  electric_field = solve_poisson(rho, grid.lambda_Poisson);
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].elecfield_F = electric_field[cell];
}

void ElectricFieldSolver::update_pic(std::vector<NeParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].computemoments();

  std::vector<double> rho(grid.Nx);
  for (int cell = 0; cell < grid.Nx; ++cell)
    rho[cell] = groups[cell].full_moments.m0 * grid.Neff_F / grid.dx;
  const auto electric_field = solve_poisson(rho, grid.lambda_Poisson);
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].elecfield_F = electric_field[cell];
}

void ElectricFieldSolver::update_from_coarse(
    std::vector<NeParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].computemoments();

  std::vector<double> rho(grid.Nx);
  for (int cell = 0; cell < grid.Nx; ++cell)
    rho[cell] = groups[cell].full_moments.m0 * grid.Neff_F / grid.dx;
  const auto electric_field = solve_poisson(rho, grid.lambda_Poisson);
  for (int cell = 0; cell < grid.Nx; ++cell)
    groups[cell].elecfield = electric_field[cell];
}

void ElectricFieldSolver::clear(std::vector<NeParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell) {
    groups[cell].computemoments();
    groups[cell].elecfield = 0.0;
    groups[cell].elecfield_F = 0.0;
  }
}

void ElectricFieldSolver::update_from_density(
    std::vector<NeParticleGroup> &groups) {
  const auto &grid = grid_;
  for (int cell = 0; cell < grid.Nx; ++cell) {
    groups[cell].computemoments();
    groups[cell].elecfield = groups[cell].rho;
    groups[cell].elecfield_F = groups[cell].rhoF;
  }
}

} // namespace coulomb
