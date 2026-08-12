#include "ParticleGroup.h"

#include <algorithm>
#include <stdexcept>

#include "RandomSampling.h"

namespace coulomb {

void ParticleGroup::set_xrange(double x1, double x2) {
  xmin = x1;
  xmax = x2;
}
void ParticleGroup::push_back(const Particle1d3d &particle) {
  vS.push_back(particle);
}

Particle1d3d &ParticleGroup::list(int index) {
  if (index < 0 || static_cast<std::size_t>(index) >= vS.size())
    throw std::out_of_range("particle index is outside the group");
  return vS[index];
}

const Particle1d3d &ParticleGroup::list(int index) const {
  if (index < 0 || static_cast<std::size_t>(index) >= vS.size())
    throw std::out_of_range("particle index is outside the group");
  return vS[index];
}

void ParticleGroup::erase(int index) {
  if (index < 0 || static_cast<std::size_t>(index) >= vS.size())
    throw std::out_of_range("particle index is outside the group");
  vS[index] = vS.back();
  vS.pop_back();
}

void ParticleGroup::computemoments() {
  moments.m0 = static_cast<double>(vS.size());
  moments.m11 = 0.0;
  moments.m12 = 0.0;
  moments.m13 = 0.0;
  moments.m21 = 0.0;
  moments.m22 = 0.0;
  moments.m23 = 0.0;
  for (const auto &particle : vS) {
    const double v1 = particle.velocity(0);
    const double v2 = particle.velocity(1);
    const double v3 = particle.velocity(2);
    moments.m11 += v1;
    moments.m12 += v2;
    moments.m13 += v3;
    moments.m21 += v1 * v1;
    moments.m22 += v2 * v2;
    moments.m23 += v3 * v3;
  }
  moments.m2 = moments.m21 + moments.m22 + moments.m23;
}

void NeParticleGroup::set_xrange(double x1, double x2) {
  xmin = x1;
  xmax = x2;
}

int NeParticleGroup::size(ParticleKind kind) const {
  return static_cast<int>(list(kind).size());
}

void NeParticleGroup::push_back(const Particle1d3d &particle,
                                ParticleKind kind) {
  list(kind).push_back(particle);
}

void NeParticleGroup::erase(int index, ParticleKind kind) {
  auto &particles = list(kind);
  if (index < 0 || static_cast<std::size_t>(index) >= particles.size())
    throw std::out_of_range("particle index is outside the selected group");
  particles[index] = particles.back();
  particles.pop_back();
}

void NeParticleGroup::clear(ParticleKind kind) { list(kind).clear(); }

void NeParticleGroup::computemoments() {
  const auto compute = [](const std::vector<Particle1d3d> &particles) {
    Moments result;
    result.m0 = static_cast<double>(particles.size());
    for (const auto &particle : particles) {
      const double v1 = particle.velocity(0);
      const double v2 = particle.velocity(1);
      const double v3 = particle.velocity(2);
      const double squared_speed = v1 * v1 + v2 * v2 + v3 * v3;
      result.m11 += v1;
      result.m12 += v2;
      result.m13 += v3;
      result.m21 += v1 * v1;
      result.m22 += v2 * v2;
      result.m23 += v3 * v3;
      result.m31 += v1 * squared_speed;
      result.m32 += v2 * squared_speed;
      result.m33 += v3 * squared_speed;
    }
    result.m2 = result.m21 + result.m22 + result.m23;
    return result;
  };
  positive_moments = compute(vSp);
  negative_moments = compute(vSn);
  full_moments = compute(vSf);
}

void NeParticleGroup::copymoments() {
  previous_positive_moments = positive_moments;
  previous_negative_moments = negative_moments;
  rho_o = rho;
  u1_o = u1;
  Tprt_o = Tprt;
}

void NeParticleGroup::setPositionRangeAndRandomizeValues(
    double x1, double x2, RandomContext &random) {
  xmin = x1;
  xmax = x2;
  for (auto kind :
       {ParticleKind::Positive, ParticleKind::Negative, ParticleKind::Full}) {
    for (auto &particle : list(kind))
      particle.set_position(RandomSampling(random).uniform() * (xmax - xmin) +
                            xmin);
  }
}

void NeParticleGroup::set_xyzrange() {
  set_xyzrange({ParticleKind::Positive, ParticleKind::Negative});
}

void NeParticleGroup::set_xyzrange(std::initializer_list<ParticleKind> kinds) {
  std::fill(xyz_minmax.begin(), xyz_minmax.end(), 0.0);
  for (const auto kind : kinds) {
    for (const auto &particle : list(kind)) {
      const auto &velocity = particle.velocity();
      for (int component = 0; component < 3; ++component) {
        xyz_minmax[2 * component] =
            std::min(xyz_minmax[2 * component], velocity[component]);
        xyz_minmax[2 * component + 1] =
            std::max(xyz_minmax[2 * component + 1], velocity[component]);
      }
    }
  }
  for (int component = 0; component < 3; ++component) {
    xyz_minmax[2 * component] -= 1e-6;
    xyz_minmax[2 * component + 1] += 1e-6;
  }
}

std::vector<Particle1d3d> &NeParticleGroup::list(ParticleKind kind) {
  switch (kind) {
  case ParticleKind::Positive:
    return vSp;
  case ParticleKind::Negative:
    return vSn;
  case ParticleKind::Full:
    return vSf;
  }
  throw std::invalid_argument("unknown particle kind");
}

const std::vector<Particle1d3d> &
NeParticleGroup::list(ParticleKind kind) const {
  switch (kind) {
  case ParticleKind::Positive:
    return vSp;
  case ParticleKind::Negative:
    return vSn;
  case ParticleKind::Full:
    return vSf;
  }
  throw std::invalid_argument("unknown particle kind");
}

Particle1d3d &NeParticleGroup::list(int index, ParticleKind kind) {
  auto &particles = list(kind);
  if (index < 0 || static_cast<std::size_t>(index) >= particles.size())
    throw std::out_of_range("particle index is outside the selected group");
  return particles[index];
}

const Particle1d3d &NeParticleGroup::list(int index, ParticleKind kind) const {
  const auto &particles = list(kind);
  if (index < 0 || static_cast<std::size_t>(index) >= particles.size())
    throw std::out_of_range("particle index is outside the selected group");
  return particles[index];
}

} // namespace coulomb
