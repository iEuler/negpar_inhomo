#include "ParticleGroupOperations.h"

#include "ParticleGroup.h"
#include "RandomSampling.h"

namespace coulomb {

void ParticleGroupOperations::merge_signed(NeParticleGroup &S_x,
                                           const NeParticleGroup &S_x_new) {
  auto &Sp = S_x_new.list(ParticleKind::Positive);
  auto &Sn = S_x_new.list(ParticleKind::Negative);
  for (const auto &Sone : Sp)
    S_x.push_back(Sone, ParticleKind::Positive);
  for (const auto &Sone : Sn)
    S_x.push_back(Sone, ParticleKind::Negative);
}

void ParticleGroupOperations::merge_full(NeParticleGroup &S_x,
                                         const NeParticleGroup &S_x_new) {
  auto &Sf = S_x_new.list(ParticleKind::Full);
  for (const auto &Sone : Sf)
    S_x.push_back(Sone, ParticleKind::Full);
}

void ParticleGroupOperations::merge(NeParticleGroup &S_x,
                                    const NeParticleGroup &S_x_new,
                                    const std::vector<ParticleKind> &parTypes) {
  for (const auto kind : parTypes) {
    auto &Sp = S_x_new.list(kind);
    for (const auto &Sone : Sp)
      S_x.push_back(Sone, kind);
  }
}

void ParticleGroupOperations::assign_positions(NeParticleGroup &S_new,
                                               double xmin, double xmax,
                                               RandomContext &random) {
  double x1 = xmin, x2 = xmax;
  auto &Sp = S_new.list(ParticleKind::Positive);
  auto &Sn = S_new.list(ParticleKind::Negative);
  auto &Sf = S_new.list(ParticleKind::Full);
  for (int kp = 0; kp < S_new.size(ParticleKind::Positive); kp++)
    Sp[kp].set_position(RandomSampling::uniform(random) * (x2 - x1) + x1);
  for (int kp = 0; kp < S_new.size(ParticleKind::Negative); kp++)
    Sn[kp].set_position(RandomSampling::uniform(random) * (x2 - x1) + x1);
  for (int kp = 0; kp < S_new.size(ParticleKind::Full); kp++)
    Sf[kp].set_position(RandomSampling::uniform(random) * (x2 - x1) + x1);
}

} // namespace coulomb
