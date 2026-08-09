#include "Particle.h"

namespace coulomb {

char particle_kind_code(ParticleKind kind) {
  switch (kind) {
    case ParticleKind::Positive: return 'p';
    case ParticleKind::Negative: return 'n';
    case ParticleKind::Full: return 'f';
  }
  throw std::invalid_argument("unknown particle kind");
}

ParticleKind particle_kind_from_code(char code) {
  switch (code) {
    case 'p': return ParticleKind::Positive;
    case 'n': return ParticleKind::Negative;
    case 'f': return ParticleKind::Full;
    default:
      throw std::invalid_argument("particle kind must be 'p', 'n', or 'f'");
  }
}

}  // namespace coulomb
