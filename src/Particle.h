#pragma once

#include <array>
#include <stdexcept>
#include <vector>

namespace coulomb {

enum class ParticleKind { Positive, Negative, Full };

ParticleKind particle_kind_from_code(char code);
char particle_kind_code(ParticleKind kind);

class Particle1d3d {
 public:
  Particle1d3d() : flag_moved(false) {}
  Particle1d3d(const std::vector<double>& velocity)
      : x(0.0), v(velocity), flag_moved(false) { validate_velocity(v); }
  Particle1d3d(double position, const std::vector<double>& velocity)
      : x(position), v(velocity), flag_moved(false) { validate_velocity(v); }
  Particle1d3d(double position, const std::vector<double>& velocity, bool moved)
      : x(position), v(velocity), flag_moved(moved) { validate_velocity(v); }

  void set_position(double position) { x = position; }
  void set_velocity(const std::vector<double>& velocity) {
    validate_velocity(velocity);
    v = velocity;
  }
  void set_velocity(const std::array<double, 3>& velocity) {
    v.assign(velocity.begin(), velocity.end());
  }
  double position() const { return x; }
  const std::vector<double>& velocity() const { return v; }
  double velocity(int component) const {
    if (component < 0 || component >= 3)
      throw std::out_of_range("velocity component is outside [0, 3)");
    return v[component];
  }

  bool flag_moved;

 private:
  static void validate_velocity(const std::vector<double>& velocity) {
    if (velocity.size() != 3)
      throw std::invalid_argument("particle velocity must have 3 components");
  }

  double x{0.0};
  std::vector<double> v{0.0, 0.0, 0.0};
};

}  // namespace coulomb
