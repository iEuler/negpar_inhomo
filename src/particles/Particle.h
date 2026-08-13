#pragma once

#include <array>
#include <stdexcept>
#include <vector>

namespace coulomb {

enum class ParticleKind { Positive, Negative, Full };

class ParticleKindCodec {
  public:
	ParticleKind decode(char code) const;
	char encode(ParticleKind kind) const;
};

class Particle1D3D {
  public:
	Particle1D3D() : flagMoved(false) {}
	Particle1D3D(const std::vector<double>& velocity)
		: x(0.0), v(velocity), flagMoved(false) {
		validateVelocity(v);
	}
	Particle1D3D(double position, const std::vector<double>& velocity)
		: x(position), v(velocity), flagMoved(false) {
		validateVelocity(v);
	}
	Particle1D3D(double position, const std::vector<double>& velocity,
				 bool moved)
		: x(position), v(velocity), flagMoved(moved) {
		validateVelocity(v);
	}

	void setPosition(double position) { x = position; }
	void setVelocity(const std::vector<double>& velocity) {
		validateVelocity(velocity);
		v = velocity;
	}
	void setVelocity(const std::array<double, 3>& velocity) {
		v.assign(velocity.begin(), velocity.end());
	}
	double position() const { return x; }
	const std::vector<double>& velocity() const { return v; }
	double velocity(int component) const {
		if (component < 0 || component >= 3)
			throw std::out_of_range("velocity component is outside [0, 3)");
		return v[component];
	}

	bool flagMoved;

  private:
	static void validateVelocity(const std::vector<double>& velocity) {
		if (velocity.size() != 3)
			throw std::invalid_argument(
				"particle velocity must have 3 components");
	}

	double x{0.0};
	std::vector<double> v{0.0, 0.0, 0.0};
};

} // namespace coulomb
