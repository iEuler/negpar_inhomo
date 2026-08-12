#pragma once

#include <initializer_list>
#include <vector>

#include "Particle.h"
#include "RandomContext.h"

namespace coulomb {

struct Moments {
	double m0{};
	double m11{}, m12{}, m13{};
	double m2{};
	double m21{}, m22{}, m23{};
	double m31{}, m32{}, m33{};
};

class ParticleGroup {
  public:
	Moments moments;
	double elecfield;

	ParticleGroup() : elecfield(0.), xmin(0.), xmax(0.) {}
	void set_xrange(double, double);
	double get_xmin() const { return xmin; }
	double get_xmax() const { return xmax; }
	int size() const { return static_cast<int>(vS.size()); }
	void push_back(const Particle1d3d&);
	void erase(int);
	std::vector<Particle1d3d>& list() { return vS; }
	const std::vector<Particle1d3d>& list() const { return vS; }
	Particle1d3d& list(int index);
	const Particle1d3d& list(int index) const;
	void computemoments();

  private:
	double xmin, xmax;
	std::vector<Particle1d3d> vS;
};

class NeParticleGroup {
  public:
	double rhoM{}, u1M{}, u2M{}, u3M{}, TprtM{};
	double T1M{}, T2M{}, T3M{};
	double rhoF{}, u1F{}, u2F{}, u3F{}, TprtF{};
	double dx_rhoM{}, dx_u1M{}, dx_TprtM{};
	double rho{}, u1{}, Tprt{};
	double rho_o{}, u1_o{}, Tprt_o{};
	double drho{}, dm1{}, denergy{};
	double drho_g{}, dm1_g{}, denergy_g{};
	Moments positive_moments;
	Moments negative_moments;
	Moments full_moments;
	Moments previous_positive_moments;
	Moments previous_negative_moments;
	double elecfield{}, elecfield_F{};
	double alpha_neg{}, alpha_pos{}, rmax{};
	std::vector<double> xyz_minmax{0., 0., 0., 0., 0., 0.};
	bool isResampled{false};

	NeParticleGroup() = default;
	void set_xrange(double, double);
	double get_xmin() const { return xmin; }
	double get_xmax() const { return xmax; }
	int size(ParticleKind kind) const;
	void push_back(const Particle1d3d&, ParticleKind kind);
	void erase(int index, ParticleKind kind);
	void clear(ParticleKind kind);
	std::vector<Particle1d3d>& list(ParticleKind kind);
	const std::vector<Particle1d3d>& list(ParticleKind kind) const;
	Particle1d3d& list(int index, ParticleKind kind);
	const Particle1d3d& list(int index, ParticleKind kind) const;
	void computemoments();
	void set_xyzrange();
	void set_xyzrange(std::initializer_list<ParticleKind> kinds);
	void copymoments();
	void reset_flag_resampled() { isResampled = false; }
	void setPositionRangeAndRandomizeValues(double xmin, double xmax,
											RandomContext& random);

  private:
	double xmin{0.}, xmax{0.};
	std::vector<Particle1d3d> vSp, vSn, vSf;
};

} // namespace coulomb
