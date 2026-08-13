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
	double elecField;

	ParticleGroup() : elecField(0.), xmin(0.), xmax(0.) {}
	void setXRange(double, double);
	double getXMin() const { return xmin; }
	double getXMax() const { return xmax; }
	int size() const { return static_cast<int>(vS.size()); }
	void pushBack(const Particle1D3D&);
	void erase(int);
	std::vector<Particle1D3D>& list() { return vS; }
	const std::vector<Particle1D3D>& list() const { return vS; }
	Particle1D3D& list(int index);
	const Particle1D3D& list(int index) const;
	void computeMoments();

  private:
	double xmin, xmax;
	std::vector<Particle1D3D> vS;
};

class NeParticleGroup {
  public:
	double rhoM{}, u1M{}, u2M{}, u3M{}, tprtM{};
	double t1M{}, t2M{}, t3M{};
	double rhoF{}, u1F{}, u2F{}, u3F{}, tprtF{};
	double dxRhoM{}, dxU1M{}, dxTprtM{};
	double rho{}, u1{}, tprt{};
	double rhoO{}, u1O{}, tprtO{};
	double drho{}, dm1{}, dEnergy{};
	double drhoG{}, dm1G{}, dEnergyG{};
	Moments positiveMoments;
	Moments negativeMoments;
	Moments fullMoments;
	Moments previousPositiveMoments;
	Moments previousNegativeMoments;
	double elecField{}, elecFieldF{};
	double alphaNeg{}, alphaPos{}, rMax{};
	std::vector<double> xyzMinMax{0., 0., 0., 0., 0., 0.};
	bool isResampled{false};

	NeParticleGroup() = default;
	void setXRange(double, double);
	double getXMin() const { return xmin; }
	double getXMax() const { return xmax; }
	int size(ParticleKind kind) const;
	void pushBack(const Particle1D3D&, ParticleKind kind);
	void erase(int index, ParticleKind kind);
	void clear(ParticleKind kind);
	std::vector<Particle1D3D>& list(ParticleKind kind);
	const std::vector<Particle1D3D>& list(ParticleKind kind) const;
	Particle1D3D& list(int index, ParticleKind kind);
	const Particle1D3D& list(int index, ParticleKind kind) const;
	void computeMoments();
	void setXyzRange();
	void setXyzRange(std::initializer_list<ParticleKind> kinds);
	void copyMoments();
	void resetFlagResampled() { isResampled = false; }
	void setPositionRangeAndRandomizeValues(double xmin, double xmax,
											RandomContext& random);

  private:
	double xmin{0.}, xmax{0.};
	std::vector<Particle1D3D> vSp, vSn, vSf;
};

} // namespace coulomb
