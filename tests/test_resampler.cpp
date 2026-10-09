#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <limits>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Constants.h"
#include "BoundedMomentCorrection.h"
#include "EffectiveWeightSelector.h"
#include "ParticlePartition.h"
#include "RandomContext.h"
#include "RandomSampling.h"
#include "Resampler.h"
#include "ResamplerHelper.h"
#include "ResamplingNumerics.h"
#include "ResamplingVelocity.h"
#include "WeightedFourierCoupling.h"

TEST_CASE("Fourier variance vanishes for empty and deterministic populations", "[resampling][research]") {
	using coulomb::resampling::WeightedFourierCoupling;
	constexpr double piValue = 3.14159265358979323846;
	const double volume = 8.0 * piValue * piValue * piValue;
	const double weight = 0.01;
	REQUIRE(WeightedFourierCoupling::particleVariance({}, weight, 0) == 0.0);
	const std::complex<double> coefficient = std::polar(10.0 * weight / volume, 0.7);
	REQUIRE(WeightedFourierCoupling::particleVariance(coefficient, weight, 10) < 1e-16);
	const auto positive = std::polar(10.0 * weight / volume, 0.7);
	const auto negative = std::polar(10.0 * weight / volume, 1.4);
	REQUIRE(WeightedFourierCoupling::optimalWeight({}, positive, negative,
		weight, weight, 10, 10, 10) < 1e-12);
}

TEST_CASE("bounded allocation preserves variance constraints and reduces kinetic cost", "[resampling][adaptive][research]") {
	const std::vector<coulomb::EffectiveWeightCell> cells{{0.8, 1.0, 1000, 1000}};
	const auto result = coulomb::EffectiveWeightSelector{}.select(
		0.001, 0.001, 0.0001, 0.005, 0.0001, 0.005, 0.2, 3.0, 0.01, 1.0, cells);
	const double x = 0.001 / result.signedWeight;
	const double y = 0.001 / result.fullWeight;
	REQUIRE(0.2 * x + 0.8 * y >= 1.0 - 1e-10);
	REQUIRE(y >= x - 1e-10);
	REQUIRE(1.23 * x + 0.5 * y <= 1.73 + 1e-10);
	REQUIRE(result.fullWeight < 0.001);
	REQUIRE(result.signedWeight > 0.001);
}

TEST_CASE("Certified quadratic envelope covers interior extrema and cross terms", "[resampling][bounds]") {
	using coulomb::resampling::ResamplingNumerics;
	for (const auto& d : {std::vector<double>{.01,-.2,0,0,2,0,0,0,0,0},
		std::vector<double>{-.2,1,-2,3,-4,5,-6,7,-8,9}}) {
		const double bound=ResamplingNumerics{}.quadraticEnvelope(.5,d);
		for (int x=0;x<=10;++x) for (int y=0;y<=10;++y) for (int z=0;z<=10;++z)
			REQUIRE(std::abs(ResamplingNumerics{}.evaluateQuadraticTaylor(x*.1-.5,y*.1-.5,z*.1-.5,d)) <= bound);
	}
	REQUIRE_THROWS_AS(ResamplingNumerics{}.quadraticEnvelope(-1,std::vector<double>(10)),std::invalid_argument);
	auto invalid=std::vector<double>(10); invalid[7]=std::numeric_limits<double>::infinity();
	REQUIRE_THROWS_AS(ResamplingNumerics{}.quadraticEnvelope(.5,invalid),std::invalid_argument);
	invalid.assign(10,0.0); invalid[0]=std::numeric_limits<double>::max();
	REQUIRE_THROWS_AS(ResamplingNumerics{}.quadraticEnvelope(.5,invalid),std::invalid_argument);
}

TEST_CASE("Certified rejection recovers analytic quadratic mass in expectation", "[resampling][bounds][statistical]") {
	coulomb::RandomContext random; random.reseed(71003);
	coulomb::resampling::ResamplerHelper helper(random);
	coulomb::NeParticleGroup result;
	const std::vector<double> d{.01,-.2,0,0,2,0,0,0,0,0};
	const double bound=coulomb::resampling::ResamplingNumerics{}.quadraticEnvelope(.5,d);
	constexpr int count=40000;
	for (int i=0;i<count;++i) {
		double x=coulomb::RandomSampling(random).uniform()-.5;
		helper.acceptBoundedSample({coulomb::pi+x,coulomb::pi,coulomb::pi},result,
			coulomb::resampling::ResamplingNumerics{}.evaluateQuadraticTaylor(x,0,0,d),bound);
	}
	const double probability=(1./12.+.01)/bound;
	REQUIRE(std::abs(result.size(coulomb::ParticleKind::Positive)-count*probability)<6*std::sqrt(count*probability*(1-probability)));
	REQUIRE(result.size(coulomb::ParticleKind::Negative)==0);
	REQUIRE_THROWS_AS(helper.acceptBoundedSample({coulomb::pi,coulomb::pi,coulomb::pi},result,2*bound,bound),std::runtime_error);
}

TEST_CASE("Periodic cell wrapping restores the analytic sphere volume", "[resampling][geometry][statistical]") {
	using namespace coulomb;
	using namespace coulomb::resampling;
	const double h=pi/8.;
	RandomContext positions,oldRandom,newRandom;
	positions.reseed(810701); oldRandom.reseed(811); newRandom.reseed(811);
	ResamplerHelper oldHelper(oldRandom),newHelper(newRandom);
	NeParticleGroup pointOld,pointNew;
	oldHelper.acceptBoundedSample({-0.1,pi,pi},pointOld,1.,1.);
	newHelper.acceptBoundedSample(ResamplingNumerics{}.wrapPeriodicSample({-0.1,pi,pi}),pointNew,1.,1.);
	REQUIRE(pointOld.size(ParticleKind::Positive)==0);
	REQUIRE(pointNew.size(ParticleKind::Positive)==1);
	constexpr int count=100000;
	int oldCount=0,newCount=0;
	auto inSphere=[](const std::vector<double>& v) {
		double radiusSquared=0.;
		for(double x:v) radiusSquared+=(x-pi)*(x-pi);
		return radiusSquared<pi*pi;
	};
	for (int i=0;i<count;++i) {
		std::vector<double> v(3);
		for(double& x:v) x=2*pi*RandomSampling(positions).uniform()-h;
		oldCount+=inSphere(v);
		newCount+=inSphere(ResamplingNumerics{}.wrapPeriodicSample(v));
	}
	const double a=h/pi;
	// Three disjoint spherical caps: pi*R*h^2 - pi*h^3/3 each.
	const double missingFraction=pi*(3*a*a-a*a*a)/8.;
	REQUIRE(newCount>oldCount);
	REQUIRE(std::abs((newCount-oldCount)-count*missingFraction)<6*std::sqrt(count*missingFraction*(1-missingFraction)));
	const double sphereFraction=pi/6.;
	REQUIRE(std::abs(newCount-count*sphereFraction)<6*std::sqrt(count*sphereFraction*(1-sphereFraction)));
	REQUIRE(ResamplingNumerics{}.wrapPeriodicSample({-.1,2*pi,4*pi+.2})[0]==Catch::Approx(2*pi-.1));
	REQUIRE(ResamplingNumerics{}.wrapPeriodicSample({-.1,2*pi,4*pi+.2})[1]==0.);
	REQUIRE_THROWS_AS(ResamplingNumerics{}.wrapPeriodicSample({0.,0.,std::numeric_limits<double>::infinity()}),std::invalid_argument);
}

namespace { // Fourier resampler fixtures

coulomb::NeParticleGroup signedFixture() {
	coulomb::NeParticleGroup particles;
	for (int component = 0; component < 3; ++component) {
		for (const double sign : {-1.0, 1.0}) {
			std::vector<double> velocity(3, 0.0);
			velocity[component] = sign;
			const coulomb::Particle1D3D anchor(velocity);
			particles.pushBack(anchor, coulomb::ParticleKind::Positive);
			particles.pushBack(anchor, coulomb::ParticleKind::Negative);
		}
	}

	for (int index = 0; index < 24; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 24.0;
		particles.pushBack(
			coulomb::Particle1D3D({0.5 * std::cos(angle), 0.5 * std::sin(angle),
								   0.25 * std::cos(2.0 * angle)}),
			coulomb::ParticleKind::Positive);
	}
	for (int index = 0; index < 8; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 8.0;
		particles.pushBack(
			coulomb::Particle1D3D({0.2 * std::cos(angle), 0.2 * std::sin(angle),
								   0.1 * std::cos(2.0 * angle)}),
			coulomb::ParticleKind::Negative);
	}
	return particles;
}

void requireSameParticles(const coulomb::NeParticleGroup& first,
						  const coulomb::NeParticleGroup& second,
						  coulomb::ParticleKind kind) {
	REQUIRE(first.size(kind) == second.size(kind));
	for (int index = 0; index < first.size(kind); ++index) {
		for (int component = 0; component < 3; ++component) {
			REQUIRE(first.list(index, kind).velocity(component) ==
					second.list(index, kind).velocity(component));
		}
	}
}

TEST_CASE("Fixed physical resampling bounds retain a translated spherical domain", "[resampling][geometry][bounds]") {
	using namespace coulomb;
	using namespace coulomb::resampling;
	auto particles=signedFixture();
	const std::array<double,3> center{1.,-.5,.25};
	for(auto kind:{ParticleKind::Positive,ParticleKind::Negative})
		for(auto& p:particles.list(kind)) {
			auto v=p.velocity();
			for(int j=0;j<3;++j) v[j]+=center[j];
			p.setVelocity(v);
		}
	FourierResamplerConfig config;
	config.frequencyCount=4; config.effectiveParticleWeight=.1;
	config.envelope=ResamplingEnvelope::CertifiedQuadratic;
	config.cellGeometry=ResamplingCellGeometry::PeriodicWrapped;
	config.fixedVelocityBounds={-1.,3.,-2.5,1.5,-1.75,2.25};
	RandomContext firstRandom,secondRandom;
	firstRandom.reseed(811307); secondRandom.reseed(811307);
	const auto first=FourierResampler(particles,config).resample(firstRandom);
	const auto second=FourierResampler(particles,config).resample(secondRandom);
	REQUIRE(first.size(ParticleKind::Positive)+first.size(ParticleKind::Negative)>0);
	for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) {
		requireSameParticles(first,second,kind);
		for(const auto& p:first.list(kind)) {
			double r2=0.;
			for(int j=0;j<3;++j) r2+=(p.velocity(j)-center[j])*(p.velocity(j)-center[j]);
			REQUIRE(r2<4.+1e-12);
		}
	}
	REQUIRE(particles.xyzMinMax==std::vector<double>{0.,0.,0.,0.,0.,0.});
	for(const auto& bounds:{std::vector<double>{0.,1.},
		std::vector<double>{0.,0.,-1.,1.,-1.,1.},
		std::vector<double>{1.,-1.,-1.,1.,-1.,1.},
		std::vector<double>{0.,std::numeric_limits<double>::infinity(),-1.,1.,-1.,1.},
		std::vector<double>{-std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),-1.,1.,-1.,1.}}) {
		auto invalid=config; invalid.fixedVelocityBounds=bounds;
		REQUIRE_THROWS_AS(FourierResampler(particles,invalid),std::invalid_argument);
	}
	auto outside=particles;
	outside.pushBack(Particle1D3D({4.,0.,0.}),ParticleKind::Positive);
	REQUIRE_THROWS_AS(FourierResampler(outside,config).resample(firstRandom),std::invalid_argument);
	auto invalid=config; invalid.cellGeometry=ResamplingCellGeometry::LegacyShifted;
	REQUIRE_THROWS_AS(FourierResampler(particles,invalid),std::invalid_argument);
}

TEST_CASE("Aligned core bounds map the physical partition radius exactly", "[resampling][geometry]") {
	using namespace coulomb;
	NeParticleGroup source;
	const std::array<double,3> center{1.,-.5,.25};
	constexpr double radius=3.;
	source.xyzMinMax={-2.,4.,-3.5,2.5,-2.75,3.25};
	for(const auto& offset:{std::array<double,3>{2.6,1.2,0.},std::array<double,3>{-2.,1.,1.},std::array<double,3>{0.,0.,0.}})
		source.pushBack(Particle1D3D({center[0]+offset[0],center[1]+offset[1],center[2]+offset[2]}),ParticleKind::Positive);
	const auto normalized=resampling::ResamplingVelocity{}.normalizeSigned(source);
	for(int i=0;i<source.size(ParticleKind::Positive);++i) {
		double physical=0.,scaled=0.;
		for(int j=0;j<3;++j) {
			physical+=std::pow(source.list(i,ParticleKind::Positive).velocity(j)-center[j],2);
			scaled+=std::pow(normalized.list(i,ParticleKind::Positive).velocity(j)-pi,2);
		}
		REQUIRE(scaled==Catch::Approx(physical*pi*pi/(radius*radius)).margin(1e-14));
		REQUIRE(scaled<pi*pi);
	}
}

TEST_CASE("negpar.unit.resampling.weighted Fourier reconstruction replays",
		  "[resampling][fourier][weighted][reproducibility]") {
	auto particles = signedFixture();
	for (int index = 0; index < 12; ++index) {
		const double angle = 2.0 * coulomb::pi * index / 12.0;
		particles.pushBack(
			coulomb::Particle1D3D({0.35 * std::cos(angle),
								   0.35 * std::sin(angle), 0.1}),
			coulomb::ParticleKind::Full);
	}
	coulomb::resampling::FourierResamplerConfig config;
	config.effectiveParticleWeight = 0.1;
	config.fullParticleWeight = 0.2;
	config.weightedCoupling = true;
	config.frequencyCount = 4;
	config.maxSamplingAttempts = 10'000;

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(440011);
	secondRandom.reseed(440011);
	const auto originalBounds = particles.xyzMinMax;
	const auto first =
		coulomb::resampling::FourierResampler(particles, config).resample(
			firstRandom);
	const auto second =
		coulomb::resampling::FourierResampler(particles, config).resample(
			secondRandom);
	requireSameParticles(first, second, coulomb::ParticleKind::Positive);
	requireSameParticles(first, second, coulomb::ParticleKind::Negative);
	REQUIRE(particles.xyzMinMax == originalBounds);
	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		for (const auto& particle : first.list(kind))
			for (int component = 0; component < 3; ++component)
				REQUIRE(std::isfinite(particle.velocity(component)));
	}
}

} // namespace

TEST_CASE("negpar.unit.resampling.Fourier resampler uses explicit RNG",
		  "[resampling][fourier]") {
	static_assert(
		std::is_constructible_v<coulomb::resampling::FourierResampler,
								const coulomb::NeParticleGroup&,
								coulomb::resampling::FourierResamplerConfig>);

	coulomb::NeParticleGroup first;
	coulomb::NeParticleGroup second;
	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(9876);
	secondRandom.reseed(9876);

	double firstBound = 1.0;
	double secondBound = 1.0;
	const std::vector<double> sample{coulomb::pi, coulomb::pi, coulomb::pi};

	coulomb::resampling::ResamplerHelper(firstRandom)
		.acceptSample(sample, first, 0.25, firstBound);
	coulomb::resampling::ResamplerHelper(secondRandom)
		.acceptSample(sample, second, 0.25, secondBound);

	REQUIRE(firstBound == secondBound);
	REQUIRE(first.size(coulomb::ParticleKind::Positive) ==
			second.size(coulomb::ParticleKind::Positive));
	REQUIRE(first.size(coulomb::ParticleKind::Negative) ==
			second.size(coulomb::ParticleKind::Negative));
}

TEST_CASE("negpar.unit.resampling.Fourier resampler configuration rejects "
		  "invalid grids",
		  "[resampling][fourier][validation]") {
	coulomb::NeParticleGroup particles;
	coulomb::resampling::FourierResamplerConfig config;

	config.effectiveParticleWeight = 0.0;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);

	config.effectiveParticleWeight = 1.0;
	config.frequencyCount = 3;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);

	config.frequencyCount = 4;
	config.maxSamplingAttempts = 0;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),
					  std::invalid_argument);
	config.maxSamplingAttempts = 1;
	config.useApproximation = false;
	config.envelope = coulomb::resampling::ResamplingEnvelope::CertifiedQuadratic;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler(particles, config),std::invalid_argument);
}

TEST_CASE("negpar.unit.resampling.weighted Fourier coupling clamps finite "
		  "frequency weights",
		  "[resampling][fourier][weighted]") {
	const std::complex<double> full{0.1, 0.2};
	const std::complex<double> positive{0.4, -0.1};
	const std::complex<double> negative{-0.2, 0.3};
	const double weight = coulomb::resampling::WeightedFourierCoupling::
		optimalWeight(full, positive, negative, 0.2, 0.1, 12, 20, 8);
	REQUIRE(std::isfinite(weight));
	REQUIRE(weight >= 0.0);
	REQUIRE(weight <= 1.0);
	const auto blended = coulomb::resampling::WeightedFourierCoupling::blend(
		full, positive - negative, weight);
	REQUIRE(std::isfinite(blended.real()));
	REQUIRE(std::isfinite(blended.imag()));
	REQUIRE(coulomb::resampling::WeightedFourierCoupling::blend(
				full, positive - negative, 0.0) == positive - negative);
	REQUIRE(coulomb::resampling::WeightedFourierCoupling::blend(
				full, positive - negative, 1.0) == full);
}

TEST_CASE("negpar.unit.resampling.partial partition keeps typed core and tail",
		  "[resampling][partial]") {
	coulomb::NeParticleGroup source;
	source.u1M = 0.0;
	source.u2M = 0.0;
	source.u3M = 0.0;
	source.tprtM = 1.0;
	source.pushBack(coulomb::Particle1D3D(0.25, {0.5, 0.0, 0.0}),
				   coulomb::ParticleKind::Positive);
	source.pushBack(coulomb::Particle1D3D(0.75, {4.0, 0.0, 0.0}),
				   coulomb::ParticleKind::Positive);
	source.pushBack(coulomb::Particle1D3D(0.5, {0.0, 0.0, 0.0}),
				   coulomb::ParticleKind::Full);

	const auto partition = coulomb::ParticlePartitioning::split(source, 3.0);
	REQUIRE(partition.core.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(partition.tail.size(coulomb::ParticleKind::Positive) == 1);
	REQUIRE(partition.core.size(coulomb::ParticleKind::Full) == 1);
	REQUIRE(partition.tail.size(coulomb::ParticleKind::Full) == 0);
	REQUIRE(partition.tail.list(0, coulomb::ParticleKind::Positive).position() ==
			Catch::Approx(0.75));
}

TEST_CASE("negpar.unit.resampling.adaptive effective weights stay bounded "
		  "for degenerate cells",
		  "[resampling][adaptive]") {
	const std::vector<coulomb::EffectiveWeightCell> cells{
		{0.0, 1.0, 0, 0}, {1.0, 1.0, 20, 1}};
	const auto selection = coulomb::EffectiveWeightSelector{}.select(
		0.001, 0.002, 0.0001, 0.005, 0.0001, 0.005, 0.205, 3.277, 0.01,
		10.0, cells);
	REQUIRE(std::isfinite(selection.signedWeight));
	REQUIRE(std::isfinite(selection.fullWeight));
	REQUIRE(selection.signedWeight >= 0.0001);
	REQUIRE(selection.signedWeight <= 0.005);
	REQUIRE(selection.fullWeight >= 0.0001);
	REQUIRE(selection.fullWeight <= 0.005);
}

TEST_CASE(
	"negpar.unit.resampling.Fourier interpolation uses every first derivative",
	"[resampling][fourier][interpolation]") {
	const std::vector<double> derivatives{1.0, 2.0, 3.0, 5.0, 0.0,
										  0.0, 0.0, 0.0, 0.0, 0.0};
	REQUIRE(coulomb::resampling::ResamplingNumerics{}.evaluateQuadraticTaylor(
				0.1, 0.2, 0.4, derivatives) == Catch::Approx(3.8));
}

TEST_CASE("negpar.unit.resampling.Fourier upper bounds cover their own "
		  "interpolation cell",
		  "[resampling][fourier][bounds]") {
	coulomb::Vector3D values(2, std::vector(2, std::vector<double>(2, 0.0)));
	values[0][0][0] = 7.0;

	coulomb::RandomContext random;
	const auto bounds =
		coulomb::resampling::ResamplerHelper(random).upperBound(values);
	for (const auto& plane : bounds)
		for (const auto& row : plane)
			for (const double bound : row)
				REQUIRE(bound == Catch::Approx(7.0));
}

TEST_CASE(
	"negpar.unit.resampling.Fourier signed resampling conserves low moments",
	"[resampling][fourier][conservation]") {
	auto particles = signedFixture();
	auto original = particles;
	original.computeMoments();

	coulomb::resampling::FourierResamplerConfig config;
	config.effectiveParticleWeight = 0.1;
	config.frequencyCount = 4;
	config.maxSamplingAttempts = 10'000;

	SECTION("exact Fourier transform") { config.useApproximation = false; }
	SECTION("approximate Fourier transform") { config.useApproximation = true; }
	SECTION("certified quadratic envelope") {
		config.useApproximation = true;
		config.envelope = coulomb::resampling::ResamplingEnvelope::CertifiedQuadratic;
	}
	SECTION("certified stratified proposal allocation") {
		config.envelope = coulomb::resampling::ResamplingEnvelope::CertifiedQuadratic;
		config.cellGeometry = coulomb::resampling::ResamplingCellGeometry::PeriodicWrapped;
		config.proposalAllocation = coulomb::resampling::ResamplingProposalAllocation::Stratified;
	}
	SECTION("certified periodic cell geometry") {
		config.useApproximation = true;
		config.envelope = coulomb::resampling::ResamplingEnvelope::CertifiedQuadratic;
		config.cellGeometry = coulomb::resampling::ResamplingCellGeometry::PeriodicWrapped;
	}

	coulomb::RandomContext firstRandom;
	coulomb::RandomContext secondRandom;
	firstRandom.reseed(20260808);
	secondRandom.reseed(20260808);

	const coulomb::resampling::FourierResampler firstResampler(particles,
															   config);
	const coulomb::resampling::FourierResampler secondResampler(particles,
																config);
	auto first = firstResampler.resample(firstRandom);
	const auto second = secondResampler.resample(secondRandom);

	requireSameParticles(first, second, coulomb::ParticleKind::Positive);
	requireSameParticles(first, second, coulomb::ParticleKind::Negative);
	REQUIRE(particles.xyzMinMax ==
			std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

	first.computeMoments();
	const double originalMass =
		config.effectiveParticleWeight *
		(original.positiveMoments.m0 - original.negativeMoments.m0);
	const double sampledMass =
		config.effectiveParticleWeight *
		(first.positiveMoments.m0 - first.negativeMoments.m0);
	const double originalEnergy =
		config.effectiveParticleWeight *
		(original.positiveMoments.m2 - original.negativeMoments.m2);
	const double sampledEnergy =
		config.effectiveParticleWeight *
		(first.positiveMoments.m2 - first.negativeMoments.m2);

	CAPTURE(originalMass, sampledMass, originalEnergy, sampledEnergy,
			first.positiveMoments.m0, first.negativeMoments.m0);
	REQUIRE(sampledMass == Catch::Approx(originalMass).margin(0.8));
	REQUIRE(sampledEnergy == Catch::Approx(originalEnergy).margin(1.0));

	for (const auto kind :
		 {coulomb::ParticleKind::Positive, coulomb::ParticleKind::Negative}) {
		for (const auto& particle : first.list(kind)) {
			for (int component = 0; component < 3; ++component) {
				REQUIRE(std::isfinite(particle.velocity(component)));
				REQUIRE(std::abs(particle.velocity(component)) <= 1.000001);
			}
		}
	}
}

TEST_CASE("Stratified proposal allocation balances prefixes and preserves cell expectations", "[resampling][stratified][statistical]") {
	using coulomb::resampling::StratifiedProposalAllocator;
	const std::array<double,7> expected{{.2,.4,1.8,0.,2.1,8.,.5}};
	std::array<double,7> sums{},squares{};
	coulomb::RandomContext random; random.reseed(831006);
	bool balanced=true, nonnegative=true;
	constexpr int replicas=40000;
	for(int r=0;r<replicas;++r) {
		StratifiedProposalAllocator allocation;
		double prefix=0.; int total=0;
		for(std::size_t j=0;j<expected.size();++j) {
			const int n=allocation.allocate(expected[j],random);
			prefix+=expected[j]; total+=n;
			balanced=balanced && total>=std::floor(prefix) && total<=std::ceil(prefix);
			nonnegative=nonnegative && n>=0;
			sums[j]+=n; squares[j]+=static_cast<double>(n)*n;
		}
	}
	REQUIRE(balanced); REQUIRE(nonnegative);
	for(std::size_t j=0;j<expected.size();++j) {
		const double mean=sums[j]/replicas;
		const double variance=(squares[j]-sums[j]*mean)/(replicas-1);
		CAPTURE(j,mean,expected[j]);
		REQUIRE(std::abs(mean-expected[j])<=5.*std::sqrt(variance/replicas)+1e-12);
	}
	StratifiedProposalAllocator invalid;
	REQUIRE_THROWS_AS(invalid.allocate(-1.,random),std::invalid_argument);
	REQUIRE_THROWS_AS(invalid.allocate(std::numeric_limits<double>::infinity(),random),std::invalid_argument);
	REQUIRE_THROWS_AS(invalid.allocate(std::numeric_limits<int>::max(),random),std::invalid_argument);
	coulomb::resampling::FourierResamplerConfig config;
	config.proposalAllocation=coulomb::resampling::ResamplingProposalAllocation::Stratified;
	REQUIRE_THROWS_AS(coulomb::resampling::FourierResampler({},config),std::invalid_argument);
}

TEST_CASE("Bounded moment correction preserves signed moments and translated support", "[resampling][correction]") {
    using namespace coulomb;
    using namespace coulomb::resampling;
    const std::array<double,3> center{1.,-.5,.25};
    auto source=signedFixture();
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) for(auto& p:source.list(kind)) {
        auto v=p.velocity(); for(std::size_t j=0;j<3;++j) v[j]+=center[j]; p.setVelocity(v);
    }
    const double weight=.05;
    const auto target=BoundedMomentCorrection::moments(source,weight);
    auto candidate=source;
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) for(auto& p:candidate.list(kind)) {
        auto v=p.velocity(); for(std::size_t j=0;j<3;++j) v[j]=center[j]+.95*(v[j]-center[j])+.01; p.setVelocity(v);
    }
    candidate.pushBack(Particle1D3D({1.,-.5,.25}),ParticleKind::Positive);
    candidate.pushBack(Particle1D3D({12.,13.,14.}),ParticleKind::Full);
    candidate.rhoM=7.;
    auto second=candidate;
    RandomContext random; random.reseed(831007);
    const auto result=BoundedMomentCorrection{}.apply(candidate,target,weight,center,3.,random);
    REQUIRE(result.success()); REQUIRE(result.removed==1);
    const auto actual=BoundedMomentCorrection::moments(candidate,weight);
    for(std::size_t j=0;j<actual.size();++j) REQUIRE(actual[j]==Catch::Approx(target[j]).margin(1e-8));
    for(auto kind:{ParticleKind::Positive,ParticleKind::Negative}) for(const auto& p:candidate.list(kind)) {
        double r2=0.; for(std::size_t j=0;j<3;++j) r2+=std::pow(p.velocity(static_cast<int>(j))-center[j],2);
        REQUIRE(r2<=9.+1e-12);
    }
    REQUIRE(candidate.rhoM==7.); requireSameParticles(candidate,second,ParticleKind::Full);
    random.reseed(831007);
    REQUIRE(BoundedMomentCorrection{}.apply(second,target,weight,center,3.,random).success());
    requireSameParticles(candidate,second,ParticleKind::Positive);
    requireSameParticles(candidate,second,ParticleKind::Negative);
}

TEST_CASE("Bounded moment correction rejects incompatible or infeasible targets transactionally", "[resampling][correction]") {
    using namespace coulomb;
    using namespace coulomb::resampling;
    auto original=signedFixture(), candidate=original;
    RandomContext random; random.reseed(831008);
    const auto target=BoundedMomentCorrection::moments(original,.05);
    auto bad=target; bad[0]+=.025;
    REQUIRE(BoundedMomentCorrection{}.apply(candidate,bad,.05,{0.,0.,0.},2.,random).status==MomentCorrectionStatus::IncompatibleMass);
    requireSameParticles(candidate,original,ParticleKind::Positive);
    requireSameParticles(candidate,original,ParticleKind::Negative);
    bad=target; bad[1]=1e6;
    REQUIRE_FALSE(BoundedMomentCorrection{}.apply(candidate,bad,.05,{0.,0.,0.},2.,random).success());
    requireSameParticles(candidate,original,ParticleKind::Positive);
    requireSameParticles(candidate,original,ParticleKind::Negative);
    NeParticleGroup zero;
    zero.pushBack(Particle1D3D({0.,0.,0.}),ParticleKind::Positive);
    zero.pushBack(Particle1D3D({0.,0.,0.}),ParticleKind::Negative);
    SignedLowMoments nonzero{}; nonzero[4]=.01;
    REQUIRE(BoundedMomentCorrection{}.apply(zero,nonzero,.05,{0.,0.,0.},1.,random).status==MomentCorrectionStatus::Singular);
    REQUIRE(zero.list(0,ParticleKind::Positive).velocity(0)==0.);
    REQUIRE_THROWS_AS(BoundedMomentCorrection{}.apply(candidate,target,.05,{0.,0.,0.},0.,random),std::invalid_argument);
    REQUIRE_THROWS_AS(BoundedMomentCorrection{}.apply(candidate,target,.05,{0.,0.,0.},.1,random),std::invalid_argument);
    bad=target; bad[1]+=.001;
    MomentCorrectionConfig restricted; restricted.maxRmsDisplacementFraction=1e-12;
    REQUIRE(BoundedMomentCorrection{}.apply(candidate,bad,.05,{0.,0.,0.},2.,random,restricted).status==MomentCorrectionStatus::DisplacementLimit);
    requireSameParticles(candidate,original,ParticleKind::Positive);
    requireSameParticles(candidate,original,ParticleKind::Negative);
}
