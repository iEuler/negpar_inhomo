#include "ParticlePartition.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ParticleGroup.h"

namespace coulomb {

ParticlePartition ParticlePartitioning::split(
	const NeParticleGroup& source, double standardDeviationCutoff) {
	if (!(standardDeviationCutoff > 0.0) ||
		!std::isfinite(standardDeviationCutoff))
		throw std::invalid_argument(
			"partial resampling cutoff must be finite and positive");
	if (!std::isfinite(source.u1M) || !std::isfinite(source.u2M) ||
		!std::isfinite(source.u3M) || !std::isfinite(source.tprtM))
		throw std::invalid_argument(
			"partial resampling requires finite Maxwellian moments");

	ParticlePartition result{source, source};
	for (const auto kind : {ParticleKind::Positive, ParticleKind::Negative,
							ParticleKind::Full}) {
		result.core.clear(kind);
		result.tail.clear(kind);
		for (const auto& particle : source.list(kind)) {
			const double dv1 = particle.velocity(0) - source.u1M;
			const double dv2 = particle.velocity(1) - source.u2M;
			const double dv3 = particle.velocity(2) - source.u3M;
			const double temperature = std::max(source.tprtM, 1e-14);
			const double normalizedSquaredDistance =
				(dv1 * dv1 + dv2 * dv2 + dv3 * dv3) / temperature;
			if (!std::isfinite(normalizedSquaredDistance))
				throw std::invalid_argument(
					"partial resampling requires finite particle velocities");
			const bool inCore =
				std::sqrt(normalizedSquaredDistance) <= standardDeviationCutoff;
			(inCore ? result.core : result.tail).pushBack(particle, kind);
		}
	}
	return result;
}

} // namespace coulomb
