#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ConfigFile.h"
#include "RunOptions.h"

namespace {

class TemporaryConfig {
  public:
	explicit TemporaryConfig(const std::string& contents) {
		path = std::filesystem::temp_directory_path() /
			   "negpar_inhomo_test_config.json";
		std::ofstream file(path);
		REQUIRE(file.good());
		file << contents;
	}

	~TemporaryConfig() {
		std::error_code error;
		std::filesystem::remove(path, error);
	}

	const std::filesystem::path& value() const { return path; }

  private:
	std::filesystem::path path;
};

coulomb::RunOptions parse(std::initializer_list<std::string> arguments) {
	std::vector<std::string> values(arguments);
	std::vector<char*> argv;
	argv.reserve(values.size());
	for (auto& value : values)
		argv.push_back(value.data());
	return coulomb::RunOptions::parse(static_cast<int>(argv.size()), argv.data());
}

} // namespace

TEST_CASE("configuration defaults preserve decoupled fixed behavior",
		  "[configuration]") {
	const auto options = parse({"negpar_inhomo"});
	REQUIRE(options.parameters.hdpCouplingMode ==
			coulomb::HdpCouplingMode::Decoupled);
	REQUIRE(options.parameters.effectiveWeightPolicy ==
			coulomb::EffectiveWeightPolicy::Fixed);
	REQUIRE_FALSE(options.parameters.weightedFourierResampling);
	REQUIRE_FALSE(options.parameters.partialResampling);
	REQUIRE(options.outputDirectory == "result");
	REQUIRE(options.effectiveConfigJson().find("\"schema_version\": 1") !=
			std::string::npos);
}

TEST_CASE("configuration file is loaded before command line overrides",
		  "[configuration]") {
	TemporaryConfig file(R"json({
  "schema_version": 1,
  "features": {
    "weighted_hdp": true,
    "weighted_fourier_resampling": true,
    "partial_resampling": true,
    "adaptive_effective_weights": true,
    "linearized_coulomb": true,
    "delta_m_sampling": false,
    "projection_mode": "maxwellian_only",
    "coulomb_bgk_hybrid": true
  },
  "resampling": {
    "partial_cutoff_standard_deviations": 2.5,
    "conserve_weighted_moments": true,
    "signed_weight_min": 0.000001,
    "signed_weight_max": 0.004,
    "full_weight_min": 0.000002,
    "full_weight_max": 0.003,
    "cpu_cost_constant": 0.2,
    "cpu_cost_collision_coefficient": 3.0
  },
  "collisions": {"bgk_strength": 0.5},
  "runtime": {
    "seed": 17,
    "steps": 8,
    "threads": 2,
    "output_directory": "config-result"
  }
})json");

	const auto options = parse({"negpar_inhomo", "--config",
							   file.value().string(), "--seed", "19",
							   "--steps", "4"});
	REQUIRE(options.parameters.hdpCouplingMode ==
			coulomb::HdpCouplingMode::VarianceWeighted);
	REQUIRE(options.parameters.effectiveWeightPolicy ==
			coulomb::EffectiveWeightPolicy::QuadraticAdaptive);
	REQUIRE(options.parameters.collisionCoupling ==
			coulomb::CollisionCoupling::Linearized);
	REQUIRE(options.parameters.projectionMode ==
			coulomb::ProjectionMode::MaxwellianOnly);
	REQUIRE(options.parameters.deltaMMode == coulomb::DeltaMMode::Disabled);
	REQUIRE(options.parameters.coulombBgkHybrid);
	REQUIRE(options.parameters.bgkStrength == 0.5);
	REQUIRE(options.parameters.partialResamplingCutoff == 2.5);
	REQUIRE(options.seed == 19);
	REQUIRE(options.steps == 4);
	REQUIRE(options.threads == 2);
	REQUIRE(options.outputDirectory ==
			(std::filesystem::path(file.value().parent_path()) /
			 "config-result")
				.lexically_normal()
				.string());
}

TEST_CASE("configuration rejects unknown keys and invalid dependencies",
		  "[configuration]") {
	TemporaryConfig unknown(R"json({
  "schema_version": 1,
  "features": {"weighted_hdp": false, "typo": true}
})json");
	REQUIRE_THROWS_AS(coulomb::ConfigFile::load(unknown.value()),
					  std::invalid_argument);

	TemporaryConfig dependency(R"json({
  "schema_version": 1,
  "features": {"partial_resampling": true}
})json");
	REQUIRE_THROWS_AS(
		coulomb::ConfigFile::load(dependency.value()),
		std::invalid_argument);
}

TEST_CASE("validate-config requires a file and does not start a run",
		  "[configuration]") {
	TemporaryConfig file(R"json({"schema_version": 1})json");
	const auto options =
		parse({"negpar_inhomo", "--validate-config", file.value().string()});
	REQUIRE(options.validateConfigOnly);
	REQUIRE(options.configPath.has_value());
	REQUIRE_THROWS(parse({"negpar_inhomo", "--validate-config"}));
}

TEST_CASE("configuration rejects unsupported schema and invalid types",
		  "[configuration][validation]") {
	TemporaryConfig unsupported(R"json({"schema_version": 2})json");
	REQUIRE_THROWS_AS(coulomb::ConfigFile::load(unsupported.value()),
					  std::invalid_argument);

	TemporaryConfig wrongType(
		R"json({"schema_version": 1, "features": {"weighted_hdp": 1}})json");
	REQUIRE_THROWS_AS(coulomb::ConfigFile::load(wrongType.value()),
					  std::invalid_argument);

	TemporaryConfig nonFinite(
		R"json({"schema_version": 1, "resampling": {"signed_weight_min": 1e999}})json");
	REQUIRE_THROWS_AS(coulomb::ConfigFile::load(nonFinite.value()),
					  std::invalid_argument);
}
