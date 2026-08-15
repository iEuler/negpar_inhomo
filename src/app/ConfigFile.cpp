#include "ConfigFile.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace coulomb {
namespace {

enum class JsonType { Object, Array, String, Number, Boolean, Null };

struct JsonValue {
	JsonType type{JsonType::Null};
	std::map<std::string, JsonValue> object;
	std::vector<JsonValue> array;
	std::string string;
	double number{};
	bool boolean{};
};

class JsonParser {
  public:
	explicit JsonParser(std::string input) : input(std::move(input)) {}

	JsonValue parse() {
		skipWhitespace();
		auto value = parseValue("$");
		skipWhitespace();
		if (position != input.size())
			fail("$", "trailing characters");
		return value;
	}

  private:
	[[noreturn]] void fail(const std::string& path,
						   const std::string& message) const {
		throw std::invalid_argument("Invalid configuration at " + path + ": " +
									message);
	}

	void skipWhitespace() {
		while (position < input.size()) {
			const char character = input[position];
			if (character != ' ' && character != '\t' && character != '\n' &&
				character != '\r')
				break;
			++position;
		}
	}

	void expect(char expected, const std::string& path) {
		if (position >= input.size() || input[position] != expected)
			fail(path, std::string("expected '") + expected + "'");
		++position;
	}

	bool consume(char expected) {
		if (position < input.size() && input[position] == expected) {
			++position;
			return true;
		}
		return false;
	}

	JsonValue parseValue(const std::string& path) {
		skipWhitespace();
		if (position >= input.size())
			fail(path, "expected a value");
		switch (input[position]) {
		case '{':
			return parseObject(path);
		case '[':
			return parseArray(path);
		case '"':
			return JsonValue{JsonType::String, {}, {}, parseString(path), 0.0,
							 false};
		case 't':
			parseLiteral("true", path);
			return JsonValue{JsonType::Boolean, {}, {}, {}, 0.0, true};
		case 'f':
			parseLiteral("false", path);
			return JsonValue{JsonType::Boolean, {}, {}, {}, 0.0, false};
		case 'n':
			parseLiteral("null", path);
			return JsonValue{};
		default:
			if (input[position] == '-' ||
				(input[position] >= '0' && input[position] <= '9'))
				return parseNumber(path);
			fail(path, "unexpected character");
		}
	}

	JsonValue parseObject(const std::string& path) {
		JsonValue result;
		result.type = JsonType::Object;
		expect('{', path);
		skipWhitespace();
		if (consume('}'))
			return result;

		while (true) {
			skipWhitespace();
			if (position >= input.size() || input[position] != '"')
				fail(path, "object keys must be strings");
			const auto key = parseString(path + ".<key>");
			skipWhitespace();
			expect(':', path + "." + key);
			auto value = parseValue(path + "." + key);
			if (!result.object.emplace(key, std::move(value)).second)
				fail(path + "." + key, "duplicate key");
			skipWhitespace();
			if (consume('}'))
				return result;
			expect(',', path);
		}
	}

	JsonValue parseArray(const std::string& path) {
		JsonValue result;
		result.type = JsonType::Array;
		expect('[', path);
		skipWhitespace();
		if (consume(']'))
			return result;
		while (true) {
			result.array.push_back(parseValue(path + "[]"));
			skipWhitespace();
			if (consume(']'))
				return result;
			expect(',', path);
		}
	}

	std::string parseString(const std::string& path) {
		expect('"', path);
		std::string result;
		while (position < input.size()) {
			const char character = input[position++];
			if (character == '"')
				return result;
			if (static_cast<unsigned char>(character) < 0x20)
				fail(path, "control character in string");
			if (character != '\\') {
				result += character;
				continue;
			}
			if (position >= input.size())
				fail(path, "unterminated escape");
			const char escaped = input[position++];
			switch (escaped) {
			case '"':
			case '\\':
			case '/':
				result += escaped;
				break;
			case 'b':
				result += '\b';
				break;
			case 'f':
				result += '\f';
				break;
			case 'n':
				result += '\n';
				break;
			case 'r':
				result += '\r';
				break;
			case 't':
				result += '\t';
				break;
			case 'u': {
				if (position + 4 > input.size())
					fail(path, "incomplete unicode escape");
				unsigned value = 0;
				for (int index = 0; index < 4; ++index) {
					const char hex = input[position++];
					value <<= 4;
					if (hex >= '0' && hex <= '9')
						value += static_cast<unsigned>(hex - '0');
					else if (hex >= 'a' && hex <= 'f')
						value += static_cast<unsigned>(hex - 'a' + 10);
					else if (hex >= 'A' && hex <= 'F')
						value += static_cast<unsigned>(hex - 'A' + 10);
					else
						fail(path, "invalid unicode escape");
				}
				if (value > 0x7f)
					fail(path, "configuration strings must use ASCII escapes");
				result += static_cast<char>(value);
				break;
			}
			default:
				fail(path, "invalid string escape");
			}
		}
		fail(path, "unterminated string");
	}

	void parseLiteral(const char* literal, const std::string& path) {
		const std::string value(literal);
		if (input.compare(position, value.size(), value) != 0)
			fail(path, "invalid literal");
		position += value.size();
	}

	JsonValue parseNumber(const std::string& path) {
		const auto start = position;
		consume('-');
		if (position >= input.size())
			fail(path, "invalid number");
		if (input[position] == '0') {
			++position;
		} else {
			if (input[position] < '1' || input[position] > '9')
				fail(path, "invalid number");
			while (position < input.size() && input[position] >= '0' &&
				   input[position] <= '9')
				++position;
		}
		if (consume('.')) {
			if (position >= input.size() || input[position] < '0' ||
				input[position] > '9')
				fail(path, "invalid number fraction");
			while (position < input.size() && input[position] >= '0' &&
				   input[position] <= '9')
				++position;
		}
		if (position < input.size() &&
			(input[position] == 'e' || input[position] == 'E')) {
			++position;
			if (position < input.size() &&
				(input[position] == '+' || input[position] == '-'))
				++position;
			if (position >= input.size() || input[position] < '0' ||
				input[position] > '9')
				fail(path, "invalid number exponent");
			while (position < input.size() && input[position] >= '0' &&
				   input[position] <= '9')
				++position;
		}
		const auto text = input.substr(start, position - start);
		std::size_t parsed = 0;
		double number{};
		try {
			number = std::stod(text, &parsed);
		}
		catch (const std::exception&) {
			fail(path, "invalid number");
		}
		if (parsed != text.size() || !std::isfinite(number))
			fail(path, "number must be finite");
		return JsonValue{JsonType::Number, {}, {}, {}, number, false};
	}

	std::string input;
	std::size_t position{};
};

const JsonValue& requireObject(const JsonValue& value, const std::string& path) {
	if (value.type != JsonType::Object)
		throw std::invalid_argument("Invalid configuration at " + path +
									": expected an object");
	return value;
}

const JsonValue* optionalMember(const JsonValue& object, const std::string& key,
								const std::string& path) {
	const auto& members = requireObject(object, path).object;
	const auto iterator = members.find(key);
	return iterator == members.end() ? nullptr : &iterator->second;
}

void rejectUnknown(const JsonValue& object,
				   const std::string& path,
				   std::initializer_list<const char*> allowed) {
	const auto& members = requireObject(object, path).object;
	for (const auto& member : members) {
		bool known = false;
		for (const auto* name : allowed)
			known = known || member.first == name;
		if (!known)
			throw std::invalid_argument("Invalid configuration at " + path +
										": unknown key '" + member.first + "'");
	}
}

bool readBool(const JsonValue& value, const std::string& path) {
	if (value.type != JsonType::Boolean)
		throw std::invalid_argument("Invalid configuration at " + path +
									": expected a boolean");
	return value.boolean;
}

double readNumber(const JsonValue& value, const std::string& path) {
	if (value.type != JsonType::Number || !std::isfinite(value.number))
		throw std::invalid_argument("Invalid configuration at " + path +
									": expected a finite number");
	return value.number;
}

std::string readString(const JsonValue& value, const std::string& path) {
	if (value.type != JsonType::String)
		throw std::invalid_argument("Invalid configuration at " + path +
									": expected a string");
	return value.string;
}

std::uint32_t readUint32(const JsonValue& value, const std::string& path) {
	const double number = readNumber(value, path);
	if (number < 0.0 || number > std::numeric_limits<std::uint32_t>::max() ||
		std::floor(number) != number)
		throw std::invalid_argument("Invalid configuration at " + path +
									": expected a uint32 value");
	return static_cast<std::uint32_t>(number);
}

int readPositiveInt(const JsonValue& value, const std::string& path) {
	const double number = readNumber(value, path);
	if (number < 1.0 || number > std::numeric_limits<int>::max() ||
		std::floor(number) != number)
		throw std::invalid_argument("Invalid configuration at " + path +
									": expected a positive integer");
	return static_cast<int>(number);
}

void validateParameters(const ParaClass& parameters) {
	const auto invalid = [](double value) {
		return !std::isfinite(value);
	};
	if (parameters.signedWeightMin <= 0.0 ||
		parameters.signedWeightMin > parameters.signedWeightMax ||
		parameters.fullWeightMin <= 0.0 ||
		parameters.fullWeightMin > parameters.fullWeightMax ||
		invalid(parameters.signedWeightMax) || invalid(parameters.fullWeightMax) ||
		invalid(parameters.partialResamplingCutoff) ||
		parameters.partialResamplingCutoff <= 0.0 ||
		invalid(parameters.cpuCostConstant) ||
		invalid(parameters.cpuCostCollisionCoefficient) ||
		parameters.cpuCostConstant < 0.0 ||
		parameters.cpuCostCollisionCoefficient < 0.0 ||
		invalid(parameters.bgkStrength) || parameters.bgkStrength < 0.0 ||
		parameters.bgkStrength > 1.0) {
		throw std::invalid_argument(
			"Invalid configuration at resampling/collisions: numeric values are "
			"outside their allowed ranges");
	}
	if ((parameters.partialResampling ||
		 parameters.effectiveWeightPolicy ==
			 EffectiveWeightPolicy::QuadraticAdaptive) &&
		!parameters.weightedFourierResampling) {
		throw std::invalid_argument(
			"Invalid configuration: partial or adaptive resampling requires "
			"features.weighted_fourier_resampling");
	}
	if (parameters.coulombBgkHybrid && parameters.bgkStrength <= 0.0)
		throw std::invalid_argument(
			"Invalid configuration: coulomb_bgk_hybrid requires a positive "
			"collisions.bgk_strength");
}

void appendEscaped(std::ostringstream& output, const std::string& value) {
	output << '"';
	for (const char character : value) {
		switch (character) {
		case '"':
			output << "\\\"";
			break;
		case '\\':
			output << "\\\\";
			break;
		case '\n':
			output << "\\n";
			break;
		case '\r':
			output << "\\r";
			break;
		case '\t':
			output << "\\t";
			break;
		default:
			output << character;
			break;
		}
	}
	output << '"';
}

void appendBool(std::ostringstream& output, bool value) {
	output << (value ? "true" : "false");
}

void appendOptionalUint(std::ostringstream& output,
						const std::optional<std::uint32_t>& value) {
	if (value)
		output << *value;
	else
		output << "null";
}

void appendOptionalInt(std::ostringstream& output,
					   const std::optional<int>& value) {
	if (value)
		output << *value;
	else
		output << "null";
}

void appendOptionalString(std::ostringstream& output,
						  const std::optional<std::string>& value) {
	if (value)
		appendEscaped(output, *value);
	else
		output << "null";
}

} // namespace

ConfigFileValues ConfigFile::load(const std::filesystem::path& path) {
	std::ifstream file(path);
	if (!file)
		throw std::invalid_argument("Unable to open configuration file '" +
									path.string() + "'");
	std::ostringstream contents;
	contents << file.rdbuf();

	const auto root = JsonParser(contents.str()).parse();
	rejectUnknown(root, "$", {"schema_version", "features", "resampling",
							  "collisions", "runtime"});
	const auto* schema = optionalMember(root, "schema_version", "$");
	if (schema == nullptr || readNumber(*schema, "$.schema_version") != 1.0)
		throw std::invalid_argument(
			"Invalid configuration at $.schema_version: expected version 1");

	ConfigFileValues result;
	if (const auto* features = optionalMember(root, "features", "$")) {
		rejectUnknown(*features, "$.features",
					  {"weighted_hdp", "weighted_fourier_resampling",
					   "partial_resampling", "adaptive_effective_weights",
					   "linearized_coulomb", "delta_m_sampling",
					   "projection_mode", "coulomb_bgk_hybrid"});
		if (const auto* value =
				optionalMember(*features, "weighted_hdp", "$.features"))
			result.parameters.hdpCouplingMode =
				readBool(*value, "$.features.weighted_hdp")
					? HdpCouplingMode::VarianceWeighted
					: HdpCouplingMode::Decoupled;
		if (const auto* value = optionalMember(
				*features, "weighted_fourier_resampling", "$.features"))
			result.parameters.weightedFourierResampling =
				readBool(*value, "$.features.weighted_fourier_resampling");
		if (const auto* value =
				optionalMember(*features, "partial_resampling", "$.features"))
			result.parameters.partialResampling =
				readBool(*value, "$.features.partial_resampling");
		if (const auto* value = optionalMember(
				*features, "adaptive_effective_weights", "$.features"))
			result.parameters.effectiveWeightPolicy =
				readBool(*value, "$.features.adaptive_effective_weights")
					? EffectiveWeightPolicy::QuadraticAdaptive
					: EffectiveWeightPolicy::Fixed;
		if (const auto* value = optionalMember(
				*features, "linearized_coulomb", "$.features"))
			result.parameters.collisionCoupling =
				readBool(*value, "$.features.linearized_coulomb")
					? CollisionCoupling::Linearized
					: CollisionCoupling::Standard;
		if (const auto* value =
				optionalMember(*features, "delta_m_sampling", "$.features"))
			result.parameters.deltaMMode =
				readBool(*value, "$.features.delta_m_sampling")
					? DeltaMMode::Enabled
					: DeltaMMode::Disabled;
		if (const auto* value =
				optionalMember(*features, "projection_mode", "$.features")) {
			const auto mode = readString(*value, "$.features.projection_mode");
			if (mode == "full_micro_macro")
				result.parameters.projectionMode = ProjectionMode::FullMicroMacro;
			else if (mode == "maxwellian_only")
				result.parameters.projectionMode = ProjectionMode::MaxwellianOnly;
			else
				throw std::invalid_argument(
					"Invalid configuration at $.features.projection_mode: "
					"expected full_micro_macro or maxwellian_only");
		}
		if (const auto* value = optionalMember(
				*features, "coulomb_bgk_hybrid", "$.features"))
			result.parameters.coulombBgkHybrid =
				readBool(*value, "$.features.coulomb_bgk_hybrid");
	}

	if (const auto* resampling = optionalMember(root, "resampling", "$")) {
		rejectUnknown(*resampling, "$.resampling",
					  {"partial_cutoff_standard_deviations",
					   "conserve_weighted_moments", "signed_weight_min",
					   "signed_weight_max", "full_weight_min", "full_weight_max",
					   "cpu_cost_constant",
					   "cpu_cost_collision_coefficient"});
		auto readOptionalNumber = [&](const char* key, double& destination) {
			if (const auto* value = optionalMember(*resampling, key, "$.resampling"))
				destination = readNumber(*value,
										 "$.resampling." + std::string(key));
		};
		readOptionalNumber("partial_cutoff_standard_deviations",
						  result.parameters.partialResamplingCutoff);
		if (const auto* value = optionalMember(
				*resampling, "conserve_weighted_moments", "$.resampling"))
			result.parameters.conserveWeightedMoments =
				readBool(*value, "$.resampling.conserve_weighted_moments");
		readOptionalNumber("signed_weight_min", result.parameters.signedWeightMin);
		readOptionalNumber("signed_weight_max", result.parameters.signedWeightMax);
		readOptionalNumber("full_weight_min", result.parameters.fullWeightMin);
		readOptionalNumber("full_weight_max", result.parameters.fullWeightMax);
		readOptionalNumber("cpu_cost_constant",
						  result.parameters.cpuCostConstant);
		readOptionalNumber("cpu_cost_collision_coefficient",
						  result.parameters.cpuCostCollisionCoefficient);
	}

	if (const auto* collisions = optionalMember(root, "collisions", "$")) {
		rejectUnknown(*collisions, "$.collisions", {"bgk_strength"});
		if (const auto* value =
				optionalMember(*collisions, "bgk_strength", "$.collisions"))
			result.parameters.bgkStrength =
				readNumber(*value, "$.collisions.bgk_strength");
	}

	if (const auto* runtime = optionalMember(root, "runtime", "$")) {
		rejectUnknown(*runtime, "$.runtime",
					  {"seed", "steps", "threads", "output_directory"});
		if (const auto* value = optionalMember(*runtime, "seed", "$.runtime")) {
			if (value->type != JsonType::Null)
				result.runtime.seed = readUint32(*value, "$.runtime.seed");
		}
		if (const auto* value = optionalMember(*runtime, "steps", "$.runtime")) {
			if (value->type != JsonType::Null)
				result.runtime.steps = readPositiveInt(*value, "$.runtime.steps");
		}
		if (const auto* value =
				optionalMember(*runtime, "threads", "$.runtime")) {
			if (value->type != JsonType::Null)
				result.runtime.threads =
					readPositiveInt(*value, "$.runtime.threads");
		}
		if (const auto* value =
				optionalMember(*runtime, "output_directory", "$.runtime")) {
			if (value->type != JsonType::Null) {
				const auto directory =
					readString(*value, "$.runtime.output_directory");
				if (directory.empty())
					throw std::invalid_argument(
						"Invalid configuration at $.runtime.output_directory: "
						"must not be empty");
				result.runtime.outputDirectory = directory;
			}
		}
	}

	validateParameters(result.parameters);
	return result;
}

std::string ConfigFile::serialize(const ConfigFileValues& values) {
	validateParameters(values.parameters);
	std::ostringstream output;
	output << std::setprecision(17) << "{\n"
		   << "  \"schema_version\": 1,\n"
		   << "  \"features\": {\n"
		   << "    \"weighted_hdp\": ";
	appendBool(output,
			   values.parameters.hdpCouplingMode ==
				   HdpCouplingMode::VarianceWeighted);
	output << ",\n    \"weighted_fourier_resampling\": ";
	appendBool(output, values.parameters.weightedFourierResampling);
	output << ",\n    \"partial_resampling\": ";
	appendBool(output, values.parameters.partialResampling);
	output << ",\n    \"adaptive_effective_weights\": ";
	appendBool(output,
			   values.parameters.effectiveWeightPolicy ==
				   EffectiveWeightPolicy::QuadraticAdaptive);
	output << ",\n    \"linearized_coulomb\": ";
	appendBool(output,
			   values.parameters.collisionCoupling ==
				   CollisionCoupling::Linearized);
	output << ",\n    \"delta_m_sampling\": ";
	appendBool(output, values.parameters.deltaMMode == DeltaMMode::Enabled);
	output << ",\n    \"projection_mode\": ";
	appendEscaped(
		output, values.parameters.projectionMode == ProjectionMode::FullMicroMacro
					? "full_micro_macro"
					: "maxwellian_only");
	output << ",\n    \"coulomb_bgk_hybrid\": ";
	appendBool(output, values.parameters.coulombBgkHybrid);
	output << "\n  },\n  \"resampling\": {\n"
		   << "    \"partial_cutoff_standard_deviations\": "
		   << values.parameters.partialResamplingCutoff
		   << ",\n    \"conserve_weighted_moments\": ";
	appendBool(output, values.parameters.conserveWeightedMoments);
	output << ",\n    \"signed_weight_min\": "
		   << values.parameters.signedWeightMin
		   << ",\n    \"signed_weight_max\": "
		   << values.parameters.signedWeightMax
		   << ",\n    \"full_weight_min\": "
		   << values.parameters.fullWeightMin
		   << ",\n    \"full_weight_max\": "
		   << values.parameters.fullWeightMax
		   << ",\n    \"cpu_cost_constant\": "
		   << values.parameters.cpuCostConstant
		   << ",\n    \"cpu_cost_collision_coefficient\": "
		   << values.parameters.cpuCostCollisionCoefficient
		   << "\n  },\n  \"collisions\": {\n"
		   << "    \"bgk_strength\": " << values.parameters.bgkStrength
		   << "\n  },\n  \"runtime\": {\n"
		   << "    \"seed\": ";
	appendOptionalUint(output, values.runtime.seed);
	output << ",\n    \"steps\": ";
	appendOptionalInt(output, values.runtime.steps);
	output << ",\n    \"threads\": ";
	appendOptionalInt(output, values.runtime.threads);
	output << ",\n    \"output_directory\": ";
	appendOptionalString(output, values.runtime.outputDirectory);
	output << "\n  }\n}\n";
	return output.str();
}

} // namespace coulomb
