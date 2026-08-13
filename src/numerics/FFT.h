#pragma once

#include <complex>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "fftw3.h"

namespace coulomb {

struct FftwBufferDeleter {
	void operator()(fftw_complex* buffer) const noexcept;
};

struct FftwPlanDeleter {
	void operator()(std::remove_pointer_t<fftw_plan>* plan) const noexcept;
};

using FftwBuffer = std::unique_ptr<fftw_complex[], FftwBufferDeleter>;
using FftwPlan =
	std::unique_ptr<std::remove_pointer_t<fftw_plan>, FftwPlanDeleter>;

using Vector3D = std::vector<std::vector<std::vector<double>>>;
using VectorComplex3D =
	std::vector<std::vector<std::vector<std::complex<double>>>>;

class TensorReshape {
  public:
	template <typename T>
	std::vector<std::vector<std::vector<T>>>
	to3D(const std::vector<T>& vec1D, size_t n1, size_t n2, size_t n3) {
		if (vec1D.size() != n1 * n2 * n3)
			throw std::invalid_argument(
				"reshape1dTo3d - output vector sizes do not match input size");

		auto vec3D = std::vector(n1, std::vector(n2, std::vector<T>(n3)));
		for (size_t kk1 = 0; kk1 < n1; kk1++) {
			for (size_t kk2 = 0; kk2 < n2; kk2++) {
				for (size_t kk3 = 0; kk3 < n3; kk3++) {
					const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
					vec3D[kk1][kk2][kk3] = vec1D[kk];
				}
			}
		}
		return vec3D;
	}

	template <typename T>
	std::vector<T> to1D(const std::vector<std::vector<std::vector<T>>>& vec3D) {
		if (vec3D.empty() || vec3D.front().empty() ||
			vec3D.front().front().empty())
			throw std::invalid_argument(
				"reshape3dTo1d - input vector is empty");
		const auto n1 = vec3D.size();
		const auto n2 = vec3D.front().size();
		const auto n3 = vec3D.front().front().size();

		for (const auto& plane : vec3D) {
			if (plane.size() != n2)
				throw std::invalid_argument(
					"reshape3dTo1d - input is not rectangular");
			for (const auto& row : plane) {
				if (row.size() != n3)
					throw std::invalid_argument(
						"reshape3dTo1d - input is not rectangular");
			}
		}

		auto vec1D = std::vector<T>(n1 * n2 * n3);
		for (size_t kk1 = 0; kk1 < n1; kk1++) {
			for (size_t kk2 = 0; kk2 < n2; kk2++) {
				for (size_t kk3 = 0; kk3 < n3; kk3++) {
					const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
					vec1D[kk] = vec3D[kk1][kk2][kk3];
				}
			}
		}
		return vec1D;
	}
};

class Fft1D {
  public:
	Fft1D(size_t n);
	~Fft1D();

	Fft1D(const Fft1D&) = delete;
	Fft1D& operator=(const Fft1D&) = delete;
	Fft1D(Fft1D&&) = delete;
	Fft1D& operator=(Fft1D&&) = delete;

	std::vector<std::complex<double>> fft(const std::vector<double>& func);
	std::vector<double> ifft(const std::vector<std::complex<double>>& funcFft);

  private:
	size_t n;
	FftwBuffer inputBuffer, transformBuffer;
	FftwPlan fftPlan, ifftPlan;
};

class Fft3D {
  public:
	Fft3D(size_t n1, size_t n2, size_t n3);
	~Fft3D();

	Fft3D(const Fft3D&) = delete;
	Fft3D& operator=(const Fft3D&) = delete;
	Fft3D(Fft3D&&) = delete;
	Fft3D& operator=(Fft3D&&) = delete;

	VectorComplex3D fft(const Vector3D& func);
	Vector3D ifft(const VectorComplex3D& funcFft);

  private:
	size_t n1, n2, n3;
	FftwBuffer inputBuffer, transformBuffer;
	FftwPlan fftPlan, ifftPlan;
};

} // namespace coulomb
