#include "FFT.h"

#include <limits>
#include <stdexcept>
#include <string>

namespace coulomb {

namespace {

int fftwDimension(size_t value) {
	if (value > static_cast<size_t>(std::numeric_limits<int>::max()))
		throw std::invalid_argument("FFTW dimensions must fit in an int");
	return static_cast<int>(value);
}

size_t checkedElementCount(size_t n1, size_t n2, size_t n3 = 1) {
	if (n1 != 0 && n2 > std::numeric_limits<size_t>::max() / n1)
		throw std::length_error("FFT buffer size overflows size_t");
	const auto planeCount = n1 * n2;
	if (planeCount != 0 && n3 > std::numeric_limits<size_t>::max() / planeCount)
		throw std::length_error("FFT buffer size overflows size_t");
	const auto elementCount = planeCount * n3;
	if (elementCount >
		std::numeric_limits<size_t>::max() / sizeof(fftw_complex))
		throw std::length_error("FFT buffer allocation size overflows size_t");
	return elementCount;
}

template <typename T>
void validateShape(const std::vector<std::vector<std::vector<T>>>& values,
				   size_t n1, size_t n2, size_t n3, const char* operation) {
	if (values.size() != n1)
		throw std::invalid_argument(std::string(operation) +
									" first dimension does not match plan");
	for (const auto& plane : values) {
		if (plane.size() != n2)
			throw std::invalid_argument(
				std::string(operation) +
				" second dimension does not match plan");
		for (const auto& row : plane) {
			if (row.size() != n3)
				throw std::invalid_argument(
					std::string(operation) +
					" third dimension does not match plan");
		}
	}
}

} // namespace

void FftwBufferDeleter::operator()(fftw_complex* buffer) const noexcept {
	if (buffer != nullptr)
		fftw_free(buffer);
}

void FftwPlanDeleter::operator()(
	std::remove_pointer_t<fftw_plan>* plan) const noexcept {
	if (plan != nullptr)
		fftw_destroy_plan(plan);
}

Fft1D::Fft1D(size_t n)
	: n(n), inputBuffer(nullptr), transformBuffer(nullptr), fftPlan(nullptr),
	  ifftPlan(nullptr) {
	if (n == 0)
		throw std::invalid_argument("Fft1D dimension must be positive");
	const int nInt = fftwDimension(n);
	const auto elementCount = checkedElementCount(n, 1);
	inputBuffer.reset(static_cast<fftw_complex*>(
		fftw_malloc(elementCount * sizeof(fftw_complex))));
	transformBuffer.reset(static_cast<fftw_complex*>(
		fftw_malloc(elementCount * sizeof(fftw_complex))));
	if (!inputBuffer || !transformBuffer)
		throw std::runtime_error("Unable to allocate FFTW 1D buffers");

	fftPlan.reset(fftw_plan_dft_1d(nInt, inputBuffer.get(),
								   transformBuffer.get(), FFTW_FORWARD,
								   FFTW_ESTIMATE));
	ifftPlan.reset(fftw_plan_dft_1d(nInt, transformBuffer.get(),
									inputBuffer.get(), FFTW_BACKWARD,
									FFTW_ESTIMATE));
	if (!fftPlan || !ifftPlan)
		throw std::runtime_error("Unable to create FFTW 1D plans");
}

Fft1D::~Fft1D() = default;

std::vector<std::complex<double>> Fft1D::fft(const std::vector<double>& func) {
	if (func.size() != n)
		throw std::invalid_argument(
			"Fft1D::fft input size does not match plan");

	for (size_t kx = 0; kx < n; kx++) {
		inputBuffer[kx][0] = func[kx];
		inputBuffer[kx][1] = 0;
	}

	fftw_execute(fftPlan.get());

	std::vector<std::complex<double>> funcFft(n);
	for (size_t kx = 0; kx < n; ++kx) {
		funcFft[kx] = {transformBuffer[kx][0], transformBuffer[kx][1]};
	}
	return funcFft;
}

std::vector<double>
Fft1D::ifft(const std::vector<std::complex<double>>& funcFft) {
	if (funcFft.size() != n)
		throw std::invalid_argument(
			"Fft1D::ifft input size does not match plan");

	for (size_t kx = 0; kx < n; kx++) {
		transformBuffer[kx][0] = funcFft[kx].real();
		transformBuffer[kx][1] = funcFft[kx].imag();
	}

	fftw_execute(ifftPlan.get());

	std::vector<double> result(n);
	for (size_t kx = 0; kx < n; ++kx) {
		result[kx] = inputBuffer[kx][0];
	}
	return result;
}

Fft3D::Fft3D(size_t n1, size_t n2, size_t n3)
	: n1(n1), n2(n2), n3(n3), inputBuffer(nullptr), transformBuffer(nullptr),
	  fftPlan(nullptr), ifftPlan(nullptr) {
	if (n1 == 0 || n2 == 0 || n3 == 0)
		throw std::invalid_argument("Fft3D dimensions must be positive");
	const int n1Int = fftwDimension(n1);
	const int n2Int = fftwDimension(n2);
	const int n3Int = fftwDimension(n3);
	const auto elementCount = checkedElementCount(n1, n2, n3);
	inputBuffer.reset(static_cast<fftw_complex*>(
		fftw_malloc(elementCount * sizeof(fftw_complex))));
	transformBuffer.reset(static_cast<fftw_complex*>(
		fftw_malloc(elementCount * sizeof(fftw_complex))));
	if (!inputBuffer || !transformBuffer)
		throw std::runtime_error("Unable to allocate FFTW 3D buffers");

	fftPlan.reset(fftw_plan_dft_3d(n1Int, n2Int, n3Int, inputBuffer.get(),
								   transformBuffer.get(), FFTW_FORWARD,
								   FFTW_ESTIMATE));
	ifftPlan.reset(fftw_plan_dft_3d(n1Int, n2Int, n3Int, transformBuffer.get(),
									inputBuffer.get(), FFTW_BACKWARD,
									FFTW_ESTIMATE));
	if (!fftPlan || !ifftPlan)
		throw std::runtime_error("Unable to create FFTW 3D plans");
}

Fft3D::~Fft3D() = default;

VectorComplex3D Fft3D::fft(const Vector3D& func) {
	validateShape(func, n1, n2, n3, "Fft3D::fft input");
	// the (i,j,k)-th element of the array with size (nx,Ny,Nz), you would use
	// the expression an_array[k + Nz * (j + Ny * i)].

	for (size_t kk1 = 0; kk1 < n1; kk1++) {
		for (size_t kk2 = 0; kk2 < n2; kk2++) {
			for (size_t kk3 = 0; kk3 < n3; kk3++) {
				const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
				inputBuffer[kk][0] = func[kk1][kk2][kk3];
				inputBuffer[kk][1] = 0;
			}
		}
	}

	fftw_execute(fftPlan.get());
	VectorComplex3D funcFft =
		std::vector(n1, std::vector(n2, std::vector<std::complex<double>>(n3)));

	for (size_t kk1 = 0; kk1 < n1; kk1++) {
		for (size_t kk2 = 0; kk2 < n2; kk2++) {
			for (size_t kk3 = 0; kk3 < n3; kk3++) {
				const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
				funcFft[kk1][kk2][kk3] = {transformBuffer[kk][0],
										  transformBuffer[kk][1]};
			}
		}
	}

	return funcFft;
}

Vector3D Fft3D::ifft(const VectorComplex3D& funcFft) {
	validateShape(funcFft, n1, n2, n3, "Fft3D::ifft input");
	for (size_t kk1 = 0; kk1 < n1; kk1++) {
		for (size_t kk2 = 0; kk2 < n2; kk2++) {
			for (size_t kk3 = 0; kk3 < n3; kk3++) {
				const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
				transformBuffer[kk][0] = funcFft[kk1][kk2][kk3].real();
				transformBuffer[kk][1] = funcFft[kk1][kk2][kk3].imag();
			}
		}
	}

	fftw_execute(ifftPlan.get());

	Vector3D result = std::vector(n1, std::vector(n2, std::vector<double>(n3)));
	for (size_t kk1 = 0; kk1 < n1; kk1++) {
		for (size_t kk2 = 0; kk2 < n2; kk2++) {
			for (size_t kk3 = 0; kk3 < n3; kk3++) {
				const auto kk = kk3 + n3 * (kk2 + n2 * kk1);
				result[kk1][kk2][kk3] = inputBuffer[kk][0];
			}
		}
	}

	return result;
}

} // namespace coulomb
