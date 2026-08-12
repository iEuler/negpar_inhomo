#include "FFT.h"

#include <limits>
#include <stdexcept>
#include <string>

namespace coulomb {

namespace {

int fftw_dimension(size_t value) {
	if (value > static_cast<size_t>(std::numeric_limits<int>::max()))
		throw std::invalid_argument("FFTW dimensions must fit in an int");
	return static_cast<int>(value);
}

size_t checked_element_count(size_t n1, size_t n2, size_t n3 = 1) {
	if (n1 != 0 && n2 > std::numeric_limits<size_t>::max() / n1)
		throw std::length_error("FFT buffer size overflows size_t");
	const auto plane_count = n1 * n2;
	if (plane_count != 0 &&
		n3 > std::numeric_limits<size_t>::max() / plane_count)
		throw std::length_error("FFT buffer size overflows size_t");
	const auto element_count = plane_count * n3;
	if (element_count >
		std::numeric_limits<size_t>::max() / sizeof(fftw_complex))
		throw std::length_error("FFT buffer allocation size overflows size_t");
	return element_count;
}

template <typename T>
void validate_shape(const std::vector<std::vector<std::vector<T>>>& values,
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

void FFTWBufferDeleter::operator()(fftw_complex* buffer) const noexcept {
	if (buffer != nullptr)
		fftw_free(buffer);
}

void FFTWPlanDeleter::operator()(
	std::remove_pointer_t<fftw_plan>* plan) const noexcept {
	if (plan != nullptr)
		fftw_destroy_plan(plan);
}

FFT1D::FFT1D(size_t n)
	: n_(n), func_(nullptr), funcFFT_(nullptr), fftPlan_(nullptr),
	  ifftPlan_(nullptr) {
	if (n_ == 0)
		throw std::invalid_argument("FFT1D dimension must be positive");
	const int n_int = fftw_dimension(n_);
	const auto element_count = checked_element_count(n_, 1);
	func_.reset(static_cast<fftw_complex*>(
		fftw_malloc(element_count * sizeof(fftw_complex))));
	funcFFT_.reset(static_cast<fftw_complex*>(
		fftw_malloc(element_count * sizeof(fftw_complex))));
	if (!func_ || !funcFFT_)
		throw std::runtime_error("Unable to allocate FFTW 1D buffers");

	fftPlan_.reset(fftw_plan_dft_1d(n_int, func_.get(), funcFFT_.get(),
									FFTW_FORWARD, FFTW_ESTIMATE));
	ifftPlan_.reset(fftw_plan_dft_1d(n_int, funcFFT_.get(), func_.get(),
									 FFTW_BACKWARD, FFTW_ESTIMATE));
	if (!fftPlan_ || !ifftPlan_)
		throw std::runtime_error("Unable to create FFTW 1D plans");
}

FFT1D::~FFT1D() = default;

std::vector<std::complex<double>> FFT1D::fft(const std::vector<double>& func) {
	if (func.size() != n_)
		throw std::invalid_argument(
			"FFT1D::fft input size does not match plan");

	for (size_t kx = 0; kx < n_; kx++) {
		func_[kx][0] = func[kx];
		func_[kx][1] = 0;
	}

	fftw_execute(fftPlan_.get());

	std::vector<std::complex<double>> funcFFT(n_);
	for (size_t kx = 0; kx < n_; ++kx) {
		funcFFT[kx] = {funcFFT_[kx][0], funcFFT_[kx][1]};
	}
	return funcFFT;
}

std::vector<double>
FFT1D::ifft(const std::vector<std::complex<double>>& funcFFT) {
	if (funcFFT.size() != n_)
		throw std::invalid_argument(
			"FFT1D::ifft input size does not match plan");

	for (size_t kx = 0; kx < n_; kx++) {
		funcFFT_[kx][0] = funcFFT[kx].real();
		funcFFT_[kx][1] = funcFFT[kx].imag();
	}

	fftw_execute(ifftPlan_.get());

	std::vector<double> func(n_);
	for (size_t kx = 0; kx < n_; ++kx) {
		func[kx] = func_[kx][0];
	}
	return func;
}

FFT3D::FFT3D(size_t n1, size_t n2, size_t n3)
	: n1_(n1), n2_(n2), n3_(n3), func_(nullptr), funcFFT_(nullptr),
	  fftPlan_(nullptr), ifftPlan_(nullptr) {
	if (n1_ == 0 || n2_ == 0 || n3_ == 0)
		throw std::invalid_argument("FFT3D dimensions must be positive");
	const int n1_int = fftw_dimension(n1_);
	const int n2_int = fftw_dimension(n2_);
	const int n3_int = fftw_dimension(n3_);
	const auto element_count = checked_element_count(n1_, n2_, n3_);
	func_.reset(static_cast<fftw_complex*>(
		fftw_malloc(element_count * sizeof(fftw_complex))));
	funcFFT_.reset(static_cast<fftw_complex*>(
		fftw_malloc(element_count * sizeof(fftw_complex))));
	if (!func_ || !funcFFT_)
		throw std::runtime_error("Unable to allocate FFTW 3D buffers");

	fftPlan_.reset(fftw_plan_dft_3d(n1_int, n2_int, n3_int, func_.get(),
									funcFFT_.get(), FFTW_FORWARD,
									FFTW_ESTIMATE));
	ifftPlan_.reset(fftw_plan_dft_3d(n1_int, n2_int, n3_int, funcFFT_.get(),
									 func_.get(), FFTW_BACKWARD,
									 FFTW_ESTIMATE));
	if (!fftPlan_ || !ifftPlan_)
		throw std::runtime_error("Unable to create FFTW 3D plans");
}

FFT3D::~FFT3D() = default;

VectorComplex3D FFT3D::fft(const Vector3D& func) {
	validate_shape(func, n1_, n2_, n3_, "FFT3D::fft input");
	// the (i,j,k)-th element of the array with size (Nx,Ny,Nz), you would use
	// the expression an_array[k + Nz * (j + Ny * i)].

	for (size_t kk1 = 0; kk1 < n1_; kk1++) {
		for (size_t kk2 = 0; kk2 < n2_; kk2++) {
			for (size_t kk3 = 0; kk3 < n3_; kk3++) {
				const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
				func_[kk][0] = func[kk1][kk2][kk3];
				func_[kk][1] = 0;
			}
		}
	}

	fftw_execute(fftPlan_.get());
	VectorComplex3D funcFFT = std::vector(
		n1_, std::vector(n2_, std::vector<std::complex<double>>(n3_)));

	for (size_t kk1 = 0; kk1 < n1_; kk1++) {
		for (size_t kk2 = 0; kk2 < n2_; kk2++) {
			for (size_t kk3 = 0; kk3 < n3_; kk3++) {
				const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
				funcFFT[kk1][kk2][kk3] = {funcFFT_[kk][0], funcFFT_[kk][1]};
			}
		}
	}

	return funcFFT;
}

Vector3D FFT3D::ifft(const VectorComplex3D& funcFFT) {
	validate_shape(funcFFT, n1_, n2_, n3_, "FFT3D::ifft input");
	for (size_t kk1 = 0; kk1 < n1_; kk1++) {
		for (size_t kk2 = 0; kk2 < n2_; kk2++) {
			for (size_t kk3 = 0; kk3 < n3_; kk3++) {
				const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
				funcFFT_[kk][0] = funcFFT[kk1][kk2][kk3].real();
				funcFFT_[kk][1] = funcFFT[kk1][kk2][kk3].imag();
			}
		}
	}

	fftw_execute(ifftPlan_.get());

	Vector3D func =
		std::vector(n1_, std::vector(n2_, std::vector<double>(n3_)));
	for (size_t kk1 = 0; kk1 < n1_; kk1++) {
		for (size_t kk2 = 0; kk2 < n2_; kk2++) {
			for (size_t kk3 = 0; kk3 < n3_; kk3++) {
				const auto kk = kk3 + n3_ * (kk2 + n2_ * kk1);
				func[kk1][kk2][kk3] = func_[kk][0];
			}
		}
	}

	return func;
}

} // namespace coulomb
