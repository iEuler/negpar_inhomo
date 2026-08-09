#pragma once

#include <complex>
#include <vector>

namespace coulomb {

using Vector3D = std::vector<std::vector<std::vector<double>>>;
using VectorComplex3D =
    std::vector<std::vector<std::vector<std::complex<double>>>>;
using VectorBool3D = std::vector<std::vector<std::vector<bool>>>;

}  // namespace coulomb
