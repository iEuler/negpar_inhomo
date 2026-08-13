#include "Histogram.h"

namespace coulomb {

void Histogram::fixedBins(const std::vector<double>& values,
						  std::vector<int>& counts, double minimum,
						  double maximum) const {
	const int valueCount = static_cast<int>(values.size());
	const int binCount = static_cast<int>(counts.size());
	const double width = (maximum - minimum) / binCount;

	for (int bin = 0; bin < binCount; ++bin)
		counts[bin] = 0;

	for (int index = 0; index < valueCount; ++index) {
		int bin = static_cast<int>((values[index] - minimum) / width);
		if (bin < 0)
			bin = 0;
		if (bin >= binCount)
			bin = binCount - 1;
		++counts[bin];
	}
}

} // namespace coulomb
