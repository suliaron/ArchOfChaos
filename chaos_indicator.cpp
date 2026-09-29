#include "chaos_indicator.h"

#include "math_utils.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chaos_indicator {

    double computeFLI(const double *dy, std::size_t n, double previous_fli, double norm_0)
    {
        const double norm = astro::norm(dy, n);

        if (norm == 0.0) {
            throw std::domain_error("Cannot compute FLI from a zero deviation vector.");
        }
        if (norm_0 == 0.0) {
            throw std::domain_error("Cannot compute FLI from a zero initial deviation vector.");
        }

        return std::max(previous_fli, std::log(norm / norm_0));
    }

    double computeLCI(double t, double t0, const double *dy, const double *dy0, std::size_t n)
    {
        const double elapsedTime = std::abs(t - t0);

        if (elapsedTime == 0.0) {
            throw std::domain_error("Cannot compute LCI at the initial time.");
        }

        const double initialNorm = astro::norm(dy0, n);
        const double currentNorm = astro::norm(dy, n);
        if (initialNorm == 0.0) {
            throw std::domain_error("Cannot compute LCI from a zero initial deviation vector.");
        }

        if (currentNorm == 0.0) {
            throw std::domain_error("Cannot compute LCI from a zero deviation vector.");
        }

        return std::log(currentNorm / initialNorm) / elapsedTime;
    }

}  // namespace chaos_indicator
