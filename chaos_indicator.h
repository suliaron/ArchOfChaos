#pragma once

#include <cstddef>

namespace chaos_indicator {
    /**
     * @brief Computes the Fast Lyapunov Indicator from a deviation vector.
     *
     * The FLI is defined as the maximum logarithmic growth of the deviation
     * vector relative to its initial norm:
     *
     * \f[
     *     \mathrm{FLI}(t)
     *     =
     *     \max_{\tau \leq t}
     *     \log
     *     \left(
     *         \frac{\|\delta(\tau)\|}
     *              {\|\delta(t_0)\|}
     *     \right).
     * \f]
     *
     * @param dy Current deviation vector.
     * @param n Number of components in the deviation vector.
     * @param previous_fli Largest FLI value obtained before the current step.
     * @param norm_0 Initial norm of the deviation vector.
     *
     * @return Updated FLI value.
     *
     * @throws std::domain_error If the current or initial deviation-vector
     *         norm is zero.
     */
    double computeFLI(const double *dy, std::size_t n, double previous_fli, double norm_0);

    /**
     * @brief Computes the Lyapunov Characteristic Indicator (LCI).
     *
     * The LCI is calculated as
     *
     * \f[
     *     \mathrm{LCI}(t)
     *     =
     *     \frac{1}{|t-t_0|}
     *     \ln
     *     \left(
     *         \frac{\|\delta(t)\|}
     *              {\|\delta(t_0)\|}
     *     \right).
     * \f]
     *
     * @param t Current physical time.
     * @param t0 Initial physical time.
     * @param dy Current deviation vector.
     * @param dy0 Initial deviation vector.
     * @param n Number of components in the deviation vectors.
     *
     * @return Current LCI value.
     *
     * @throws std::domain_error If the elapsed time is zero or if either
     *         deviation-vector norm is zero.
     */
    double computeLCI(double t, double t0, const double *dy, const double *dy0, std::size_t n);

}  // namespace chaos_indicator
