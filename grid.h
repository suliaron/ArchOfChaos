#pragma once

#include <cstdint>  // std::uint32_t
#include <iosfwd>   // std::ostream
#include <string>   // std::string

/**
 * @brief Iterates over a two-dimensional (a,e) parameter grid.
 *
 * The iterator traverses the grid in row-major order: the semimajor axis
 * is incremented first, followed by the eccentricity.
 *
 * Both endpoints of the semimajor-axis and eccentricity intervals are
 * included in the grid.
 */
class GridIterator {
   public:
    /**
     * @brief Constructs a two-dimensional parameter-grid iterator.
     *
     * @param a0 Minimum semimajor axis.
     * @param a1 Maximum semimajor axis.
     * @param Na Number of semimajor-axis intervals.
     * @param e0 Minimum eccentricity.
     * @param e1 Maximum eccentricity.
     * @param Ne Number of eccentricity intervals.
     *
     * @throws std::invalid_argument If @p Na or @p Ne is zero.
     */
    GridIterator(double a0, double a1, std::uint32_t Na, double e0, double e1, std::uint32_t Ne);

    /**
     * @brief Returns the current semimajor axis.
     *
     * @return Current value of a.
     */
    double a() const noexcept;

    /**
     * @brief Returns the current eccentricity.
     *
     * @return Current value of e.
     */
    double e() const noexcept;

    /**
     * @brief Returns the table header.
     *
     * @return Header string.
     */
    std::string header() const;

    /**
     * @brief Advances the iterator to the next grid point.
     *
     * @return True if the next grid point exists, false otherwise.
     */
    bool next() noexcept;

    /**
     * @brief Prints the current grid progress.
     *
     * @param os Output stream.
     */
    void printProgress(std::ostream &os) const;

   private:
    /// Minimum semimajor axis.
    double a0_;

    /// Minimum eccentricity.
    double e0_;

    /// Number of semimajor-axis intervals.
    std::uint32_t Na_;

    /// Number of eccentricity intervals.
    std::uint32_t Ne_;

    /// Current semimajor-axis index.
    std::uint32_t ia_;

    /// Current eccentricity index.
    std::uint32_t ie_;

    /// Semimajor-axis step size.
    double da_;

    /// Eccentricity step size.
    double de_;
};