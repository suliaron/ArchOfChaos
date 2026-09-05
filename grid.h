#pragma once

#include "astro_types.h"  /**< astro::OrbitalElements. */

#include <cstddef>  /**< std::size_t   */
#include <cstdint>  /**< std::uint32_t */
#include <iosfwd>   /**< std::ostream  */
#include <iostream> /**< std::cout. */
#include <ostream>  /**< std::ostream. */
#include <string>   /**< std::string   */
#include <vector>   /**< std::vector container for grid axes and indices. */

/**
 * @brief Identifies an orbital element that can be used as a grid axis.
 *
 * The orbital phase can be represented either by the time of pericenter
 * passage (tau) or by the mean anomaly (M). These two representations
 * are mutually exclusive for a given grid definition.
 */
enum class OrbitalElement {
    SEMIMAJOR_AXIS,              /**< Semimajor axis a [AU]. */
    ECCENTRICITY,                /**< Eccentricity e. */
    INCLINATION,                 /**< Inclination i [deg]. */
    ARGUMENT_OF_PERICENTER,      /**< Argument of pericenter omega [deg]. */
    LONGITUDE_OF_ASCENDING_NODE, /**< Longitude of ascending node Omega [deg]. */
    PERICENTER_TIME,             /**< Time of pericenter passage tau [day]. */
    MEAN_ANOMALY                 /**< Mean anomaly M [deg]. */
};

/**
 * @brief Defines one axis of an orbital-element grid.
 *
 * A grid axis specifies the orbital element to be varied, its minimum
 * and maximum values, and the number of intervals between the endpoints.
 *
 * Both endpoints are included. Therefore, an axis containing
 * @p nIntervals intervals contains @p nIntervals + 1 grid points.
 */
struct GridAxis {
    /// Orbital element represented by this grid axis.
    OrbitalElement element;
    /// Minimum value of the orbital element.
    double min;
    /// Maximum value of the orbital element.
    double max;
    /**
     * @brief Number of intervals along the grid axis.
     *
     * The corresponding number of grid points is nIntervals + 1.
     */
    std::uint32_t nIntervals;
};

/**
 * @brief Converts an orbital-element name to an OrbitalElement value.
 *
 * The orbital-element name must exactly match one of the supported
 * input-file identifiers:
 *
 *     a, e, i, omega, Omega, tau, M
 *
 * The comparison is case-sensitive because omega and Omega represent
 * different orbital elements.
 *
 * @param name Orbital-element identifier used in the input file.
 *
 * @return Corresponding orbital-element type.
 *
 * @throws std::runtime_error If @p name is not a supported orbital element.
 */
OrbitalElement orbitalElementFromString(const std::string &name);

/**
 * @brief Returns the input-file identifier of an orbital element.
 *
 * @param element Orbital-element type.
 *
 * @return Orbital-element identifier.
 */
const char *orbitalElementToString(OrbitalElement element) noexcept;

/**
 * @brief Returns the input unit of an orbital element.
 *
 * @param element Orbital-element type.
 *
 * @return Unit of the orbital element.
 */
const char *orbitalElementUnit(OrbitalElement element) noexcept;

/**
 * @brief Checks whether an orbital element is periodic.
 *
 * The angular orbital elements omega, Omega, and M are periodic
 * with a period of 360 degrees.
 *
 * @param element Orbital-element type.
 *
 * @return True if the orbital element is periodic, false otherwise.
 */
bool isPeriodic(OrbitalElement element) noexcept;

/**
 * @brief Checks whether a grid axis covers one complete angular period.
 *
 * A full-period axis is a periodic orbital element whose range extends
 * from 0 to 360 degrees. For such an axis, the upper endpoint is not
 * included because 360 degrees is equivalent to 0 degrees.
 *
 * @param axis Grid-axis definition.
 *
 * @return True if the axis covers the full interval [0, 360] degrees,
 *         false otherwise.
 */
bool isFullPeriod(const GridAxis &axis) noexcept;

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
     * @brief Constructs an iterator for an arbitrary orbital-element grid.
     *
     * The grid axes are stored in the order in which they are specified.
     * The first axis varies fastest during iteration.
     *
     * All grid indices are initialized to zero, so the iterator initially
     * represents the first grid point.
     *
     * @param axes Definitions of the orbital-element grid axes.
     */
    explicit GridIterator(const std::vector<GridAxis> &axes);

    /**
     * @brief Returns the number of grid axes.
     *
     * @return Number of orbital-element grid axes.
     */
    std::size_t getAxisCount() const noexcept
    {
        return axes_.size();
    }

    /**
     * @brief Returns the current value along a grid axis.
     *
     * The value is calculated from the minimum axis value, the number of
     * intervals, and the current index as
     *
     * @f[
     * q = q_{\min}
     *     + j \frac{q_{\max} - q_{\min}}{N},
     * @f]
     *
     * where @f$j@f$ is the current grid index and @f$N@f$ is the number
     * of intervals along the selected axis.
     *
     * @param axisIndex Zero-based index of the grid axis.
     *
     * @return Current value along the selected grid axis.
     *
     * @throws std::out_of_range If @p axisIndex is outside the range of
     *         available grid axes.
     */
    double getValue(std::size_t axisIndex) const;

    /**
     * @brief Returns the current grid value of an orbital element.
     *
     * Searches the grid axes for the specified orbital element and returns
     * its value at the current grid point.
     *
     * @param element Orbital element whose current grid value is requested.
     *
     * @return Current value of the specified orbital element.
     *
     * @throws std::out_of_range If the orbital element is not present
     *         among the grid axes.
     */
    double getValue(OrbitalElement element) const;

    /**
     * @brief Returns the number of grid points along an axis.
     *
     * Normally, both endpoints of a grid interval are included, so an axis
     * with N intervals contains N + 1 grid points.
     *
     * For a periodic axis covering the complete interval from 0 to 360
     * degrees, the upper endpoint is excluded because it is equivalent to
     * the lower endpoint. Such an axis therefore contains N grid points.
     *
     * @param axisIndex Zero-based index of the grid axis.
     *
     * @return Number of distinct grid points along the selected axis.
     *
     * @throws std::out_of_range If @p axisIndex is outside the range of
     *         available grid axes.
     */
    std::size_t getPointCount(std::size_t axisIndex) const;

    /**
     * @brief Returns the total number of points in the grid.
     *
     * The total number of grid points is the product of the number of
     * points along all grid axes. Periodic full-period axes are handled
     * according to getPointCount().
     *
     * @return Total number of grid points.
     */
    std::size_t getTotalPointCount() const;

    /**
     * @brief Checks whether an orbital element is present as a grid axis.
     *
     * @param element Orbital element to search for.
     *
     * @return True if the orbital element is present among the grid axes,
     *         false otherwise.
     */
    bool hasAxis(OrbitalElement element) const noexcept;

    /**
     * @brief Applies the current grid values to a set of orbital elements.
     *
     * Grid-controlled orbital elements are replaced by their current grid
     * values, while orbital elements that are not grid axes remain unchanged.
     *
     * Angular grid values i, omega, and Omega are specified in degrees and
     * are converted to radians before being stored in @p elements.
     *
     * The time of pericenter passage tau is copied directly in days.
     *
     * Mean anomaly M is not stored in astro::OrbitalElements and is therefore
     * not applied by this function. If M is a grid axis, its current value
     * must be handled separately by the caller.
     *
     * @param elements Orbital elements to modify.
     */
    void apply(astro::OrbitalElements& elements) const;

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
     * For a multidimensional grid, the first axis is the fastest-varying
     * axis. When an axis reaches its last grid point, it is reset to zero
     * and the next axis is incremented.
     *
     * The number of points along each axis is determined by getPointCount(),
     * so full-period angular axes automatically exclude the duplicated
     * 360-degree endpoint.
     *
     * @return True if the iterator was advanced to another grid point,
     *         false if the current point is the last point of the grid.
     */
    bool next();

    /**
     * @brief Prints the current grid-integration progress.
     *
     * Displays the current grid-point number, the total number of grid
     * points, a graphical progress bar, and the completion percentage.
     *
     * The progress calculation is independent of the number of grid axes.
     *
     * @param os Output stream. Defaults to std::cout.
     */
    void printProgress(std::ostream &os = std::cout) const;

    /**
     * @brief Prints the current multidimensional grid point.
     *
     * Displays the zero-based linear grid-point index and, for every
     * grid axis, the current axis index, orbital-element name, and
     * corresponding grid value.
     *
     * This function is intended primarily for testing and debugging
     * the multidimensional grid traversal.
     *
     * @param os Output stream. Defaults to std::cout.
     */
    void printCurrentPoint(std::ostream &os = std::cout) const;

   private:
    /**
     * @brief Returns the grid-axis index of an orbital element.
     *
     * @param element Orbital element to search for.
     *
     * @return Zero-based index of the corresponding grid axis.
     *
     * @throws std::out_of_range If the orbital element is not present
     *         among the grid axes.
     */
    std::size_t getAxisIndex(OrbitalElement element) const;

    /**
     * @brief Definitions of the orbital-element grid axes.
     *
     * The axes are stored in the same order in which they were supplied
     * to the constructor. The first axis is the fastest-varying axis.
     */
    std::vector<GridAxis> axes_;

    /**
     * @brief Current index along each grid axis.
     *
     * The size of this vector is identical to the number of grid axes.
     * All indices are initialized to zero.
     */
    std::vector<std::uint32_t> indices_;
    /**
     * @brief Zero-based linear index of the current grid point.
     *
     * The first grid point has index zero. The value is incremented
     * whenever next() successfully advances the iterator.
     */
    std::size_t currentPoint_;

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