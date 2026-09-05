#include "grid.h"
#include "math_utils.h" /**< astro::toRad. */

#include <cmath>     /**< std::abs. */
#include <cstddef>   /**< std::size_t  */
#include <iomanip>   /**< std::fixed, std::setprecision, std::setw */
#include <ostream>   /**< std::ostream */
#include <stdexcept> /**< std::invalid_argument, std::out_of_range, std::runtime_error. */

#include "grid.h"

#include <stdexcept>

namespace {

    constexpr double FULL_PERIOD_DEG = 360.0;
    constexpr double ANGLE_EPS       = 1.0e-12;

}  // namespace

/**
 * @brief Converts an orbital-element name to an OrbitalElement value.
 *
 * @param name Orbital-element identifier used in the input file.
 * @return Corresponding orbital-element type.
 *
 * @throws std::runtime_error If @p name is not a supported orbital element.
 */
OrbitalElement orbitalElementFromString(const std::string &name)
{
    if (name == "a") {
        return OrbitalElement::SEMIMAJOR_AXIS;
    }

    if (name == "e") {
        return OrbitalElement::ECCENTRICITY;
    }

    if (name == "i") {
        return OrbitalElement::INCLINATION;
    }

    if (name == "omega") {
        return OrbitalElement::ARGUMENT_OF_PERICENTER;
    }

    if (name == "Omega") {
        return OrbitalElement::LONGITUDE_OF_ASCENDING_NODE;
    }

    if (name == "tau") {
        return OrbitalElement::PERICENTER_TIME;
    }

    if (name == "M") {
        return OrbitalElement::MEAN_ANOMALY;
    }

    throw std::runtime_error("Unknown orbital element in grid definition: " + name);
}

/**
 * @brief Returns the input-file identifier of an orbital element.
 *
 * @param element Orbital-element type.
 * @return Orbital-element identifier.
 */
const char *orbitalElementToString(OrbitalElement element) noexcept
{
    switch (element) {
        case OrbitalElement::SEMIMAJOR_AXIS:
            return "a";

        case OrbitalElement::ECCENTRICITY:
            return "e";

        case OrbitalElement::INCLINATION:
            return "i";

        case OrbitalElement::ARGUMENT_OF_PERICENTER:
            return "omega";

        case OrbitalElement::LONGITUDE_OF_ASCENDING_NODE:
            return "Omega";

        case OrbitalElement::PERICENTER_TIME:
            return "tau";

        case OrbitalElement::MEAN_ANOMALY:
            return "M";
    }

    return "UNKNOWN";
}

const char *orbitalElementUnit(OrbitalElement element) noexcept
{
    switch (element) {
        case OrbitalElement::SEMIMAJOR_AXIS:
            return "AU";

        case OrbitalElement::ECCENTRICITY:
            return "-";

        case OrbitalElement::INCLINATION:
        case OrbitalElement::ARGUMENT_OF_PERICENTER:
        case OrbitalElement::LONGITUDE_OF_ASCENDING_NODE:
        case OrbitalElement::MEAN_ANOMALY:
            return "deg";

        case OrbitalElement::PERICENTER_TIME:
            return "day";
    }

    return "UNKNOWN";
}

bool isPeriodic(OrbitalElement element) noexcept
{
    switch (element) {
        case OrbitalElement::ARGUMENT_OF_PERICENTER:
        case OrbitalElement::LONGITUDE_OF_ASCENDING_NODE:
        case OrbitalElement::MEAN_ANOMALY:
            return true;

        case OrbitalElement::SEMIMAJOR_AXIS:
        case OrbitalElement::ECCENTRICITY:
        case OrbitalElement::INCLINATION:
        case OrbitalElement::PERICENTER_TIME:
            return false;
    }

    return false;
}

bool isFullPeriod(const GridAxis &axis) noexcept
{
    if (!isPeriodic(axis.element)) {
        return false;
    }

    return std::abs(axis.min) <= ANGLE_EPS && std::abs(axis.max - FULL_PERIOD_DEG) <= ANGLE_EPS;
}

GridIterator::GridIterator(double a0, double a1, std::uint32_t Na, double e0, double e1, std::uint32_t Ne)
    : a0_(a0), e0_(e0), Na_(Na), Ne_(Ne), ia_(0), ie_(0), da_(0.0), de_(0.0)
{
    if (Na_ == 0) {
        throw std::invalid_argument("Number of semimajor-axis intervals must be greater than zero.");
    }

    if (Ne_ == 0) {
        throw std::invalid_argument("Number of eccentricity intervals must be greater than zero.");
    }

    da_ = (a1 - a0_) / static_cast<double>(Na_);
    de_ = (e1 - e0_) / static_cast<double>(Ne_);

    currentPoint_ = 0;
}

GridIterator::GridIterator(const std::vector<GridAxis> &axes) : axes_(axes), indices_(axes.size(), 0), currentPoint_(0)
{
}

double GridIterator::getValue(std::size_t axisIndex) const
{
    if (axisIndex >= axes_.size()) {
        throw std::out_of_range("Grid axis index is out of range.");
    }

    const GridAxis &axis = axes_[axisIndex];

    const double step = (axis.max - axis.min) / static_cast<double>(axis.nIntervals);

    return axis.min + static_cast<double>(indices_[axisIndex]) * step;
}

double GridIterator::getValue(OrbitalElement element) const
{
    return getValue(getAxisIndex(element));
}

std::size_t GridIterator::getAxisIndex(OrbitalElement element) const
{
    for (std::size_t axisIndex = 0; axisIndex < axes_.size(); ++axisIndex) {
        if (axes_[axisIndex].element == element) {
            return axisIndex;
        }
    }

    throw std::out_of_range("Orbital element '" + std::string(orbitalElementToString(element)) +
                            "' is not a grid axis.");
}

std::size_t GridIterator::getPointCount(std::size_t axisIndex) const
{
    if (axisIndex >= axes_.size()) {
        throw std::out_of_range("Grid axis index is out of range.");
    }

    const GridAxis &axis = axes_[axisIndex];

    if (isFullPeriod(axis)) {
        return static_cast<std::size_t>(axis.nIntervals);
    }

    return static_cast<std::size_t>(axis.nIntervals) + 1;
}

std::size_t GridIterator::getTotalPointCount() const
{
    // General multidimensional grid.
    if (!axes_.empty()) {
        std::size_t total = 1;

        for (std::size_t axisIndex = 0; axisIndex < axes_.size(); ++axisIndex) {
            total *= getPointCount(axisIndex);
        }

        return total;
    }

    // Legacy two-dimensional (a,e) grid.
    return (static_cast<std::size_t>(Na_) + 1) * (static_cast<std::size_t>(Ne_) + 1);
}

bool GridIterator::hasAxis(OrbitalElement element) const noexcept
{
    for (const GridAxis &axis : axes_) {
        if (axis.element == element) {
            return true;
        }
    }

    return false;
}

void GridIterator::apply(astro::OrbitalElements &elements) const
{
    for (std::size_t axisIndex = 0; axisIndex < axes_.size(); ++axisIndex) {
        const GridAxis &axis  = axes_[axisIndex];
        const double    value = getValue(axisIndex);

        switch (axis.element) {
            case OrbitalElement::SEMIMAJOR_AXIS:
                elements.a = value;
                break;

            case OrbitalElement::ECCENTRICITY:
                elements.e = value;
                break;

            case OrbitalElement::INCLINATION:
                elements.i = astro::toRad(value);
                break;

            case OrbitalElement::ARGUMENT_OF_PERICENTER:
                elements.omega = astro::toRad(value);
                break;

            case OrbitalElement::LONGITUDE_OF_ASCENDING_NODE:
                elements.Omega = astro::toRad(value);
                break;

            case OrbitalElement::PERICENTER_TIME:
                elements.tau = value;
                break;

            case OrbitalElement::MEAN_ANOMALY:
                // Mean anomaly is not stored in astro::OrbitalElements.
                // It is handled separately by the caller.
                break;
        }
    }
}

double GridIterator::a() const noexcept
{
    return a0_ + static_cast<double>(ia_) * da_;
}

double GridIterator::e() const noexcept
{
    return e0_ + static_cast<double>(ie_) * de_;
}

std::string GridIterator::header() const
{
    return "  a         e";
}

bool GridIterator::next()
{
    // ---------------------------------------------------------------------
    // General multidimensional grid.
    // ---------------------------------------------------------------------

    if (!axes_.empty()) {
        for (std::size_t axisIndex = 0; axisIndex < indices_.size(); ++axisIndex) {
            const std::size_t pointCount = getPointCount(axisIndex);

            // If the current axis can still be incremented,
            // advance it and reset all faster-varying axes.
            if (static_cast<std::size_t>(indices_[axisIndex]) + 1 < pointCount) {
                for (std::size_t i = 0; i < axisIndex; ++i) {
                    indices_[i] = 0;
                }

                ++indices_[axisIndex];
                ++currentPoint_;

                return true;
            }
        }

        // No axis can be incremented: the current point is the
        // last point of the multidimensional grid.
        return false;
    }

    // ---------------------------------------------------------------------
    // Legacy two-dimensional (a,e) grid.
    //
    // This branch is kept temporarily until runGrid() is converted to the
    // general GridAxis-based interface.
    // ---------------------------------------------------------------------

    if (ia_ < Na_) {
        ++ia_;
        ++currentPoint_;

        return true;
    }

    ia_ = 0;

    if (ie_ < Ne_) {
        ++ie_;
        ++currentPoint_;

        return true;
    }

    return false;
}

void GridIterator::printProgress(std::ostream &os) const
{
    constexpr std::size_t BAR_WIDTH = 40;

    const std::size_t total = getTotalPointCount();

    if (total == 0) {
        return;
    }

    // currentPoint_ is zero-based, while the displayed point number
    // is one-based.
    const std::size_t current  = currentPoint_ + 1;
    const double      progress = static_cast<double>(current) / static_cast<double>(total);
    const std::size_t filled   = static_cast<std::size_t>(progress * static_cast<double>(BAR_WIDTH));
    os << '\r' << '[';

    for (std::size_t i = 0; i < BAR_WIDTH; ++i) {
        os << (i < filled ? '#' : '.');
    }

    os << "] " << std::setw(6) << std::fixed << std::setprecision(2) << progress * 100.0 << " %"
       << "  (" << current << '/' << total << ')' << std::flush;
}

void GridIterator::printCurrentPoint(std::ostream &os) const
{
    os << "point " << currentPoint_ << " : ";

    for (std::size_t axisIndex = 0; axisIndex < axes_.size(); ++axisIndex) {
        const GridAxis &axis = axes_[axisIndex];

        os << orbitalElementToString(axis.element) << '[' << indices_[axisIndex] << "]=" << getValue(axisIndex);
        if (axisIndex + 1 < axes_.size()) {
            os << "  ";
        }
    }

    os << '\n';
}
