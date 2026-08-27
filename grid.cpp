#include "grid.h"

#include <cstddef>    // std::size_t
#include <iomanip>    // std::fixed, std::setprecision, std::setw
#include <ostream>    // std::ostream
#include <stdexcept>  // std::invalid_argument

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

bool GridIterator::next() noexcept
{
    if (ia_ < Na_) {
        ++ia_;
        return true;
    }

    ia_ = 0;

    if (ie_ < Ne_) {
        ++ie_;
        return true;
    }

    return false;
}

void GridIterator::printProgress(std::ostream &os) const
{
    constexpr std::size_t kBarWidth = 40;

    const std::size_t total = (static_cast<std::size_t>(Na_) + 1) * (static_cast<std::size_t>(Ne_) + 1);
    const std::size_t current =
        static_cast<std::size_t>(ie_) * (static_cast<std::size_t>(Na_) + 1) + static_cast<std::size_t>(ia_) + 1;
    const double progress = static_cast<double>(current) / static_cast<double>(total);
    const std::size_t filled = static_cast<std::size_t>(progress * static_cast<double>(kBarWidth));

    os << '\r' << '(' << std::setw(4) << ia_ << ", " << std::setw(4) << ie_ << ") [";

    for (std::size_t i = 0; i < kBarWidth; ++i) {
        os << (i < filled ? '#' : '.');
    }

    os << "] " << std::fixed << std::setw(6) << std::setprecision(2) << progress * 100.0 << " %" << std::flush;
}