#include <cstddef>    // std::size_t
#include <cstdint>    // std::uint32_t
#include <exception>  // std::exception
#include <iomanip>    // std::fixed, std::setprecision, std::setw
#include <iostream>   // std::cout, std::cerr
#include <memory>     // std::unique_ptr, std::make_unique
#include <new>        // std::bad_alloc
#include <stdexcept>  // std::invalid_argument, std::runtime_error

/**
 * @brief Stores the input parameters of the simulation.
 *
 * This structure contains the initial epoch, the mass parameter,
 * the definition of the (a,e) parameter grid, and the fixed orbital
 * elements used to generate the initial conditions.
 */
struct InitData {
    /// Initial epoch [day].
    double t0 = 0.0;

    /// Mass parameter.
    double mu = 0.0;

    /// Minimum semimajor axis.
    double a0 = 0.0;

    /// Maximum semimajor axis.
    double a1 = 0.0;

    /// Number of semimajor-axis intervals.
    std::uint32_t Na = 0;

    /// Minimum eccentricity.
    double e0 = 0.0;

    /// Maximum eccentricity.
    double e1 = 0.0;

    /// Number of eccentricity intervals.
    std::uint32_t Ne = 0;

    /// Argument of periapsis [rad].
    double omega = 0.0;

    /// Time of periapsis passage [day].
    double tau = 0.0;

    /// Constructs an empty initialization data structure.
    InitData() = default;

    /// Destroys the initialization data structure.
    ~InitData() = default;
};

/**
 * @brief Iterates over a two-dimensional (a,e) parameter grid.
 *
 * The iterator traverses the grid in row-major order: the semimajor axis
 * is incremented first, followed by the eccentricity.
 */
class GridIterator {
public:
    /**
     * @brief Constructs a grid iterator.
     *
     * @param data Grid definition.
     */
    explicit GridIterator(const InitData& data)
        : data_(data),
        ia_(0),
        ie_(0),
        da_((data.a1 - data.a0) / static_cast<double>(data.Na)),
        de_((data.e1 - data.e0) / static_cast<double>(data.Ne))
    {}

    /**
     * @brief Returns the current semimajor axis.
     *
     * @return Current value of a.
     */
    double a() const
    {
        return data_.a0 + static_cast<double>(ia_) * da_;
    }

    /**
     * @brief Returns the current eccentricity.
     *
     * @return Current value of e.
     */
    double e() const
    {
        return data_.e0 + static_cast<double>(ie_) * de_;
    }

    /**
     * @brief Returns the table header.
     *
     * @return Header string.
     */
    std::string header() const
    {
        return "  a         e";
    }

    /**
     * @brief Advances the iterator to the next grid point.
     *
     * @return True if the next grid point exists, false otherwise.
     */
    bool next()
    {
        if (ia_ < data_.Na) {
            ++ia_;
            return true;
        }

        ia_ = 0;

        if (ie_ < data_.Ne) {
            ++ie_;
            return true;
        }

        return false;
    }

private:
    /// Grid definition.
    const InitData& data_;

    /// Current semimajor-axis index.
    std::uint32_t ia_;

    /// Current eccentricity index.
    std::uint32_t ie_;

    /// Semimajor-axis step size.
    double da_;

    /// Eccentricity step size.
    double de_;
};

/**
 * @brief Allocates an array of n elements of type T.
 *
 * @tparam T Element type.
 * @param n Number of elements to allocate.
 * @return std::unique_ptr<T[]> Pointer to the allocated array.
 *
 * @throws std::invalid_argument If n is zero.
 * @throws std::bad_alloc If memory allocation fails.
 */
template <typename T>
std::unique_ptr<T[]> allocate_array(std::size_t n)
{
    if (n == 0) {
        throw std::invalid_argument("The array size must be greater than zero.");
    }

    return std::make_unique<T[]>(n);
}

void getInitialCondition(double a, double e, double *x)
{
}

int main()
{
    try {
        InitData init;
        init.a0 = 1.0;
        init.a1 = 2.0;
        init.Na = 4;

        init.e0 = 0.0;
        init.e1 = 0.2;
        init.Ne = 2;

        GridIterator grid(init);

        std::cout << grid.header() << '\n';
        do {
            std::cout << std::fixed << std::setprecision(6)
                << std::setw(10) << grid.a()
                << std::setw(10) << grid.e() << '\n';
        } while (grid.next());

        return EXIT_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }
    catch (...) {
        std::cerr << "Unknown error.\n";
    }

    return EXIT_FAILURE;
}