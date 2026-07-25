#include <algorithm>  // std::copy, std::copy_n
#include <cstddef>    // std::size_t
#include <cstdint>    // std::uint32_t
#include <exception>  // std::exception
#include <format>     // std::format
#include <iomanip>    // std::fixed, std::setprecision, std::setw
#include <iostream>   // std::cout, std::cerr
#include <memory>     // std::unique_ptr, std::make_unique
#include <new>        // std::bad_alloc
#include <stdexcept>  // std::invalid_argument, std::runtime_error

static constexpr const char *PROGRAM_NAME    = "archofchaos";
static constexpr const char *PROGRAM_VERSION = "1.1.0";

static constexpr double PI     = 3.1415926535897932;
static constexpr double TWO_PI = 2.0 * PI;

static constexpr double kG  = 1.720209895e-2;
static constexpr double kG2 = 2.9591220828559110e-4;

typedef double var_t;
/* prototype of the right-hand side of the diff. eq. system */
typedef void rhs_t(double x, const double *y, double *dy, void *par);

/**
 * @brief Returns the square of @p x.
 *
 * Computes @c x*x using a single evaluation of the argument, avoiding the
 * multiple-evaluation pitfalls of function-like macros (e.g. @c SQR(x)).
 *
 * @tparam T  A type supporting multiplication (@c operator*).
 * @param x   Input value.
 * @return    The value @c x*x.
 *
 * @note This function does not check for overflow/underflow.
 */
template <typename T>
inline T sqr(T x) noexcept
{
    return (x * x);
}

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
    explicit GridIterator(const InitData &data)
        : data_(data),
          ia_(0),
          ie_(0),
          da_((data.a1 - data.a0) / static_cast<double>(data.Na)),
          de_((data.e1 - data.e0) / static_cast<double>(data.Ne))
    {
    }

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
    const InitData &data_;

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

namespace model {
    struct CRTBP2DParams {
        double mu;
    };

    void CRTBP2Dfun(double t, const double *y, double *dydx, void *par)
    {
        static const CRTBP2DParams *p  = static_cast<const CRTBP2DParams *>(par);
        static double               mu = (p->mu);

        double r1   = std::sqrt(sqr(y[0] + mu) + sqr(y[1]));
        double r2   = std::sqrt(sqr(y[0] - 1.0 + mu) + sqr(y[1]));
        double r1_3 = 1.0 / (r1 * r1 * r1);
        double r2_3 = 1.0 / (r2 * r2 * r2);

        // equations of motion in the rotating frame

        // dx/dt
        dydx[0] = y[2];

        // dy/dt
        dydx[1] = y[3];

        // dpx/dt
        dydx[2] = 2 * y[3] + y[0] - (1 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1 + mu) * r2_3;

        // dpy/dt
        dydx[3] = -2 * y[2] + y[1] * (1 - (1 - mu) * r1_3 - mu * r2_3);
    }
}  // namespace model

namespace ode_integrator {
    var_t rkf54(var_t t, var_t h, var_t *h_nxt, var_t *y_in, var_t *y_out, int n_var, var_t relTol, var_t absTol,
                rhs_t *fun, void *par)
    {
        static const double Bi[] = {17.0 / 192.0, 0.0, 64.0 / 231.0, 2187.0 / 8960.0, 2875.0 / 8448.0, 1.0 / 20.0, 0.0};

        static const double Ci[] = {0.0, 1.0 / 8.0, 1.0 / 4.0, 4.0 / 9.0, 4.0 / 5.0, 1.0, 1.0};

        static const double Aij[][6] = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 8.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 1.0 / 4.0, 0.0, 0.0, 0.0, 0.0},
            {196.0 / 729.0, -320.0 / 729.0, 448.0 / 729.0, 0.0, 0.0, 0.0},
            {836.0 / 2875.0, 64.0 / 575.0, -13376.0 / 20125.0, 21384.0 / 20125.0, 0.0, 0.0},
            {-73.0 / 48.0, 0.0, 1312.0 / 231.0, -2025.0 / 448.0, 2875.0 / 2112.0, 0.0},
            {17.0 / 192.0, 0.0, 64.0 / 231.0, 2187.0 / 8960.0, 2875.0 / 8448.0, 1.0 / 20.0}};

        static auto dy = allocate_array<var_t>(7 * n_var);  ///< Workspace storing the 7 Runge–Kutta stage derivatives.
        static auto y  = allocate_array<double>(n_var);     ///< Allocate memory for y

        double t0    = t;
        double h_did = 0.0;
        var_t  temax = 0.0;

        fun(t0, y_in, dy.get(), par);
        do {
            temax = 0.0;
            for (int k = 1; k < 7; k++) {
                t = t0 + Ci[k] * h;

                for (int n = 0; n < n_var; n++) {
                    y[n] = y_in[n];

                    for (int l = 0; l < k; l++)
                        y[n] += h * Aij[k][l] * dy[l * n_var + n];
                }
                fun(t, y.get(), dy.get() + k * n_var, par);
            }

            for (int n = 0; n < n_var; ++n) {
                y_out[n] = y_in[n];

                for (int k = 0; k < 7; ++k) {
                    y_out[n] += h * Bi[k] * dy[k * n_var + n];
                }
                const var_t err = h * std::fabs(dy[5 * n_var + n] - dy[6 * n_var + n]) / 60.0;
                const var_t tol = absTol + relTol * std::max(std::fabs(y_in[n]), std::fabs(y_out[n]));

                if (err / tol > temax) {
                    temax = err / tol;
                }
            }
            h_did = h;
            h     = 0.9 * h * std::pow(1.0 / temax, 1.0 / 5.0);
        } while (temax > 1.0);
        *h_nxt = h;
        return h_did;
    }
}  // namespace ode_integrator

namespace {
    void print_version()
    {
        std::cout << PROGRAM_NAME << " version " << PROGRAM_VERSION << '\n';
    }

    void getInitialCondition(double mu, double a, double e, double *x)
    {
        x[0] = a * (1.0 - e) - mu;                                      /// x_0
        x[1] = 0.0;                                                     /// y_0
        x[2] = 0.0;                                                     /// vx_0
        x[3] = std::sqrt((1 - mu) / a * (1.0 + e) / (1.0 - e)) - x[0];  /// vy_0
    }

    void printInitialCondition(const GridIterator &grid, const double *x)
    {
        static bool first = true;

        if (first) {
            std::cout << grid.header() << "         x"
                      << "         y"
                      << "         vx"
                      << "        vy" << '\n';
            first = false;
        }
        std::cout << std::fixed << std::setprecision(6) << std::setw(10) << grid.a() << std::setw(10) << grid.e()
                  << std::setw(10) << x[0] << std::setw(10) << x[1] << std::setw(10) << x[2] << std::setw(10) << x[3]
                  << '\n';
    }

}  // namespace

int main()
{
    try {
        print_version();

        model::CRTBP2DParams param;
        param.mu = 1.0e-3;

        InitData init;
        init.t0 = 0.0;
        init.a0 = 1.0;
        init.a1 = 2.0;
        init.Na = 4;

        init.e0 = 0.0;
        init.e1 = 0.2;
        init.Ne = 2;

        GridIterator grid(init);

        std::unique_ptr<double[]> y_in =
            allocate_array<double>(4);  // Allocate an array of 4 doubles for initial conditions
        std::unique_ptr<double[]> y_out =
            allocate_array<double>(4);  // Allocate an array of 4 doubles for intermedient results

        int    n_var = 4;  // Number of variables in the system (x, y, vx, vy)
        double h     = 0.1;
        // double tol = 1.0e-16;
        double relTol = 1.0e-6;
        double absTol = 1.0e-10;
        do {
            double a = grid.a();
            double e = grid.e();
            getInitialCondition(param.mu, a, e, y_in.get());  // Get initial conditions for the current (a,e)
            printInitialCondition(grid, y_in.get());          // Print the initial conditions

            double t     = init.t0;
            double h_nxt = h;
            double h_did = 0;
            do {
                h_did = ode_integrator::rkf54(t, h, &h_nxt, y_in.get(), y_out.get(), n_var, relTol, absTol,
                                              model::CRTBP2Dfun, (void *)&param);
                std::copy_n(y_out.get(), n_var, y_in.get());
                std::cout << std::fixed << std::setprecision(6) << std::setw(10) << t << std::setw(10) << y_out[0]
                          << std::setw(10) << y_out[1] << std::setw(10) << y_out[2] << std::setw(10) << y_out[3] << '\n';
                h = h_nxt;
                t += h_did;
            } while (t < 10.0);
        } while (grid.next());

        return EXIT_SUCCESS;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }

    return EXIT_FAILURE;
}