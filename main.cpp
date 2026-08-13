#include "math_utils.h"  // astro::sqr
#include "model.h"       // Model base class
#include "crtbp.h"

#include <algorithm>  /**< std::copy, std::copy_n, std::remove_if. */
#include <chrono>     /**< Provides time measurement utilities. */
#include <cmath>      /**< Mathematical functions (std::sqrt, std::abs). */
#include <cctype>     /**< Character classification functions (std::isspace). */
#include <cstddef>    /**< std::size_t. */
#include <cstdint>    /**< Fixed-width integer types (std::uint32_t). */
#include <exception>  /**< std::exception base class. */
#include <filesystem> /**< std::filesystem::path for portable path handling. */
#include <format>     /**< std::format. */
#include <fstream>    /**< File stream classes (std::ifstream, std::ofstream). */
#include <iomanip>    /**< Output manipulators (std::fixed, std::setprecision, std::setw). */
#include <iostream>   /**< Standard input/output streams (std::cout, std::cerr). */
#include <memory>     /**< Smart pointers (std::unique_ptr, std::make_unique). */
#include <new>        /**< std::bad_alloc. */
#include <sstream>    /**< String stream classes (std::istringstream, std::ostringstream). */
#include <stdexcept>  /**< Standard exception classes (std::invalid_argument, std::runtime_error). */
#include <string>     /**< std::string class. */

namespace fs = std::filesystem;

static constexpr const char *PROGRAM_NAME    = "archofchaos";
static constexpr const char *PROGRAM_VERSION = "1.1";

constexpr double TIME_EPS = 1.0e-15;

typedef double var_t;

/**
 * @brief Specifies the overall computation mode.
 */
enum class RunMode {
    ORBIT,     /**< Integrate a single orbit and save the state evolution. */
    INDICATOR, /**< Integrate a single orbit and save an indicator evolution. */
    GRID       /**< Integrate orbits over a parameter grid and save final indicator values. */
};
/**
 * @brief Returns the name of a run mode.
 *
 * @param mode Run mode.
 * @return Name of the run mode.
 */
/**
 * @brief Returns the name of a run mode.
 *
 * @param mode Run mode.
 * @return Name of the run mode.
 */
const char *RunModeToString(RunMode mode)
{
    switch (mode) {
        case RunMode::ORBIT:
            return "ORBIT";

        case RunMode::INDICATOR:
            return "INDICATOR";

        case RunMode::GRID:
            return "GRID";
    }

    return "UNKNOWN";
}

/**
 * @brief Specifies the chaos indicator to be computed.
 */
enum class IndicatorType {
    NONE, /**< No chaos indicator is computed. */
    FLI,  /**< Fast Lyapunov Indicator. */
    LCI,  /**< Lyapunov Characteristic Indicator. */
    RLI   /**< Relative Lyapunov Indicator. */
};
/**
 * @brief Returns the name of a chaos indicator.
 *
 * @param indicator Indicator type.
 * @return Name of the indicator.
 */
const char *IndicatorTypeToString(IndicatorType indicator)
{
    switch (indicator) {
        case IndicatorType::NONE:
            return "NONE";

        case IndicatorType::FLI:
            return "FLI";

        case IndicatorType::LCI:
            return "LCI";

        case IndicatorType::RLI:
            return "RLI";
    }

    return "UNKNOWN";
}

/**
 * @brief Command-line options.
 *
 * Stores the input and output file names, directories, and full paths,
 * together with flags controlling the display of help, version, and
 * verbose information.
 */
struct CommandLineOptions {
    /** Input file name without directory path. */
    std::string input_file;

    /** Full path of the input file. */
    std::string input_path;

    /** Input directory. */
    std::string input_dir;

    /** Output file name without directory path. */
    std::string output_file;

    /** Full path of the output file. */
    std::string output_path;

    /** Output directory. */
    std::string output_dir;

    /** Display help message. */
    bool show_help = false;

    /** Display program version. */
    bool show_version = false;

    /** Display verbose output. */
    bool verbose = false;

    /**
     * @brief Prints the command-line options.
     *
     * Prints the input and output file names, directories, and full paths,
     * together with the command-line flags in a human-readable format.
     *
     * @param os Output stream. Defaults to std::cout.
     */
    void Print(std::ostream &os = std::cout) const
    {
        os << "----------------------------------------\n";
        os << "Command-line options\n";
        os << "----------------------------------------\n";
        os << "Input file    : " << input_file << '\n'
           << "Input dir     : " << input_dir << '\n'
           << "Input path    : " << input_path << '\n'
           << "Output file   : " << output_file << '\n'
           << "Output dir    : " << output_dir << '\n'
           << "Output path   : " << output_path << '\n'
           << "Show help     : " << std::boolalpha << show_help << '\n'
           << "Show version  : " << std::boolalpha << show_version << '\n'
           << "Verbose       : " << std::boolalpha << verbose << '\n';
    }
};
/**
 * @brief Stores the input parameters of the simulation.
 *
 * This structure contains the initial epoch, the mass parameter,
 * the definition of the (a,e) parameter grid, and the fixed orbital
 * elements used to generate the initial conditions.
 */
struct InitData {
    /**
     * @brief Overall computation mode.
     */
    RunMode run_mode = RunMode::ORBIT;

    /**
     * @brief Chaos indicator to be computed.
     *
     * Must be IndicatorType::NONE in ORBIT mode.
     */
    IndicatorType indicator = IndicatorType::NONE;
    /// Mass parameter.
    double mu = 0.0;

    /// Initial epoch [time unit].
    double t0 = 0.0;
    /// Length of integration [time unit].
    double T = 0.0;
    /**
     * @brief Output time interval.
     *
     * Specifies the time interval between consecutive output records
     * in ORBIT and INDICATOR modes.
     */
    double output_dt = 0.0;
    /// Minimum semimajor axis.
    double a0 = 0.0;
    /// Maximum semimajor axis.
    double a1 = 0.0;
    /// Number of intervals along the semimajor-axis dimension.
    std::uint32_t Na = 0;

    /// Minimum eccentricity.
    double e0 = 0.0;
    /// Maximum eccentricity.
    double e1 = 0.0;
    /// Number of intervals along the eccentricity dimension.
    std::uint32_t Ne = 0;

    /**
     * @brief Initial condition for the variational equations.
     */
    double dy[4] = {1.0, 0.0, 0.0, 0.0};

    /// Constructs an empty initialization data structure.
    InitData() = default;

    /// Destroys the initialization data structure.
    ~InitData() = default;

    void Print(std::ostream &os = std::cout) const
    {
        os << std::scientific << std::setprecision(6) << std::showpos;

        constexpr int W = 16;

        os << "----------------------------------------\n";
        os << "InitData contents\n";
        os << "----------------------------------------\n";

        os << "run mode  : " << RunModeToString(run_mode) << '\n';
        os << "indicator : " << IndicatorTypeToString(indicator) << '\n';
        os << "output dt : " << std::setw(W) << output_dt << " [time unit]\n";
        os << "mu = " << std::setw(W) << mu << "\n";
        os << "t0 = " << std::setw(W) << t0 << " [time unit]\n";
        os << "T  = " << std::setw(W) << T << " [time unit]\n";
        os << "a0 = " << std::setw(W) << a0 << " [distance unit]\n";
        os << "a1 = " << std::setw(W) << a1 << " [distance unit]\n";
        os << "Na = " << std::setw(W) << Na << '\n';
        os << "e0 = " << std::setw(W) << e0 << '\n';
        os << "e1 = " << std::setw(W) << e1 << '\n';
        os << "Ne = " << std::setw(W) << Ne << '\n';
        if (indicator != IndicatorType::NONE) {
            os << "dy1 = " << std::setw(W) << dy[0] << '\n';
            os << "dy2 = " << std::setw(W) << dy[1] << '\n';
            os << "dy3 = " << std::setw(W) << dy[2] << '\n';
            os << "dy4 = " << std::setw(W) << dy[3] << '\n';
        }
    }
};

/**
 * @brief Adaptive step-size control parameters.
 *
 * Stores the current, accepted, proposed, and allowed integration
 * step sizes used by the adaptive Runge–Kutta integrator.
 */
struct StepControl {
    double h     = 0.0;  ///< Current integration step size.
    double h_nxt = 0.0;  ///< Proposed step size for the next step.
    double h_did = 0.0;  ///< Accepted step size of the current step.
    double h_max = 0.0;  ///< Maximum allowed step size.
    double h_min = 0.0;  ///< Minimum allowed step size.

    uint32_t n_int = 0;  ///< Number of integration steps taken.
    uint32_t n_tst = 0;  ///< Test
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
    std::string Header() const
    {
        return "  a         e";
    }

    /**
     * @brief Advances the iterator to the next grid point.
     *
     * @return True if the next grid point exists, false otherwise.
     */
    bool Next()
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

    /**
     * @brief Prints the current progress.
     *
     * @param os Output stream.
     */
    void PrintProgress(std::ostream &os = std::cout) const
    {
        constexpr std::size_t kBarWidth = 40;

        const std::size_t total    = (data_.Na + 1) * (data_.Ne + 1);
        const std::size_t current  = ie_ * (data_.Na + 1) + ia_ + 1;
        const double      progress = static_cast<double>(current) / static_cast<double>(total);
        const std::size_t filled   = static_cast<std::size_t>(progress * kBarWidth);

        os << '\r' << '(' << std::setw(4) << ia_ << ", " << std::setw(4) << ie_ << ") [";

        for (std::size_t i = 0; i < kBarWidth; ++i) {
            os << (i < filled ? '#' : '.');
        }

        os << "] " << std::fixed << std::setw(6) << std::setprecision(2) << progress * 100.0 << " %" << std::flush;
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
std::unique_ptr<T[]> AllocateArray(std::size_t n)
{
    if (n == 0) {
        throw std::invalid_argument("The array size must be greater than zero.");
    }

    return std::make_unique<T[]>(n);
}

// namespace model {
//     /**
//      * @brief Parameters of the planar circular restricted three-body problem.
//      */
//     struct CRTBP2DParams {
//         /**
//          * @brief Gravitational mass parameter.
//          */
//         double mu;
//
//         /**
//          * @brief Constructs the parameter set.
//          *
//          * @param mu Gravitational mass parameter.
//          */
//         explicit CRTBP2DParams(double mu) : mu(mu)
//         {
//         }
//     };
//
//     /**
//      * @brief Computes the right-hand side of the planar circular restricted three-body problem.
//      *
//      * Evaluates the first-order equations of motion in the rotating (synodic)
//      * reference frame. The state vector is defined as
//      * y = (x, y, vx, vy), and the output vector contains the corresponding
//      * time derivatives
//      * dydt = (dx/dt, dy/dt, dvx/dt, dvy/dt).
//      *
//      * @param t Current time (unused, as the equations are autonomous).
//      * @param y Input state vector (x, y, vx, vy).
//      * @param dydt Output time derivatives of the state vector.
//      * @param par Pointer to a @c CRTBP2DParams structure containing the mass parameter.
//      */
//     void CRTBP2Dfun(double t, const double *y, double *dydt, void *par)
//     {
//         (void)t;  // Suppress unused parameter warning.
//
//         const auto  *p  = static_cast<const CRTBP2DParams *>(par);
//         const double mu = p->mu;
//
//         const double r1 = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
//         const double r2 = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));
//
//         const double r1_3 = 1.0 / (r1 * r1 * r1);
//         const double r2_3 = 1.0 / (r2 * r2 * r2);
//
//         // Equations of motion in the rotating frame.
//         dydt[0] = y[2];                                                                                 ///< dx/dt
//         dydt[1] = y[3];                                                                                 ///< dy/dt
//         dydt[2] = 2.0 * y[3] + y[0] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;  ///< dvx/dt
//         dydt[3] = -2.0 * y[2] + y[1] * (1.0 - (1.0 - mu) * r1_3 - mu * r2_3);                           ///< dvy/dt
//     }
//
//     /**
//      * @brief Computes the right-hand side of the planar CRTBP and its variational equations.
//      *
//      * Evaluates the equations of motion and one deviation vector in the rotating
//      * (synodic) reference frame. The state vector is defined as
//      * y = (x, y, vx, vy, dx, dy, dvx, dvy).
//      *
//      * @param t Current time (unused, as the equations are autonomous).
//      * @param y Input state vector and deviation vector.
//      * @param dydt Output derivatives of the state and deviation vectors.
//      * @param par Pointer to a @c CRTBP2DParams structure containing the mass parameter.
//      */
//     void CRTBP2DVariationalFun(double t, const double *y, double *dydt, void *par)
//     {
//         (void)t;  // Suppress unused parameter warning.
//
//         const auto  *p  = static_cast<const CRTBP2DParams *>(par);
//         const double mu = p->mu;
//
//         const double r1 = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
//         const double r2 = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));
//
//         const double r1_3 = 1.0 / (r1 * r1 * r1);
//         const double r2_3 = 1.0 / (r2 * r2 * r2);
//         const double r1_5 = r1_3 / (r1 * r1);
//         const double r2_5 = r2_3 / (r2 * r2);
//
//         // Equations of motion in the rotating frame.
//         dydt[0] = y[2];
//         dydt[1] = y[3];
//         dydt[2] = 2.0 * y[3] + y[0] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;
//         dydt[3] = -2.0 * y[2] + y[1] * (1.0 - (1.0 - mu) * r1_3 - mu * r2_3);
//
//         const double O_xx = 1.0 - (1.0 - mu) * r1_3 - mu * r2_3 + 3.0 * (1.0 - mu) * astro::sqr(y[0] + mu) * r1_5 +
//                             3.0 * mu * astro::sqr(y[0] - 1.0 + mu) * r2_5;
//
//         const double O_xy = 3.0 * (1.0 - mu) * (y[0] + mu) * y[1] * r1_5 + 3.0 * mu * (y[0] - 1.0 + mu) * y[1] *
//         r2_5;
//
//         const double O_yy = 1.0 - (1.0 - mu) * r1_3 - mu * r2_3 + 3.0 * (1.0 - mu) * astro::sqr(y[1]) * r1_5 +
//                             3.0 * mu * astro::sqr(y[1]) * r2_5;
//
//         // Variational equations for one deviation vector.
//         dydt[4] = y[6];
//         dydt[5] = y[7];
//         dydt[6] = O_xx * y[4] + O_xy * y[5] + 2.0 * y[7];
//         dydt[7] = O_xy * y[4] + O_yy * y[5] - 2.0 * y[6];
//     }
// }  // namespace model

namespace ode_integrator {
    void rkf54(var_t &t, StepControl &step, var_t relTol, var_t absTol, Model &model, void *par)
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

        const std::size_t n_var = model.getNVar();
        double           *y_in  = model.getY();

        static std::size_t              allocated_n_var = 0;
        static std::unique_ptr<var_t[]> dy;
        static std::unique_ptr<var_t[]> y;
        static std::unique_ptr<var_t[]> y_out;

        if (n_var != allocated_n_var) {
            dy    = AllocateArray<var_t>(7 * n_var);
            y     = AllocateArray<var_t>(n_var);
            y_out = AllocateArray<var_t>(n_var);

            allocated_n_var = n_var;
        }

        double t0    = t;
        var_t  temax = 0.0;

        model.f(t0, y_in, dy.get(), par);
        do {
            temax = 0.0;
            for (std::size_t k = 1; k < 7; k++) {
                t = t0 + Ci[k] * step.h;

                for (std::size_t n = 0; n < n_var; n++) {
                    y[n] = y_in[n];

                    for (std::size_t l = 0; l < k; l++)
                        y[n] += step.h * Aij[k][l] * dy[l * n_var + n];
                }
                model.f(t, y.get(), dy.get() + k * n_var, par);
            }

            for (std::size_t n = 0; n < n_var; ++n) {
                y_out[n] = y_in[n];

                for (std::size_t k = 0; k < 7; ++k) {
                    y_out[n] += step.h * Bi[k] * dy[k * n_var + n];
                }
                const var_t err = step.h * std::fabs(dy[5 * n_var + n] - dy[6 * n_var + n]) / 60.0;
                const var_t tol = absTol + relTol * std::max(std::fabs(y_in[n]), std::fabs(y_out[n]));

                if (err / tol > temax) {
                    temax = err / tol;
                }
            }
            step.h_did = step.h;
            step.h     = 0.9 * step.h_did * std::pow(1.0 / temax, 1.0 / 5.0);
        } while (temax > 1.0);
        step.h_nxt = step.h;
        // Advance to the end of the accepted integration step.
        t = t0 + step.h_did;
        // Copy the accepted solution to the model state vector.
        std::copy_n(y_out.get(), n_var, model.getY());
    }
}  // namespace ode_integrator

namespace {

    /**
     * @brief Parses command-line arguments.
     *
     * Processes the input and output file options and stores the file name,
     * directory, and absolute path separately.
     *
     * The -i and -o options may specify either a file name in the current
     * working directory or a file name including a directory path.
     *
     * @param argc Number of command-line arguments.
     * @param argv Array of command-line arguments.
     * @param opt Structure receiving the parsed command-line options.
     *
     * @throws std::runtime_error If an option is missing its argument or an
     *         unknown command-line option is encountered.
     */
    void ParseCommandLine(int argc, char *argv[], CommandLineOptions &opt)
    {
        for (int i = 1; i < argc; ++i) {
            const std::string arg(argv[i]);

            if (arg == "-i") {
                if (++i >= argc) {
                    throw std::runtime_error("Missing argument after '-i'.");
                }

                const fs::path p = fs::absolute(fs::path(argv[i]));

                opt.input_file = p.filename().string();
                opt.input_dir  = p.parent_path().string();
                opt.input_path = p.string();

            } else if (arg == "-o") {
                if (++i >= argc) {
                    throw std::runtime_error("Missing argument after '-o'.");
                }

                const fs::path p = fs::absolute(fs::path(argv[i]));

                opt.output_file = p.filename().string();
                opt.output_dir  = p.parent_path().string();
                opt.output_path = p.string();

            } else if (arg == "-h" || arg == "--help") {
                opt.show_help = true;

            } else if (arg == "-v" || arg == "--version") {
                opt.show_version = true;

            } else if (arg == "--verbose") {
                opt.verbose = true;

            } else {
                throw std::runtime_error("Unknown command-line option: " + arg);
            }
        }
    }
    /**
     * @brief Removes all whitespace characters from a string.
     *
     * @param text String to be modified.
     */
    void RemoveSpaces(std::string &text)
    {
        text.erase(std::remove_if(text.begin(), text.end(), ::isspace), text.end());
    }

    /**
     * @brief Converts a string to a run mode.
     *
     * @param text Run mode as text.
     *
     * @return Corresponding run mode.
     *
     * @throws std::runtime_error If the run mode is unknown.
     */
    RunMode ParseRunMode(const std::string &text)
    {
        std::string mode(text);

        std::transform(mode.begin(), mode.end(), mode.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        if (mode == "ORBIT") {
            return RunMode::ORBIT;
        }
        if (mode == "INDICATOR") {
            return RunMode::INDICATOR;
        }
        if (mode == "GRID") {
            return RunMode::GRID;
        }

        throw std::runtime_error("Unknown run mode: " + text);
    }

    /**
     * @brief Converts a string to a chaos-indicator type.
     *
     * @param text Indicator type as text.
     *
     * @return Corresponding indicator type.
     *
     * @throws std::runtime_error If the indicator type is unknown.
     */
    IndicatorType ParseIndicatorType(const std::string &text)
    {
        std::string indicator(text);

        std::transform(indicator.begin(), indicator.end(), indicator.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

        if (indicator == "NONE") {
            return IndicatorType::NONE;
        }
        if (indicator == "FLI") {
            return IndicatorType::FLI;
        }
        if (indicator == "LCI") {
            return IndicatorType::LCI;
        }
        if (indicator == "RLI") {
            return IndicatorType::RLI;
        }

        throw std::runtime_error("Unknown indicator type: " + text);
    }

    /**
     * @brief Parses a single line of the initialization file.
     *
     * Comment lines and empty lines are ignored.
     *
     * @param line Input line.
     * @param data Initialization data.
     *
     * @throws std::runtime_error If the line has an invalid format or
     *         contains an unknown keyword.
     */
    void ParseLine(const std::string &line, InitData &data)
    {
        std::string text = line;

        const std::size_t comment = text.find('#');
        if (comment != std::string::npos) {
            text.erase(comment);
        }

        RemoveSpaces(text);

        if (text.empty()) {
            return;
        }

        const std::size_t pos = text.find('=');

        if (pos == std::string::npos) {
            throw std::runtime_error("Missing '=' in initialization file.");
        }

        const std::string key   = text.substr(0, pos);
        const std::string value = text.substr(pos + 1);

        if (key == "mode") {
            data.run_mode = ParseRunMode(value);
            return;
        }

        if (key == "indicator") {
            data.indicator = ParseIndicatorType(value);
            return;
        }

        std::istringstream is(value);
        if (key == "mu") {
            is >> data.mu;
        } else if (key == "t0") {
            is >> data.t0;
        } else if (key == "T") {
            is >> data.T;
        } else if (key == "output_dt") {
            is >> data.output_dt;
        } else if (key == "a0") {
            is >> data.a0;
        } else if (key == "a1") {
            is >> data.a1;
        } else if (key == "Na") {
            is >> data.Na;
        } else if (key == "e0") {
            is >> data.e0;
        } else if (key == "e1") {
            is >> data.e1;
        } else if (key == "Ne") {
            is >> data.Ne;
        } else if (key == "dy1") {
            is >> data.dy[0];
        } else if (key == "dy2") {
            is >> data.dy[1];
        } else if (key == "dy3") {
            is >> data.dy[2];
        } else if (key == "dy4") {
            is >> data.dy[3];
        } else {
            throw std::runtime_error("Unknown keyword: " + key);
        }

        if (!is || !is.eof()) {
            throw std::runtime_error("Invalid value for '" + key + "'.");
        }
    }

    /**
     * @brief Validates the input data.
     *
     * Checks the consistency of the run mode, indicator type, output interval,
     * integration interval, and grid parameters.
     *
     * @param data Input data to validate.
     *
     * @throws std::runtime_error If the input data are inconsistent.
     */
    void ValidateInitData(const InitData &data)
    {
        if (data.T <= data.t0) {
            throw std::runtime_error("T must be greater than t0.");
        }

        switch (data.run_mode) {
            case RunMode::ORBIT:
                if (data.indicator != IndicatorType::NONE) {
                    throw std::runtime_error("ORBIT mode requires indicator = NONE.");
                }

                if (data.output_dt <= 0.0) {
                    throw std::runtime_error("ORBIT mode requires output_dt > 0.");
                }
                break;

            case RunMode::INDICATOR:
                if (data.indicator == IndicatorType::NONE) {
                    throw std::runtime_error("INDICATOR mode requires an indicator.");
                }

                if (data.output_dt <= 0.0) {
                    throw std::runtime_error("INDICATOR mode requires output_dt > 0.");
                }
                break;

            case RunMode::GRID:
                if (data.indicator == IndicatorType::NONE) {
                    throw std::runtime_error("GRID mode requires an indicator.");
                }

                if (data.Na == 0) {
                    throw std::runtime_error("GRID mode requires Na > 0.");
                }

                if (data.Ne == 0) {
                    throw std::runtime_error("GRID mode requires Ne > 0.");
                }

                if (data.a1 < data.a0) {
                    throw std::runtime_error("GRID mode requires a1 >= a0.");
                }

                if (data.e1 < data.e0) {
                    throw std::runtime_error("GRID mode requires e1 >= e0.");
                }
                break;
        }
    }

    /**
     * @brief Opens the output stream.
     *
     * If an output file is specified on the command line, the file is opened
     * using its full path. Otherwise, the standard output stream is returned.
     *
     * @param opt Command-line options.
     * @param fout Output file stream.
     *
     * @return Pointer to the selected output stream.
     *
     * @throws std::runtime_error If the output file cannot be opened.
     */
    std::ostream *OpenOutputStream(const CommandLineOptions &opt, std::ofstream &fout)
    {
        if (opt.output_path.empty()) {
            return &std::cout;
        }

        fout.open(opt.output_path);

        if (!fout) {
            throw std::runtime_error("Cannot open output file: " + opt.output_path);
        }

        return &fout;
    }

    /**
     * @brief Reads the initialization file.
     *
     * @param file_name Initialization file.
     *
     * @return Initialization data.
     *
     * @throws std::runtime_error If the file cannot be read or
     *         contains invalid data.
     */
    InitData ReadInitData(const std::string &file_name)
    {
        std::ifstream file(file_name);

        if (!file) {
            throw std::runtime_error("Cannot open initialization file.");
        }

        InitData data;

        std::string line;

        while (std::getline(file, line)) {
            ParseLine(line, data);
        }

        ValidateInitData(data);

        return data;
    }

    /**
     * @brief Initializes the state vector at the L4 Lagrange point.
     *
     * Places the massless body at the triangular L4 equilibrium point
     * with zero velocity in the rotating reference frame.
     *
     * @param mu Mass parameter.
     * @param y State vector [x, y, vx, vy].
     */
    void InitializeL4(double mu, double *y)
    {
        y[0] = 0.5 - mu;
        y[1] = std::sqrt(3.0) / 2.0;
        y[2] = 0.0;
        y[3] = 0.0;
    }

    /**
     * @brief Computes the initial state vector.
     *
     * Initializes the orbital state and, if present, the deviation vector.
     *
     * @param mu Mass parameter.
     * @param dy Initial deviation vector.
     * @param a Semimajor axis.
     * @param e Eccentricity.
     * @param y State vector.
     * @param n_var Number of variables in the state vector.
     */
    void GetInitialCondition(double mu, const double *dy, double a, double e, double *y, std::size_t n_var)
    {
        y[0] = a * (1.0 - e) - mu;                                      /// x_0
        y[1] = 0.0;                                                     /// y_0
        y[2] = 0.0;                                                     /// vx_0
        y[3] = std::sqrt((1 - mu) / a * (1.0 + e) / (1.0 - e)) - y[0];  /// vy_0

        if (n_var > 4) {
            y[4] = dy[0];  /// dy_0
            y[5] = dy[1];  /// dy_1
            y[6] = dy[2];  /// dy_2
            y[7] = dy[3];  /// dy_3
        }
    }

    /**
     * @brief Limits the current integration step size.
     *
     * Ensures that the current step does not extend beyond the final
     * integration time. If @c step.h_max is positive, the step size is
     * also limited to the specified maximum value.
     *
     * @param t Current integration time.
     * @param T Final integration time.
     * @param step Step-size control parameters.
     */
    void LimitStep(double t, double T, StepControl &step)
    {
        if (t + step.h > T) {
            step.h = T - t;
        }

        if (step.h_max > 0.0 && step.h > step.h_max) {
            step.h = step.h_max;
        }
    }

    /**
     * @brief Checks whether all elements of a state vector are finite.
     *
     * Returns true if all components of the state vector are finite, false otherwise.
     *
     * @param y State vector.
     * @param n Number of elements in the state vector.
     */
    bool CheckFinite(const double *y, std::size_t n)
    {
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(y[i])) {
                return false;
            }
        }
        return true;
    }
}  // namespace

namespace chaos_indicator {
    /**
     * @brief Computes the Fast Lyapunov Indicator from a deviation vector.
     *
     * The FLI is defined as the maximum logarithmic norm of the deviation vector
     * attained up to the current integration time:
     *
     * FLI(t) = max(FLI_previous, log(||delta y(t)||)).
     *
     * The deviation vector is stored in elements y[4], ..., y[7].
     *
     * @param y State vector containing the orbit and the deviation vector.
     * @param previous_fli Largest FLI value obtained before the current step.
     * @return Updated FLI value.
     *
     * @throws std::domain_error If the deviation-vector norm is zero.
     */
    double ComputeFLI(const double *y, double previous_fli)
    {
        const double norm = std::sqrt(astro::sqr(y[4]) + astro::sqr(y[5]) + astro::sqr(y[6]) + astro::sqr(y[7]));

        return std::max(previous_fli, std::log(norm));
    }

    /**
     * @brief Computes the Lyapunov Characteristic Indicator.
     *
     * The LCI is the finite-time approximation of the largest Lyapunov exponent:
     *
     * LCI(t) = log(||delta y(t)|| / ||delta y(t0)||) / (t - t0).
     *
     * The current deviation vector is stored in elements y[4], ..., y[7].
     *
     * @param t Current integration time.
     * @param t0 Initial integration time.
     * @param y State vector containing the orbit and the current deviation vector.
     * @param initial_deviation Initial deviation vector with four components.
     * @return Current value of the Lyapunov Characteristic Indicator.
     *
     * @throws std::domain_error If the elapsed time or either deviation-vector norm
     *         is zero.
     */
    double ComputeLCI(double t, double t0, const double *y, const double *initial_deviation)
    {
        const double elapsed_time = t - t0;

        if (elapsed_time == 0.0) {
            throw std::domain_error("Cannot compute LCI at the initial time.");
        }

        const double current_norm =
            std::sqrt(astro::sqr(y[4]) + astro::sqr(y[5]) + astro::sqr(y[6]) + astro::sqr(y[7]));

        const double initial_norm = std::sqrt(astro::sqr(initial_deviation[0]) + astro::sqr(initial_deviation[1]) +
                                              astro::sqr(initial_deviation[2]) + astro::sqr(initial_deviation[3]));

        if (current_norm == 0.0 || initial_norm == 0.0) {
            throw std::domain_error("Cannot compute LCI from a zero deviation vector.");
        }

        return std::log(current_norm / initial_norm) / elapsed_time;
    }
}  // namespace chaos_indicator

namespace print {
    /**
     * @brief Prints the program version.
     *
     * Displays the program name together with its version number.
     */
    void Version()
    {
        std::cout << PROGRAM_NAME << " version " << PROGRAM_VERSION << '\n';
    }

    /**
     * @brief Prints the command-line help.
     *
     * Displays the program usage together with the supported
     * command-line options.
     */
    void Help()
    {
        std::cout << "Arch of Chaos\n";
        std::cout << "=============\n\n";

        std::cout << "Usage:\n";
        std::cout << "  archofchaos -i <input file> -o <output file>\n\n";

        std::cout << "Options:\n";
        std::cout << "  -i <file>   Input file.\n";
        std::cout << "  -o <file>   Output file.\n";
        std::cout << "  -h          Display this help message.\n";
        std::cout << "  -v          Display program version.\n\n";

        std::cout << "Examples:\n";
        std::cout << "  archofchaos -i input.txt -o output.txt\n";
        std::cout << "  archofchaos -i data/init.txt -o results/orbit.txt\n";
    }

    /**
     * @brief Prints all input parameters.
     *
     * Displays the command-line options together with the values read from
     * the input file in a human-readable format.
     *
     * @param options Parsed command-line options.
     * @param init Input data read from the initialization file.
     */
    void InputData(const CommandLineOptions &options, const InitData &init, std::ostream &os = std::cout)
    {
        os << '\n';
        os << "============================================================\n";
        os << "Input parameters\n";
        os << "============================================================\n\n";
        options.Print(os);
        os << '\n';
        init.Print(os);
        os << "============================================================\n";
    }

}  // namespace print

int main(int argc, char *argv[])
{
    CommandLineOptions opt;
    try {
        ParseCommandLine(argc, argv, opt);
        if (opt.show_version) {
            print::Version();
            return EXIT_SUCCESS;
        }
        if (opt.show_help) {
            print::Help();
            return EXIT_SUCCESS;
        }

        const InitData init = ReadInitData(opt.input_path);
        if (opt.verbose) {
            print::InputData(opt, init);
        }

        CRTBP2D model(init.mu);
        switch (init.indicator) {
            case IndicatorType::NONE:
                model.setNVar(4);  // Equations of motion: x, y, vx, vy.
                model.setFunction(&Model::fun);
                break;

            case IndicatorType::FLI:
                model.setNVar(8);  // Orbit and variational equations: 4 state + 4 deviation variables.
                model.setFunction(&Model::varfun);
                break;

            case IndicatorType::LCI:
            case IndicatorType::RLI:
                throw std::runtime_error("LCI and RLI indicators are not yet implemented.");

            default:
                throw std::runtime_error("Unknown indicator type.");
        }

        /* If the -o option is specified then the output is written to a file */
        /** Output file stream. */
        std::ofstream fout;
        /** Output stream used by the program. */
        std::ostream *out = OpenOutputStream(opt, fout);

        double     relTol     = 1.0e-6;
        double     absTol     = 1.0e-10;
        const auto start_time = std::chrono::steady_clock::now();

        switch (init.run_mode) {
            case RunMode::ORBIT: {
                double t = init.t0;

                StepControl step;
                step.h     = 0.1;
                step.h_max = 0.1;
                step.h_nxt = step.h;
                step.n_tst = 0;
                step.n_int = 0;

                // GetInitialCondition(model.getParams()->mu, init.dy, init.a0, init.e0, model.getY(), model.getNVar());
                InitializeL4(model.getParams()->mu, model.getY());  // Initialize at L4 point
                model.getY()[0] += 2.0e-2;
                // Save the initial state.
                model.printState(*out, t, model.getY());
                double next_output = init.t0 + init.output_dt;

                while (t < init.T - TIME_EPS) {
                    if (step.n_tst % 10 == 0) {
                        if (!CheckFinite(model.getY(), model.getNVar())) {
                            throw std::runtime_error("Non-finite state encountered during orbit integration.");
                        }
                    }

                    const double target_time = std::min(next_output, init.T);

                    // Force the integrator to stop exactly at the next output time.
                    LimitStep(t, target_time, step);

                    ode_integrator::rkf54(t, step, relTol, absTol, model, model.getParams());

                    ++step.n_int;
                    ++step.n_tst;

                    if (t >= next_output - TIME_EPS) {
                        model.printState(*out, t, model.getY());
                        next_output += init.output_dt;
                    }
                } /* while */

                break;
            }
            case RunMode::INDICATOR: {
                // ...
                break;
            }

            case RunMode::GRID: {
                // ...
                break;
            }
        }

        // do {
        //     t = init.t0;
        //     StepControl step;  // Step control structure for adaptive integration
        //     step.h     = 0.1;
        //     step.h_max = 0.1;
        //     step.h_nxt = step.h;
        //     double fli = 0.0;  // Initialize the Fast Lyapunov Indicator

        //    double a = grid.a();
        //    double e = grid.e();
        //    // Get initial conditions for the current (a,e)
        //    GetInitialCondition(model.getParams()->mu, init.dy, a, e, model.getY(), model.getNVar());
        //    // InitializeL4(param, y_in.get());  // Initialize at L4 point)
        //    // y_in[0] += 3.0e-2;
        //    // y_in[4] = 1.0;

        //    step.n_tst = 0;
        //    step.n_int = 0;
        //    grid.PrintProgress();

        //    do {
        //        if (step.n_tst % 10 == 0) {
        //            if (!CheckFinite(model.getY(), model.getNVar())) {
        //                std::cerr << " Error at grid point. Proceed to the next grid point.\n";
        //                break;
        //            }
        //        }

        //        LimitStep(t, init.T, step);  // Limit the step size to not exceed the final time

        //        ode_integrator::rkf54(t, step, relTol, absTol, model, model.getParams());

        //        switch (init.compute) {
        //            case ComputeMode::ORBIT:
        //                model.printState(*out, t, model.getY());  // Print the state (time and variables)
        //                break;
        //            case ComputeMode::FLI:
        //                fli = chaos_indicator::ComputeFLI(model.getY(), fli);
        //                break;
        //            case ComputeMode::LCI:
        //            case ComputeMode::RLI:
        //                throw std::runtime_error("LCI and RLI modes are not yet implemented.");
        //                break;
        //            default:
        //                throw std::runtime_error("Unknown computation mode.");
        //        }
        //        step.n_int++;
        //        step.n_tst++;

        //    } while (t < init.T);  // std::abs(init.T - t) > TIME_EPS);
        //    print::State(*out, grid, init, t, model.getY(), model.getNVar(),
        //                 fli);  // Print the state (time and variables)
        //} while (grid.Next());
        const auto                          end_time     = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_time = end_time - start_time;
        std::cout << "\n" << "Total runtime : " << elapsed_time.count() << " s\n";

        return EXIT_SUCCESS;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }

    return EXIT_FAILURE;
}
