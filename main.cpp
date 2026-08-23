#include "init_data.h"   // InitData, RunMode, IndicatorType
#include "math_utils.h"  // astro::sqr
#include "model.h"       // Model base class
#include "orbit.h"       // astro::orbit::elementsToState
#include "crtbp.h"

#include <algorithm>  /**< std::copy, std::copy_n, std::remove_if. */
#include <chrono>     /**< Provides time measurement utilities. */
#include <cmath>      /**< Mathematical functions (std::sqrt, std::abs). */
#include <cctype>     /**< Character classification functions (std::isspace). */
#include <cstddef>    /**< std::size_t. */
#include <cstdint>    /**< Fixed-width integer types (std::uint32_t). */
#include <exception>  /**< std::exception base class. */
#include <filesystem> /**< std::filesystem::path for portable path handling. */
#include <fstream>    /**< File stream classes (std::ifstream, std::ofstream). */
#include <iomanip>    /**< Output manipulators (std::fixed, std::setprecision, std::setw). */
#include <iostream>   /**< Standard input/output streams (std::cout, std::cerr). */
#include <limits>     /**< std::numeric_limits  */
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
     */
    GridIterator(double a0, double a1, std::uint32_t Na, double e0, double e1, std::uint32_t Ne)
        : a0_(a0),
          e0_(e0),
          Na_(Na),
          Ne_(Ne),
          ia_(0),
          ie_(0),
          da_((a1 - a0) / static_cast<double>(Na)),
          de_((e1 - e0) / static_cast<double>(Ne))
    {
    }

    /**
     * @brief Returns the current semimajor axis.
     *
     * @return Current value of a.
     */
    double a() const noexcept
    {
        return a0_ + static_cast<double>(ia_) * da_;
    }

    /**
     * @brief Returns the current eccentricity.
     *
     * @return Current value of e.
     */
    double e() const noexcept
    {
        return e0_ + static_cast<double>(ie_) * de_;
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
    bool Next() noexcept
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

    /**
     * @brief Prints the current grid progress.
     *
     * @param os Output stream.
     */
    void PrintProgress(std::ostream &os = std::cout) const
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

namespace ode_integrator {
    void rkf54(Model &model, void *par, StepControl &step, var_t relTol, var_t absTol)
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

        double t0    = model.getT();
        var_t  temax = 0.0;

        model.f(t0, y_in, dy.get(), par);
        do {
            temax = 0.0;
            for (std::size_t k = 1; k < 7; k++) {
                double t = t0 + Ci[k] * step.h;

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
        model.setT(t0 + step.h_did);
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

        if (norm == 0.0) {
            throw std::domain_error("Cannot compute FLI from a zero deviation vector.");
        }

        return std::max(previous_fli, std::log(norm));
    }

    /**
     * @brief Computes the Lyapunov Characteristic Indicator (LCI).
     *
     * The current deviation vector is stored in elements y[4], ..., y[7].
     * The initial deviation vector is stored in elements y0[0], ..., y0[3].
     *
     * The LCI is calculated as
     *
     * \f[
     * \mathrm{LCI}(t) =
     * \frac{1}{t-t_0}
     * \ln\left(
     * \frac{\|\delta(t)\|}{\|\delta(t_0)\|}
     * \right).
     * \f]
     *
     * @param t Current integration time.
     * @param t0 Initial integration time.
     * @param y State vector containing the orbit and the current deviation vector.
     * @param y0 Initial deviation vector with four components.
     *
     * @return Current value of the Lyapunov Characteristic Indicator.
     *
     * @throws std::domain_error If the elapsed time is zero or if either
     *         deviation-vector norm is zero.
     */
    double ComputeLCI(double t, double t0, const double *y, const double *y0)
    {
        const double elapsed_time = t - t0;

        if (elapsed_time == 0.0) {
            throw std::domain_error("Cannot compute LCI at the initial time.");
        }

        const double initial_norm =
            std::sqrt(astro::sqr(y0[0]) + astro::sqr(y0[1]) + astro::sqr(y0[2]) + astro::sqr(y0[3]));

        const double current_norm =
            std::sqrt(astro::sqr(y[4]) + astro::sqr(y[5]) + astro::sqr(y[6]) + astro::sqr(y[7]));

        if (initial_norm == 0.0 || current_norm == 0.0) {
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

StepControl createStepControl()
{
    StepControl step;

    step.h     = 0.1;
    step.h_max = 0.1;
    step.h_nxt = step.h;
    step.n_tst = 0;
    step.n_int = 0;

    return step;
}

void runOrbit(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol, double absTol,
              double mu_13, double n)
{
    // Orbital elements -> heliocentric inertial Cartesian state.
    const astro::State state = astro::calcState(mu_13, init.getT0(), init.getElements());
    // Heliocentric inertial state -> normalized rotating CRTBP state.
    model.InertialToCRTBP(state, init.getA2(), n);
    // Save the initial state.
    model.printState(out, model.getT(), model.getY());

    double next_output = init.getT0() + init.getOutputDt();
    if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
        model.VelocityToHamiltonian();
    }

    while (model.getT() < init.getT() - TIME_EPS) {
        if (step.n_tst % 10 == 0) {
            if (!CheckFinite(model.getY(), model.getNVar())) {
                throw std::runtime_error("Non-finite state encountered during orbit integration.");
            }
        }

        const double target_time = std::min(next_output, init.getT());
        // Force the integrator to stop exactly at the next output time.
        LimitStep(model.getT(), target_time, step);

        ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);

        ++step.n_int;
        ++step.n_tst;

        const bool output_time = model.getT() >= next_output - TIME_EPS;
        const bool final_time = model.getT() >= init.getT() - TIME_EPS;
        if (output_time || final_time) {
            if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
                double y_out[4];

                model.HamiltonianToNewtonian(y_out);
                model.printState(out, model.getT(), y_out);
            } else {
                model.printState(out, model.getT(), model.getY());
            }
            if (output_time) {
                next_output += init.getOutputDt();
            }
        }
    }
}

void runIndicator(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol,
                  double absTol, double mu_13, double n)
{
    constexpr int W = 18;

    // Orbital elements -> heliocentric inertial Cartesian state.
    const astro::State state = astro::calcState(mu_13, init.getT0(), init.getElements());
    // Heliocentric inertial state -> normalized rotating CRTBP state.
    model.InertialToCRTBP(state, init.getA2(), n);
    // Convert the orbital state to Hamiltonian canonical variables if required.
    if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
        model.VelocityToHamiltonian();
    }

    // Set the initial deviation vector.
    std::copy_n(init.getDy(), 4, model.getY() + 4);
    // Current indicator value.
    double indicator_value = 0.0;

    switch (model.getIndicator()) {
        case Model::IndicatorType::FLI:
            indicator_value = chaos_indicator::ComputeFLI(model.getY(), 1.0);

            out << std::left << std::setw(W) << "t" << std::setw(W) << "FLI" << '\n';
            out << std::right << std::scientific << std::setprecision(10);
            // FLI is defined at the initial time.
            out << std::setw(W) << model.getT() << std::setw(W) << indicator_value << '\n';
            break;

        case Model::IndicatorType::LCI:
            out << std::left << std::setw(W) << "t" << std::setw(W) << "LCI" << '\n';
            out << std::right << std::scientific << std::setprecision(10);
            // LCI is not defined at the initial time.
            break;

        case Model::IndicatorType::RLI:
            throw std::runtime_error("RLI indicator is not yet implemented.");

        case Model::IndicatorType::NONE:
            throw std::runtime_error("INDICATOR mode requires a chaos indicator.");

        default:
            throw std::runtime_error("Unknown chaos indicator.");
    }

    double next_output = init.getT0() + init.getOutputDt();

    while (model.getT() < init.getT() - TIME_EPS) {
        // Check the numerical state periodically.
        if (step.n_tst % 10 == 0) {
            if (!CheckFinite(model.getY(), model.getNVar())) {
                throw std::runtime_error("Non-finite state encountered during indicator integration.");
            }
        }

        // Stop exactly at the next output time or at the final time.
        const double target_time = std::min(next_output, init.getT());

        LimitStep(model.getT(), target_time, step);

        // Integrate the orbit and the variational equations.
        ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);

        // Update indicators that depend on all intermediate integration steps.
        switch (model.getIndicator()) {
            case Model::IndicatorType::FLI:
                indicator_value = chaos_indicator::ComputeFLI(model.getY(), indicator_value);
                break;

            case Model::IndicatorType::LCI:
                indicator_value = chaos_indicator::ComputeLCI(model.getT(), init.getT0(), model.getY(), init.getDy());
                break;

            default:
                break;
        }

        ++step.n_int;
        ++step.n_tst;

        const bool output_time = model.getT() >= next_output - TIME_EPS;
        const bool final_time = model.getT() >= init.getT() - TIME_EPS;
        if (output_time || final_time) {
            out << std::setw(W) << model.getT() << std::setw(W) << indicator_value << '\n';
            if (output_time) {
                next_output += init.getOutputDt();
            }
        }
    }
}

void runGrid(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol, double absTol,
             double mu_13, double n)
{
    // TODO: Move the current GRID branch implementation here.
}

void run(const InitData &init, const CommandLineOptions &opt)
{
    constexpr double relTol = 1.0e-6;
    constexpr double absTol = 1.0e-10;

    const double mu = init.getM2() / (init.getM1() + init.getM2());
    const double mu_12 = astro::sqr(astro::k) * (init.getM1() + init.getM2());
    const double mu_13 = astro::sqr(astro::k) * init.getM1();
    const double n = std::sqrt(mu_12 / astro::cube(init.getA2()));

    CRTBP2D model(init.getT0(), mu, init.getFormalism(), init.getIndicator());
    StepControl step = createStepControl();

    std::ofstream fout;
    std::ostream* out = OpenOutputStream(opt, fout);

    switch (init.getRunMode()) {
        case RunMode::ORBIT:
            runOrbit(model, init, *out, step, relTol, absTol, mu_13, n);
            break;

        case RunMode::INDICATOR:
            runIndicator(model, init, *out, step, relTol, absTol, mu_13, n);
            break;

        case RunMode::GRID:
            runGrid(model, init, *out, step, relTol, absTol, mu_13, n);
            break;

        default:
            throw std::runtime_error("Unknown run mode.");
    }
}

int main(int argc, char *argv[])
{
    try {
        CommandLineOptions opt;
        ParseCommandLine(argc, argv, opt);

        if (opt.show_version) {
            print::Version();
            return EXIT_SUCCESS;
        }
        if (opt.show_help) {
            print::Help();
            return EXIT_SUCCESS;
        }

        const InitData init(opt.input_path);

        if (opt.verbose) {
            init.Print(std::cout);
        }

        const auto start_time = std::chrono::steady_clock::now();
        run(init, opt);
        const auto end_time = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_time = end_time - start_time;
        std::cout << "\nTotal runtime : " << elapsed_time.count() << " s\n";
        return EXIT_SUCCESS;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }

    return EXIT_FAILURE;
}

// int main(int argc, char *argv[])
//{
//     CommandLineOptions opt;
//     try {
//         ParseCommandLine(argc, argv, opt);
//         if (opt.show_version) {
//             print::Version();
//             return EXIT_SUCCESS;
//         }
//         if (opt.show_help) {
//             print::Help();
//             return EXIT_SUCCESS;
//         }
//
//         const InitData init(opt.input_path);
//         if (opt.verbose) {
//             init.Print(std::cout);
//         }
//
//         const double mu = init.getM2() / (init.getM1() + init.getM2());
//         CRTBP2D      model(init.getT0(), mu, init.getFormalism(), init.getIndicator());
//
//         /* If the -o option is specified then the output is written to a file */
//         /** Output file stream. */
//         std::ofstream fout;
//         /** Output stream used by the program. */
//         std::ostream *out = OpenOutputStream(opt, fout);
//
//         double     relTol     = 1.0e-6;
//         double     absTol     = 1.0e-10;
//         const auto start_time = std::chrono::steady_clock::now();
//
//         StepControl step;
//         step.h     = 0.1;
//         step.h_max = 0.1;
//         step.h_nxt = step.h;
//         step.n_tst = 0;
//         step.n_int = 0;
//
//         // Gravitational parameter of the P1-P2 relative orbit.
//         const double mu_12 = astro::sqr(astro::k) * (init.getM1() + init.getM2());
//         // Gravitational parameter of the P1-P3 heliocentric orbit.
//         const double mu_13 = astro::sqr(astro::k) * init.getM1();
//         // Mean motion of the P1-P2 system.
//         const double n = std::sqrt(mu_12 / astro::cube(init.getA2()));
//
//         switch (init.getRunMode()) {
//             case RunMode::ORBIT: {
//                 // Orbital elements -> heliocentric inertial Cartesian state.
//                 const astro::State state = astro::calcState(mu_13, init.getT0(), init.getElements());
//                 // Heliocentric inertial state -> normalized rotating CRTBP state.
//                 model.InertialToCRTBP(state, init.getA2(), n);
//
//                 // Save the initial state.
//                 model.printState(*out, model.getT(), model.getY());
//                 double next_output = init.getT0() + init.getOutputDt();
//
//                 if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
//                     model.VelocityToHamiltonian();
//                 }
//
//                 while (model.getT() < init.getT() - TIME_EPS) {
//                     if (step.n_tst % 10 == 0) {
//                         if (!CheckFinite(model.getY(), model.getNVar())) {
//                             throw std::runtime_error("Non-finite state encountered during orbit integration.");
//                         }
//                     }
//                     const double target_time = std::min(next_output, init.getT());
//                     // Force the integrator to stop exactly at the next output time.
//                     LimitStep(model.getT(), target_time, step);
//
//                     ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);
//
//                     ++step.n_int;
//                     ++step.n_tst;
//
//                     const bool output_time = model.getT() >= next_output - TIME_EPS;
//                     const bool final_time  = model.getT() >= init.getT() - TIME_EPS;
//
//                     if (output_time || final_time) {
//                         if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
//                             double y_out[4];
//
//                             model.HamiltonianToNewtonian(y_out);
//                             model.printState(*out, model.getT(), y_out);
//                         } else {
//                             model.printState(*out, model.getT(), model.getY());
//                         }
//                         next_output += init.getOutputDt();
//                     }
//                 } /* while */
//
//                 break;
//             } /* ORBIT     */
//             case RunMode::INDICATOR: {
//                 // Orbital elements -> heliocentric inertial Cartesian state.
//                 const astro::State state = astro::calcState(mu_13, init.getT0(), init.getElements());
//
//                 // Heliocentric inertial state -> normalized rotating CRTBP state.
//                 model.InertialToCRTBP(state, init.getA2(), n);
//
//                 // Convert the orbital state to Hamiltonian canonical variables if required.
//                 if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
//                     model.VelocityToHamiltonian();
//                 }
//
//                 // Set the initial deviation vector.
//                 std::copy_n(init.getDy(), 4, model.getY() + 4);
//
//                 constexpr int W = 18;
//
//                 // Current indicator value.
//                 double indicator_value = 0.0;
//
//                 switch (model.getIndicator()) {
//                     case Model::IndicatorType::FLI:
//                         indicator_value = chaos_indicator::ComputeFLI(model.getY(), 1.0);
//                         *out << std::left << std::setw(W) << "t" << std::setw(W) << "FLI" << '\n';
//                         *out << std::right << std::scientific << std::setprecision(10);
//                         // FLI is defined at the initial time.
//                         *out << std::setw(W) << model.getT() << std::setw(W) << indicator_value << '\n';
//                         break;
//
//                     case Model::IndicatorType::LCI:
//                         *out << std::left << std::setw(W) << "t" << std::setw(W) << "LCI" << '\n';
//                         *out << std::right << std::scientific << std::setprecision(10);
//                         // LCI is not defined at the initial time.
//                         break;
//
//                     case Model::IndicatorType::RLI:
//                         throw std::runtime_error("RLI indicator is not yet implemented.");
//
//                     case Model::IndicatorType::NONE:
//                         throw std::runtime_error("INDICATOR mode requires a chaos indicator.");
//
//                     default:
//                         throw std::runtime_error("Unknown chaos indicator.");
//                 }
//
//                 double next_output = init.getT0() + init.getOutputDt();
//
//                 while (model.getT() < init.getT() - TIME_EPS) {
//                     // Check the numerical state periodically.
//                     if (step.n_tst % 10 == 0) {
//                         if (!CheckFinite(model.getY(), model.getNVar())) {
//                             throw std::runtime_error("Non-finite state encountered during indicator integration.");
//                         }
//                     }
//
//                     // Stop exactly at the next output time or at the final time.
//                     const double target_time = std::min(next_output, init.getT());
//
//                     LimitStep(model.getT(), target_time, step);
//
//                     // Integrate the orbit and the variational equations.
//                     ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);
//
//                     // Update indicators that depend on all intermediate integration steps.
//                     switch (model.getIndicator()) {
//                         case Model::IndicatorType::FLI:
//                             indicator_value = chaos_indicator::ComputeFLI(model.getY(), indicator_value);
//                             break;
//
//                         case Model::IndicatorType::LCI:
//                             indicator_value =
//                                 chaos_indicator::ComputeLCI(model.getT(), init.getT0(), model.getY(), init.getDy());
//                             break;
//
//                         default:
//                             break;
//                     }
//
//                     ++step.n_int;
//                     ++step.n_tst;
//
//                     const bool output_time = model.getT() >= next_output - TIME_EPS;
//                     const bool final_time  = model.getT() >= init.getT() - TIME_EPS;
//                     if (output_time || final_time) {
//                         *out << std::setw(W) << model.getT() << std::setw(W) << indicator_value << '\n';
//                         if (output_time) {
//                             next_output += init.getOutputDt();
//                         }
//                     }
//                 }
//
//                 break;
//             } /* INDICATOR */
//             case RunMode::GRID: {
//                 // ...
//                 break;
//             } /* GRID      */
//         }
//
//         const auto                          end_time     = std::chrono::steady_clock::now();
//         const std::chrono::duration<double> elapsed_time = end_time - start_time;
//         std::cout << "\n" << "Total runtime : " << elapsed_time.count() << " s\n";
//
//         return EXIT_SUCCESS;
//     } catch (const std::exception &e) {
//         std::cerr << "Error: " << e.what() << '\n';
//     } catch (...) {
//         std::cerr << "Unknown error.\n";
//     }
//
//     return EXIT_FAILURE;
// }
