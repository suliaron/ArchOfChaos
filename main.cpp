#include "crtbp.h"       // CRTBP2D model
#include "grid.h"        // GridIterator
#include "init_data.h"   // InitData, RunMode
#include "math_utils.h"  // astro::sqr, astro::cube
#include "model.h"       // Model base class
#include "orbit.h"       // Orbital-element/state transformations

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
    void print(std::ostream &os = std::cout) const
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
std::unique_ptr<T[]> allocateArray(std::size_t n)
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
            dy    = allocateArray<var_t>(7 * n_var);
            y     = allocateArray<var_t>(n_var);
            y_out = allocateArray<var_t>(n_var);

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
    void parseCommandLine(int argc, char *argv[], CommandLineOptions &opt)
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
    std::ostream *openOutputStream(const CommandLineOptions &opt, std::ofstream &fout)
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
    void limitStep(double t, double T, StepControl &step)
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
    bool checkFinite(const double *y, std::size_t n)
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
    double computeFLI(const double *y, double previous_fli)
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
    double computeLCI(double t, double t0, const double *y, const double *y0)
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
    void version()
    {
        std::cout << PROGRAM_NAME << " version " << PROGRAM_VERSION << '\n';
    }

    /**
     * @brief Prints the command-line help.
     *
     * Displays the program usage together with the supported
     * command-line options.
     */
    void help()
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
    void inputData(const CommandLineOptions &options, const InitData &init, std::ostream &os = std::cout)
    {
        os << '\n';
        os << "============================================================\n";
        os << "Input parameters\n";
        os << "============================================================\n\n";
        options.print(os);
        os << '\n';
        init.print(os);
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

/**
 * @brief Integrates a single CRTBP orbit and writes its state evolution.
 *
 * Constructs the initial heliocentric inertial state from the input orbital
 * elements, transforms it to dimensionless rotating CRTBP coordinates, and
 * integrates the equations of motion over the specified dimensionless time.
 *
 * The numerical integration uses dimensionless CRTBP time, while the time
 * written to the output is expressed in physical days.
 *
 * If enabled, an additional test output file is generated containing the
 * state transformed back to the P1-centered inertial reference frame.
 *
 * @param model CRTBP model.
 * @param init Initialization data.
 * @param out Output stream for the CRTBP orbit.
 * @param step Adaptive integration step-control parameters.
 * @param relTol Relative integration tolerance.
 * @param absTol Absolute integration tolerance.
 * @param mu13 Gravitational parameter of the P1-P3 heliocentric orbit
 *             [AU^3/day^2].
 * @param n Mean motion of the P1-P2 system [rad/day].
 * @param tDimless Dimensionless integration duration.
 * @param outputDtDimless Dimensionless output time interval.
 *
 * @throws std::runtime_error If a non-finite state is encountered or the
 *         inertial test-output file cannot be opened.
 */
void runOrbit(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol, double absTol,
              double mu13, double n, double tDimless, double outputDtDimless)
{
    // Copy the input orbital elements.
    astro::OrbitalElements elements = init.getElements();
    // Calculate the pericenter passage time from the selected
    // orbital-phase input.
    elements.tau = init.calc_tau(mu13, elements.a);
    // Orbital elements -> heliocentric inertial Cartesian state.
    const astro::State state = astro::calcState(mu13, init.getT0(), elements);
    // Heliocentric inertial state -> dimensionless rotating CRTBP state.
    model.inertialToCRTBP(state, init.getA2(), n);

#if 1  // Test CRTBP -> inertial transformation
    std::ofstream inertialOut("output_inertial.txt");
    if (!inertialOut) {
        throw std::runtime_error("Cannot open output_inertial.txt.");
    }
    inertialOut << "# t [day]            x [AU]             y [AU]"
                   "             vx [AU/day]         vy [AU/day]\n";
    inertialOut << std::scientific << std::setprecision(10) << std::showpos;
    // Transform the initial CRTBP state back to the
    // P1-centered inertial reference frame.
    const astro::State initialInertial = model.crtbpToInertial(model.getY(), init.getA2(), n);
    inertialOut << init.getT0() << ' ' << initialInertial.r.x << ' ' << initialInertial.r.y << ' '
                << initialInertial.v.x << ' ' << initialInertial.v.y << '\n';
#endif

    // Save the initial state using physical time [day].
    model.printState(out, init.getT0(), model.getY());

    // Convert the orbital state to Hamiltonian canonical variables if required.
    if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
        model.velocityToHamiltonian();
    }

    // First output epoch in dimensionless CRTBP time.
    double nextOutputDimless = outputDtDimless;

    while (model.getT() < tDimless - TIME_EPS) {
        if (step.n_tst % 10 == 0) {
            if (!checkFinite(model.getY(), model.getNVar())) {
                throw std::runtime_error("Non-finite state encountered during orbit integration.");
            }
        }

        const double targetTime = std::min(nextOutputDimless, tDimless);
        // Force the integrator to stop exactly at the next output
        // or final dimensionless time.
        limitStep(model.getT(), targetTime, step);

        ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);

        ++step.n_int;
        ++step.n_tst;

        const bool outputTimeReached = model.getT() >= nextOutputDimless - TIME_EPS;
        const bool finalTimeReached = model.getT() >= tDimless - TIME_EPS;

        if (outputTimeReached || finalTimeReached) {
            const double tDay = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);

            if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
                double y_out[4];

                // Hamiltonian -> Newtonian state.
                model.hamiltonianToNewtonian(y_out);

                // Save the Newtonian CRTBP state.
                model.printState(out, tDay, y_out);

#if 1  // Test CRTBP -> inertial transformation
                const astro::State inertialState = model.crtbpToInertial(y_out, init.getA2(), n);
                inertialOut << tDay << ' ' << inertialState.r.x << ' ' << inertialState.r.y << ' ' << inertialState.v.x
                            << ' ' << inertialState.v.y << '\n';
#endif

            } else {
                // Save the Newtonian CRTBP state.
                model.printState(out, tDay, model.getY());

#if 1  // Test CRTBP -> inertial transformation
                const astro::State inertialState = model.crtbpToInertial(model.getY(), init.getA2(), n);
                inertialOut << tDay << ' ' << inertialState.r.x << ' ' << inertialState.r.y << ' ' << inertialState.v.x
                            << ' ' << inertialState.v.y << '\n';
#endif
            }

            if (outputTimeReached) {
                nextOutputDimless += outputDtDimless;
            }
        }
    } /* while */
}

// void runOrbit(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol, double
// absTol,
//               double mu_13, double n)
//{
//     // Copy the input orbital elements.
//     astro::OrbitalElements elements = init.getElements();
//     // Calculate the pericenter passage time from the selected
//     // orbital-phase input (tau or mean anomaly M).
//     elements.tau = init.calc_tau(mu_13, elements.a);
//
//     // Diagnosztikaként ideiglenesen közvetlenül elé tenném :
//     std::cout << std::scientific << std::setprecision(15) << "Calculated tau = " << elements.tau << '\n';
//
//     // Orbital elements -> heliocentric inertial Cartesian state.
//     const astro::State state = astro::calcState(mu_13, init.getT0(), elements);
//
//     // Heliocentric inertial state -> normalized rotating CRTBP state.
//     model.inertialToCRTBP(state, init.getA2(), n);
//     // Save the initial state.
//     model.printState(out, model.getT(), model.getY());
//
//     double next_output = init.getT0() + init.getOutputDt();
//     if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
//         model.velocityToHamiltonian();
//     }
//
//     while (model.getT() < init.getT() - TIME_EPS) {
//         if (step.n_tst % 10 == 0) {
//             if (!checkFinite(model.getY(), model.getNVar())) {
//                 throw std::runtime_error("Non-finite state encountered during orbit integration.");
//             }
//         }
//
//         const double target_time = std::min(next_output, init.getT());
//         // Force the integrator to stop exactly at the next output time.
//         limitStep(model.getT(), target_time, step);
//
//         ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);
//
//         ++step.n_int;
//         ++step.n_tst;
//
//         const bool output_time = model.getT() >= next_output - TIME_EPS;
//         const bool final_time  = model.getT() >= init.getT() - TIME_EPS;
//         if (output_time || final_time) {
//             if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
//                 double y_out[4];
//
//                 model.hamiltonianToNewtonian(y_out);
//                 model.printState(out, model.getT(), y_out);
//             } else {
//                 model.printState(out, model.getT(), model.getY());
//             }
//             if (output_time) {
//                 next_output += init.getOutputDt();
//             }
//         }
//     }
// }

void runIndicator(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol,
                  double absTol, double mu_13, double n)
{
    constexpr int W = 18;

    // Orbital elements -> heliocentric inertial Cartesian state.
    const astro::State state = astro::calcState(mu_13, init.getT0(), init.getElements());
    // Heliocentric inertial state -> normalized rotating CRTBP state.
    model.inertialToCRTBP(state, init.getA2(), n);
    // Convert the orbital state to Hamiltonian canonical variables if required.
    if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
        model.velocityToHamiltonian();
    }

    // Set the initial deviation vector.
    std::copy_n(init.getDy(), 4, model.getY() + 4);
    // Current indicator value.
    double indicator_value = 0.0;

    switch (model.getIndicator()) {
        case Model::IndicatorType::FLI:
            indicator_value = chaos_indicator::computeFLI(model.getY(), 1.0);

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
            if (!checkFinite(model.getY(), model.getNVar())) {
                throw std::runtime_error("Non-finite state encountered during indicator integration.");
            }
        }

        // Stop exactly at the next output time or at the final time.
        const double target_time = std::min(next_output, init.getT());

        limitStep(model.getT(), target_time, step);

        // Integrate the orbit and the variational equations.
        ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);

        // Update indicators that depend on all intermediate integration steps.
        switch (model.getIndicator()) {
            case Model::IndicatorType::FLI:
                indicator_value = chaos_indicator::computeFLI(model.getY(), indicator_value);
                break;

            case Model::IndicatorType::LCI:
                indicator_value = chaos_indicator::computeLCI(model.getT(), init.getT0(), model.getY(), init.getDy());
                break;

            default:
                break;
        }

        ++step.n_int;
        ++step.n_tst;

        const bool output_time = model.getT() >= next_output - TIME_EPS;
        const bool final_time  = model.getT() >= init.getT() - TIME_EPS;
        if (output_time || final_time) {
            out << std::setw(W) << model.getT() << std::setw(W) << indicator_value << '\n';
            if (output_time) {
                next_output += init.getOutputDt();
            }
        }
    }
}

/**
 * @brief Integrates a two-dimensional (a,e) grid and computes a chaos indicator.
 *
 * For each grid point, the semimajor axis and eccentricity are inserted into
 * the osculating orbital elements. The corresponding heliocentric inertial
 * state is transformed to the normalized rotating CRTBP frame and integrated
 * up to the final time.
 *
 * The initial deviation vector is copied directly from the input data in both
 * Newtonian and Hamiltonian formulations. In the Hamiltonian formulation only
 * the orbital state is transformed to canonical variables.
 *
 * The final FLI or LCI value is written for every grid point.
 *
 * @param model CRTBP model.
 * @param init Initialization data.
 * @param out Output stream.
 * @param step Adaptive integration step control.
 * @param relTol Relative integration tolerance.
 * @param absTol Absolute integration tolerance.
 * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit.
 * @param n Mean motion of the P1-P2 system.
 */
void runGrid(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double relTol, double absTol,
             double mu_13, double n)
{
    constexpr int W = 18;

    GridIterator grid(init.getA0(), init.getA1(), init.getNa(), init.getE0(), init.getE1(), init.getNe());
    // Check the selected chaos indicator.
    switch (model.getIndicator()) {
        case Model::IndicatorType::FLI:
        case Model::IndicatorType::LCI:
            break;

        case Model::IndicatorType::RLI:
            throw std::runtime_error("RLI indicator is not yet implemented.");

        case Model::IndicatorType::NONE:
            throw std::runtime_error("GRID mode requires a chaos indicator.");

        default:
            throw std::runtime_error("Unknown chaos indicator.");
    }

    // Write the output header.
    out << std::left << std::setw(W) << "a" << std::setw(W) << "e" << std::setw(W)
        << Model::indicatorTypeToString(model.getIndicator()) << '\n';

    out << std::right << std::scientific << std::setprecision(10);

    // Save the initial step-control values.
    const double initial_h     = step.h;
    const double initial_h_max = step.h_max;
    const double initial_h_min = step.h_min;

    do {
        // Reset the model time.
        model.setT(init.getT0());

        // Reset the adaptive step-size control.
        step.h     = initial_h;
        step.h_nxt = initial_h;
        step.h_did = 0.0;
        step.h_max = initial_h_max;
        step.h_min = initial_h_min;
        step.n_tst = 0;
        step.n_int = 0;

        const double a = grid.a();
        const double e = grid.e();

        // Construct the orbital elements for the current grid point.
        astro::OrbitalElements elements = init.getElements();

        elements.a = a;
        elements.e = e;

        // Orbital elements -> heliocentric inertial Cartesian state.
        const astro::State state = astro::calcState(mu_13, init.getT0(), elements);
        // Heliocentric inertial state -> normalized rotating CRTBP state.
        model.inertialToCRTBP(state, init.getA2(), n);
        // Convert only the orbital state to Hamiltonian canonical variables.
        if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
            model.velocityToHamiltonian();
        }
        // Set the initial deviation vector exactly as specified
        // in the input file.
        std::copy_n(init.getDy(), 4, model.getY() + 4);

        double indicator_value = 0.0;
        // Initialize the FLI.
        if (model.getIndicator() == Model::IndicatorType::FLI) {
            indicator_value = chaos_indicator::computeFLI(model.getY(), 1.0);
        }
        bool valid = true;

        // Integrate the current grid point.
        while (model.getT() < init.getT() - TIME_EPS) {
            // Check the numerical state periodically.
            if (step.n_tst % 10 == 0) {
                if (!checkFinite(model.getY(), model.getNVar())) {
                    valid = false;
                    break;
                }
            }

            // Force the last integration step to end exactly at T.
            limitStep(model.getT(), init.getT(), step);

            // Integrate the orbit and variational equations.
            ode_integrator::rkf54(model, model.getParams(), step, relTol, absTol);

            ++step.n_int;
            ++step.n_tst;

            // FLI is a running maximum and must therefore be updated
            // after every accepted integration step.
            if (model.getIndicator() == Model::IndicatorType::FLI) {
                indicator_value = chaos_indicator::computeFLI(model.getY(), indicator_value);
            }
        }

        if (valid) {
            // LCI only needs to be evaluated at the final time.
            if (model.getIndicator() == Model::IndicatorType::LCI) {
                indicator_value = chaos_indicator::computeLCI(model.getT(), init.getT0(), model.getY(), init.getDy());
            }
        } else {
            indicator_value = std::numeric_limits<double>::quiet_NaN();
            std::cerr << "\nNon-finite state at grid point "
                      << "(a = " << a << ", e = " << e << "). Proceeding to the next grid point.\n";
        }

        // Write the final indicator value.
        out << std::setw(W) << a << std::setw(W) << e << std::setw(W) << indicator_value << '\n';

        // Display the grid progress.
        grid.printProgress(std::cerr);
    } while (grid.next());
    std::cerr << '\n';
}

void run(const InitData &init, const CommandLineOptions &opt)
{
    constexpr double relTol = 1.0e-6;
    constexpr double absTol = 1.0e-10;

    const double mu    = init.getM2() / (init.getM1() + init.getM2());
    const double mu_12 = astro::sqr(astro::k) * (init.getM1() + init.getM2());
    const double mu_13 = astro::sqr(astro::k) * init.getM1();
    const double n     = std::sqrt(mu_12 / astro::cube(init.getA2()));

    CRTBP2D model(mu, init.getFormalism(), init.getIndicator());
    // Dimensionless integration duration.
    const double TDimless = CRTBP2D::toDimlessTime(init.getT(), n);
    // Dimensionless output time interval.
    const double outputDtDimless = CRTBP2D::toDimlessTime(init.getOutputDt(), n);
    StepControl  step            = createStepControl();

    std::ofstream fout;
    std::ostream *out = openOutputStream(opt, fout);

    switch (init.getRunMode()) {
        case RunMode::ORBIT:
            runOrbit(model, init, *out, step, relTol, absTol, mu_13, n, TDimless, outputDtDimless);
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
        parseCommandLine(argc, argv, opt);

        if (opt.show_version) {
            print::version();
            return EXIT_SUCCESS;
        }
        if (opt.show_help) {
            print::help();
            return EXIT_SUCCESS;
        }

        const InitData init(opt.input_path);

        if (opt.verbose) {
            init.print(std::cout);
        }

        const auto start_time = std::chrono::steady_clock::now();

        run(init, opt);

        const auto                          end_time     = std::chrono::steady_clock::now();
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
