#include "crtbp.h"       // CRTBP2D model
#include "grid.h"        // GridIterator
#include "init_data.h"   // InitData, RunMode
#include "math_utils.h"  // astro::sqr, astro::cube
#include "model.h"       // Model base class
#include "orbit.h"       // Orbital-element/state transformations
#include "version.h"     // Program name and version information

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
#include <vector>     // std::vector container.

namespace fs = std::filesystem;

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
     * @brief Converts a floating-point value to a compact string for use
     *        in automatically generated output file names.
     *
     * Decimal points, minus signs, and scientific notation are preserved.
     *
     * @param value Numerical value to convert.
     *
     * @return Compact textual representation of @p value.
     */
    std::string filenameNumber(double value)
    {
        std::ostringstream os;
        os << std::setprecision(8) << std::defaultfloat << value;
        return os.str();
    }

    /**
     * @brief Builds an automatic output file name for ORBIT or INDICATOR mode.
     *
     * The file name contains the computation type, mathematical formalism,
     * integration-duration mode and value, output time interval, P1-P2
     * separation, and all fixed orbital elements.
     *
     * For ORBIT mode, the file name starts with "ORBIT".
     * For INDICATOR mode, it starts with the selected chaos-indicator name.
     *
     * Angular orbital elements are written in degrees.
     *
     * Examples:
     *
     *     ORBIT_NEWTONIAN_nP-100_dt-10_a2-5.2026_
     *     a-5.2_e-0.1_i-0_omega-60_Omega-0_M-30.txt
     *
     *     LCI_NEWTONIAN_nP-100_dt-10_a2-5.2026_
     *     a-5.2_e-0.1_i-0_omega-60_Omega-0_M-30.txt
     *
     * @param init Initialization data.
     *
     * @return Automatically generated output file name.
     *
     * @throws std::runtime_error If the function is called outside ORBIT or
     *         INDICATOR mode, or if the integration-duration or orbital-phase
     *         input mode is unknown.
     */
    std::string buildTimeSeriesOutputFileName(const InitData &init)
    {
        if (init.getRunMode() != RunMode::ORBIT && init.getRunMode() != RunMode::INDICATOR) {
            throw std::runtime_error(
                "buildTimeSeriesOutputFileName() requires "
                "ORBIT or INDICATOR mode.");
        }

        std::ostringstream name;

        // ---------------------------------------------------------------------
        // Computation type.
        // ---------------------------------------------------------------------
        if (init.getRunMode() == RunMode::ORBIT) {
            name << "ORBIT";
        } else {
            name << Model::indicatorTypeToString(init.getIndicator());
        }

        // ---------------------------------------------------------------------
        // Mathematical formalism.
        // ---------------------------------------------------------------------
        name << '_' << Model::formalismToString(init.getFormalism());

        // ---------------------------------------------------------------------
        // Integration duration.
        // ---------------------------------------------------------------------
        if (init.usesPhysicalIntegrationTime()) {
            name << "_T-" << filenameNumber(init.getT());
        } else if (init.usesOrbitalPeriods()) {
            name << "_nP-" << filenameNumber(init.getNPeriods());
        } else {
            throw std::runtime_error("Unknown integration-duration input mode.");
        }

        // ---------------------------------------------------------------------
        // Output time interval.
        // ---------------------------------------------------------------------
        name << "_dt-" << filenameNumber(init.getOutputDt());

        // ---------------------------------------------------------------------
        // P1-P2 separation.
        // ---------------------------------------------------------------------
        name << "_a2-" << filenameNumber(init.getA2());

        // ---------------------------------------------------------------------
        // Fixed orbital elements.
        // ---------------------------------------------------------------------
        const astro::OrbitalElements &elements = init.getElements();
        name << "_a-" << filenameNumber(elements.a);
        name << "_e-" << filenameNumber(elements.e);
        name << "_i-" << filenameNumber(astro::toDeg(elements.i));
        name << "_omega-" << filenameNumber(astro::toDeg(elements.omega));
        name << "_Omega-" << filenameNumber(astro::toDeg(elements.Omega));

        // ---------------------------------------------------------------------
        // Orbital phase.
        // ---------------------------------------------------------------------
        if (init.usesFixedTau()) {
            name << "_tau-" << filenameNumber(elements.tau);
        } else if (init.usesFixedMeanAnomaly()) {
            name << "_M-" << filenameNumber(astro::toDeg(init.getMeanAnomaly()));
        } else {
            throw std::runtime_error("Unknown orbital-phase input mode.");
        }

        // ---------------------------------------------------------------------
        // File extension.
        // ---------------------------------------------------------------------
        name << ".txt";

        return name.str();
    }

    /**
     * @brief Builds an automatic output file name for GRID mode.
     *
     * The file name contains the selected chaos indicator, mathematical
     * formalism, complete grid definitions, integration-duration mode,
     * primary-body separation, and all fixed orbital elements.
     *
     * Each grid axis is represented by its orbital-element name, lower and
     * upper limits, and number of intervals.
     *
     * Example:
     *
     *     LCI_NEWTONIAN_grid-a-4.8to5.8-N250_M-30to330-N30_
     *     nP-100_a2-5.2026_e-0_i-0_omega-60_Omega-0.txt
     *
     * Angular fixed orbital elements are written in degrees.
     *
     * @param init Initialization data.
     *
     * @return Automatically generated output file name.
     *
     * @throws std::runtime_error If the function is called outside GRID mode.
     */
    std::string buildGridOutputFileName(const InitData &init)
    {
        if (init.getRunMode() != RunMode::GRID) {
            throw std::runtime_error("buildGridOutputFileName() requires GRID mode.");
        }

        std::ostringstream name;
        // ---------------------------------------------------------------------
        // Chaos indicator and mathematical formalism.
        // ---------------------------------------------------------------------
        name << Model::indicatorTypeToString(init.getIndicator()) << '_'
             << Model::formalismToString(init.getFormalism());

        // ---------------------------------------------------------------------
        // Grid definitions.
        // ---------------------------------------------------------------------

        name << "_grid";
        for (const GridAxis &axis : init.getGridAxes()) {
            name << '_' << orbitalElementToString(axis.element) << '-' << filenameNumber(axis.min) << "to"
                 << filenameNumber(axis.max) << "-N" << axis.nIntervals;
        }

        // ---------------------------------------------------------------------
        // Integration duration.
        // ---------------------------------------------------------------------
        if (init.usesPhysicalIntegrationTime()) {
            name << "_T-" << filenameNumber(init.getT());

        } else if (init.usesOrbitalPeriods()) {
            name << "_nP-" << filenameNumber(init.getNPeriods());

        } else {
            throw std::runtime_error("Unknown integration-duration input mode.");
        }

        // ---------------------------------------------------------------------
        // P1-P2 separation.
        // ---------------------------------------------------------------------
        name << "_a2-" << filenameNumber(init.getA2());

        // ---------------------------------------------------------------------
        // Fixed orbital elements.
        // ---------------------------------------------------------------------
        const astro::OrbitalElements &elements   = init.getElements();
        const auto                   &gridAxes   = init.getGridAxes();
        const auto                    isGridAxis = [&gridAxes](OrbitalElement element) {
            return std::any_of(gridAxes.begin(), gridAxes.end(),
                               [element](const GridAxis &axis) { return axis.element == element; });
        };

        if (!isGridAxis(OrbitalElement::SEMIMAJOR_AXIS)) {
            name << "_a-" << filenameNumber(elements.a);
        }

        if (!isGridAxis(OrbitalElement::ECCENTRICITY)) {
            name << "_e-" << filenameNumber(elements.e);
        }

        if (!isGridAxis(OrbitalElement::INCLINATION)) {
            name << "_i-" << filenameNumber(astro::toDeg(elements.i));
        }

        if (!isGridAxis(OrbitalElement::ARGUMENT_OF_PERICENTER)) {
            name << "_omega-" << filenameNumber(astro::toDeg(elements.omega));
        }

        if (!isGridAxis(OrbitalElement::LONGITUDE_OF_ASCENDING_NODE)) {
            name << "_Omega-" << filenameNumber(astro::toDeg(elements.Omega));
        }

        // ---------------------------------------------------------------------
        // Fixed orbital phase.
        // ---------------------------------------------------------------------
        if (!isGridAxis(OrbitalElement::PERICENTER_TIME) && !isGridAxis(OrbitalElement::MEAN_ANOMALY)) {
            if (init.usesFixedTau()) {
                name << "_tau-" << filenameNumber(elements.tau);

            } else if (init.usesFixedMeanAnomaly()) {
                name << "_M-" << filenameNumber(astro::toDeg(init.getMeanAnomaly()));
            }
        }

        // ---------------------------------------------------------------------
        // File extension.
        // ---------------------------------------------------------------------
        name << ".txt";

        return name.str();
    }

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
     * @brief Opens the program output stream.
     *
     * If an output file is explicitly specified with the -o command-line
     * option, that file is used. Otherwise, an output file name is generated
     * automatically from the initialization parameters.
     *
     * If verbose output is enabled, the absolute path of the opened output
     * file is printed to the standard output.
     *
     * Automatic output-file naming is currently implemented for GRID mode.
     *
     * @param opt Parsed command-line options.
     * @param init Initialization data.
     * @param fout Output file stream.
     *
     * @return Pointer to the opened output stream.
     *
     * @throws std::runtime_error If an automatic output file name cannot be
     *         generated or the output file cannot be opened.
     */
    std::ostream *openOutputStream(const CommandLineOptions &opt, const InitData &init, std::ofstream &fout)
    {
        fs::path outputPath;

        // Use the explicitly specified output file if -o was given.
        if (!opt.output_path.empty()) {
            outputPath = opt.output_path;
        } else {
            // Generate the output file name automatically.
            switch (init.getRunMode()) {
                case RunMode::ORBIT:
                case RunMode::INDICATOR:
                    outputPath = buildTimeSeriesOutputFileName(init);
                    break;

                case RunMode::GRID:
                    outputPath = buildGridOutputFileName(init);
                    break;

                default:
                    throw std::runtime_error("Unknown run mode while generating output file name.");
            }
        }

        // Convert the output path to an absolute path.
        outputPath = fs::absolute(outputPath);

        // Open the output file.
        fout.open(outputPath);
        if (!fout) {
            throw std::runtime_error("Cannot open output file: " + outputPath.string());
        }

        // Display the actual output path in verbose mode.
        if (opt.verbose) {
            std::cout << "Output file    : " << outputPath.string() << '\n';
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
     * @brief Prints program identification and version information.
     *
     * Displays the program name, version number, author, affiliation,
     * description, and copyright information.
     */
    void version()
    {
        std::cout << program::name << '\n'
                  << "Version     : " << program::version << '\n'
                  << "Author      : " << program::author << '\n'
                  << "Affiliation : " << program::affiliation << '\n'
                  << "Description : " << program::description << '\n'
                  << program::copyright << '\n';
    }

    /**
     * @brief Prints the command-line help.
     *
     * Displays the program name and version, command-line syntax,
     * supported options, and usage examples.
     */
    void help()
    {
        std::cout << program::name << " " << program::version << "\n";

        std::cout << "========================================\n\n";

        std::cout << program::description << "\n\n";

        std::cout << "Usage:\n";
        std::cout << "  " << program::executable << " -i <input file> -o <output file>\n\n";

        std::cout << "Options:\n";
        std::cout << "  -i <file>     Input file.\n";
        std::cout << "  -o <file>     Output file.\n";
        std::cout << "  -h, --help    Display this help message.\n";
        std::cout << "  -v, --version Display program version information.\n";
        std::cout << "  --verbose     Display detailed input information.\n\n";

        std::cout << "Examples:\n";
        std::cout << "  " << program::executable << " -i input.txt -o output.txt\n";

        std::cout << "  " << program::executable << " -i data/input.txt -o results/output.txt\n";
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

    /**
     * @brief Writes the program identification header to an output stream.
     *
     * The header contains the program name and version number and is written
     * as a comment so that numerical data-processing and plotting programs
     * can ignore it.
     *
     * @param os Output stream.
     */
    void outputHeader(std::ostream &os)
    {
        os << "# " << program::name << ' ' << program::version << '\n';
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
 * Constructs the initial heliocentric inertial state, transforms it to the
 * dimensionless rotating CRTBP frame, and integrates the equations of motion.
 * Numerical integration tolerances are obtained from the initialization data.
 *
 * @param model CRTBP model.
 * @param init Initialization data.
 * @param out Output stream.
 * @param step Adaptive integration step-control parameters.
 * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
 *             [AU^3/day^2].
 * @param n Mean motion of the P1-P2 system [rad/day].
 * @param outputDtDimless Dimensionless output time interval.
 *
 * @throws std::runtime_error If a non-finite state is encountered.
 */
void runOrbit(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double mu_13, double n,
              double outputDtDimless)
{
    // Copy the input orbital elements.
    astro::OrbitalElements elements = init.getElements();

    // Calculate the physical integration duration for the current orbit.
    const double duration = init.calcIntegrationDuration(mu_13, elements.a);

    // Convert the physical integration duration to dimensionless CRTBP time.
    const double tDimless = CRTBP2D::toDimlessTime(duration, n);

    // Calculate the pericenter passage time from the selected
    // orbital-phase input.
    elements.tau = init.calc_tau(mu_13, elements.a);

    // Orbital elements -> heliocentric inertial Cartesian state.
    const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

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

        ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

        ++step.n_int;
        ++step.n_tst;

        const bool outputTimeReached = model.getT() >= nextOutputDimless - TIME_EPS;
        const bool finalTimeReached  = model.getT() >= tDimless - TIME_EPS;
        if (outputTimeReached || finalTimeReached) {
            const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);

            if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
                double y_out[4];

                // Hamiltonian -> Newtonian state.
                model.hamiltonianToNewtonian(y_out);

                // Save the Newtonian CRTBP state.
                model.printState(out, physicalTime, y_out);

#if 1  // Test CRTBP -> inertial transformation
                const astro::State inertialState = model.crtbpToInertial(y_out, init.getA2(), n);
                inertialOut << physicalTime << ' ' << inertialState.r.x << ' ' << inertialState.r.y << ' '
                            << inertialState.v.x << ' ' << inertialState.v.y << '\n';
#endif

            } else {
                // Save the Newtonian CRTBP state.
                model.printState(out, physicalTime, model.getY());

#if 1  // Test CRTBP -> inertial transformation
                const astro::State inertialState = model.crtbpToInertial(model.getY(), init.getA2(), n);
                inertialOut << physicalTime << ' ' << inertialState.r.x << ' ' << inertialState.r.y << ' '
                            << inertialState.v.x << ' ' << inertialState.v.y << '\n';
#endif
            }

            if (outputTimeReached) {
                nextOutputDimless += outputDtDimless;
            }
        }
    } /* while */
}

/**
 * @brief Integrates a single CRTBP orbit and computes a chaos indicator.
 *
 * Initializes the orbital state and deviation vector, integrates the
 * equations of motion and variational equations in dimensionless CRTBP time,
 * and writes the selected chaos indicator at the requested physical output
 * times.
 *
 * Physical input times are specified in days, while the numerical integration
 * is performed using dimensionless CRTBP time.
 *
 * @param model CRTBP model.
 * @param init Initialization data.
 * @param out Output stream.
 * @param step Adaptive integration step-control parameters.
 * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
 *             [AU^3/day^2].
 * @param n Mean motion of the P1-P2 system [rad/day].
 * @param outputDtDimless Dimensionless output time interval.
 *
 * @throws std::runtime_error If a non-finite state is encountered.
 */
void runIndicator(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double mu_13, double n,
                  double outputDtDimless)
{
    constexpr int W = 18;

    // Copy the input orbital elements.
    astro::OrbitalElements elements = init.getElements();

    // Calculate the physical integration duration for the current orbit.
    const double duration = init.calcIntegrationDuration(mu_13, elements.a);

    // Convert the physical integration duration to dimensionless CRTBP time.
    const double tDimless = CRTBP2D::toDimlessTime(duration, n);

    // Calculate the pericenter passage time from the selected
    // orbital-phase input.
    elements.tau = init.calc_tau(mu_13, elements.a);

    // Orbital elements -> heliocentric inertial Cartesian state.
    const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

    // Heliocentric inertial state -> dimensionless rotating CRTBP state.
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
            indicator_value = chaos_indicator::computeFLI(model.getY(), indicator_value);
            out << std::left << std::setw(W) << "t [day]" << std::setw(W) << "FLI" << '\n';
            out << std::right << std::scientific << std::setprecision(10);

            // FLI is defined at the initial physical time.
            out << std::setw(W) << init.getT0() << std::setw(W) << indicator_value << '\n';
            break;

        case Model::IndicatorType::LCI:
            out << std::left << std::setw(W) << "t [day]" << std::setw(W) << "LCI [1/day]" << '\n';
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

    // First output epoch in dimensionless CRTBP time.
    double nextOutputDimless = outputDtDimless;

    while (model.getT() < tDimless - TIME_EPS) {
        // Check the numerical state periodically.
        if (step.n_tst % 10 == 0) {
            if (!checkFinite(model.getY(), model.getNVar())) {
                throw std::runtime_error("Non-finite state encountered during indicator integration.");
            }
        }

        // Stop exactly at the next output time or at the final time.
        const double targetTime = std::min(nextOutputDimless, tDimless);

        limitStep(model.getT(), targetTime, step);

        // Integrate the orbit and the variational equations.
        ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

        // Physical time corresponding to the current
        // dimensionless CRTBP time.
        const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);

        // Update indicators that depend on all intermediate integration steps.
        switch (model.getIndicator()) {
            case Model::IndicatorType::FLI:
                indicator_value = chaos_indicator::computeFLI(model.getY(), indicator_value);
                break;

            case Model::IndicatorType::LCI:
                indicator_value = chaos_indicator::computeLCI(model.getT(), 0.0, model.getY(), init.getDy());
                break;

            default:
                break;
        }

        ++step.n_int;
        ++step.n_tst;

        const bool outputTimeReached = model.getT() >= nextOutputDimless - TIME_EPS;
        const bool finalTimeReached  = model.getT() >= tDimless - TIME_EPS;
        if (outputTimeReached || finalTimeReached) {
            out << std::setw(W) << physicalTime << std::setw(W) << indicator_value << '\n';

            if (outputTimeReached) {
                nextOutputDimless += outputDtDimless;
            }
        }
    }
}

/**
 * @brief Computes a chaos indicator over a two-dimensional orbital grid.
 *
 * Integrates the CRTBP equations and variational equations for every
 * semimajor-axis and eccentricity pair of the grid and writes the final
 * value of the selected chaos indicator.
 *
 * Physical input times are specified in days, while the numerical
 * integration is performed using dimensionless CRTBP time.
 *
 * @param model CRTBP model.
 * @param init Initialization data.
 * @param out Output stream.
 * @param step Adaptive integration step-control parameters.
 * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
 *              [AU^3/day^2].
 * @param n Mean motion of the P1-P2 system [rad/day].
 */
void runGrid(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double mu_13, double n)
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
        // Reset the dimensionless CRTBP time.
        model.setT(0.0);

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
        // Calculate the physical integration duration for the current
        // semimajor axis.
        const double duration = init.calcIntegrationDuration(mu_13, a);

        // Convert the physical integration duration to dimensionless CRTBP time.
        const double tDimless = CRTBP2D::toDimlessTime(duration, n);

        // Construct the orbital elements for the current grid point.
        astro::OrbitalElements elements = init.getElements();

        elements.a = a;
        elements.e = e;

        // Calculate the pericenter passage time for the current
        // semimajor axis.
        elements.tau = init.calc_tau(mu_13, elements.a);

        // Orbital elements -> heliocentric inertial Cartesian state.
        const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

        // Heliocentric inertial state -> dimensionless rotating CRTBP state.
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

        // Integrate the current grid point in dimensionless CRTBP time.
        while (model.getT() < tDimless - TIME_EPS) {
            // Check the numerical state periodically.
            if (step.n_tst % 10 == 0) {
                if (!checkFinite(model.getY(), model.getNVar())) {
                    valid = false;
                    break;
                }
            }

            // Force the last integration step to end exactly at
            // the final dimensionless time.
            limitStep(model.getT(), tDimless, step);

            // Integrate the orbit and variational equations.
            ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

            ++step.n_int;
            ++step.n_tst;

            // FLI is a running maximum and must therefore be updated
            // after every accepted integration step.
            if (model.getIndicator() == Model::IndicatorType::FLI) {
                indicator_value = chaos_indicator::computeFLI(model.getY(), indicator_value);
            }
        } /* while (model.getT() < tDimless - TIME_EPS) */

        if (valid) {
            // LCI only needs to be evaluated at the final time.
            if (model.getIndicator() == Model::IndicatorType::LCI) {
                const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);
                indicator_value = chaos_indicator::computeLCI(physicalTime, init.getT0(), model.getY(), init.getDy());
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

/**
 * @brief Computes a chaos indicator over a multidimensional
 *        orbital-element grid.
 *
 * Integrates the CRTBP equations and variational equations for every
 * point of an arbitrary orbital-element grid and writes the final value
 * of the selected chaos indicator.
 *
 * Fixed orbital elements are obtained from InitData, while grid-controlled
 * elements are replaced by the current values supplied by GridIterator.
 *
 * The orbital phase may be specified either by tau or M, as a fixed input
 * value or as a grid axis. If M is a grid axis, the corresponding time of
 * pericenter passage is calculated from the current semimajor axis.
 *
 * If the integration duration is specified by nPeriods, the physical
 * duration is recalculated independently at every grid point from the
 * current semimajor axis of P3.
 *
 * Physical times are expressed in days, while numerical integration is
 * performed in dimensionless CRTBP time.
 *
 * @param model CRTBP model.
 * @param init Initialization data.
 * @param out Output stream.
 * @param step Adaptive integration step-control parameters.
 * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
 *              [AU^3/day^2].
 * @param n Mean motion of the P1-P2 system [rad/day].
 *
 * @throws std::runtime_error If the selected chaos indicator is invalid.
 */
void runGridGeneral(CRTBP2D &model, const InitData &init, std::ostream &out, StepControl &step, double mu_13, double n)
{
    constexpr int W = 18;

    // Construct the multidimensional orbital-element grid iterator.
    GridIterator grid(init.getGridAxes());

    // ---------------------------------------------------------------------
    // Check the selected chaos indicator.
    // ---------------------------------------------------------------------

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

    // ---------------------------------------------------------------------
    // Write the output header.
    // ---------------------------------------------------------------------

    // Write the grid-axis column headers.
    for (const GridAxis &axis : init.getGridAxes()) {
        const std::string label =
            std::string(orbitalElementToString(axis.element)) + " [" + orbitalElementUnit(axis.element) + "]";

        out << std::left << std::setw(W) << label;
    }

    // Write the chaos-indicator column header.
    switch (model.getIndicator()) {
        case Model::IndicatorType::FLI:
            out << std::setw(W) << "FLI";
            break;

        case Model::IndicatorType::LCI:
            out << std::setw(W) << "LCI [1/day]";
            break;

        default:
            break;
    }
    out << '\n';
    out << std::right << std::scientific << std::setprecision(10);

    // ---------------------------------------------------------------------
    // Save the initial adaptive step-control values.
    // ---------------------------------------------------------------------

    const double initial_h     = step.h;
    const double initial_h_max = step.h_max;
    const double initial_h_min = step.h_min;

    // ---------------------------------------------------------------------
    // Iterate over all grid points.
    // ---------------------------------------------------------------------
    do {
        // Reset the dimensionless CRTBP integration time.
        model.setT(0.0);

        // Reset the adaptive step-size control.
        step.h     = initial_h;
        step.h_nxt = initial_h;
        step.h_did = 0.0;
        step.h_max = initial_h_max;
        step.h_min = initial_h_min;
        step.n_tst = 0;
        step.n_int = 0;

        // -----------------------------------------------------------------
        // Construct the orbital elements for the current grid point.
        // -----------------------------------------------------------------

        // Start from the fixed orbital elements specified in the input file.
        astro::OrbitalElements elements = init.getElements();

        // Apply all grid-controlled orbital elements that are stored
        // directly in astro::OrbitalElements.
        grid.apply(elements);

        // -----------------------------------------------------------------
        // Resolve the orbital phase.
        // -----------------------------------------------------------------

        if (grid.hasAxis(OrbitalElement::MEAN_ANOMALY)) {
            // Mean anomaly is specified in degrees in the grid.
            const double meanAnomaly = astro::toRad(grid.getValue(OrbitalElement::MEAN_ANOMALY));

            // Keplerian mean motion of P3 [rad/day].
            const double n_3 = std::sqrt(mu_13 / astro::cube(elements.a));

            // M(t0) = n_3 * (t0 - tau)
            //
            // therefore
            //
            // tau = t0 - M(t0) / n_3.
            elements.tau = init.getT0() - meanAnomaly / n_3;

        } else if (grid.hasAxis(OrbitalElement::PERICENTER_TIME)) {
            // tau has already been assigned by grid.apply().
            // No additional conversion is required.

        } else {
            // The orbital phase is fixed in the input file.
            //
            // calc_tau() returns the specified tau directly or calculates
            // it from the fixed mean anomaly using the current semimajor axis.
            elements.tau = init.calc_tau(mu_13, elements.a);
        }

        // -----------------------------------------------------------------
        // Calculate the integration duration for the current grid point.
        // -----------------------------------------------------------------

        const double duration = init.calcIntegrationDuration(mu_13, elements.a);
        const double tDimless = CRTBP2D::toDimlessTime(duration, n);

        // -----------------------------------------------------------------
        // Construct the initial CRTBP state.
        // -----------------------------------------------------------------

        // Orbital elements -> heliocentric inertial Cartesian state.
        const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

        // Heliocentric inertial state ->
        // dimensionless rotating CRTBP state.
        model.inertialToCRTBP(state, init.getA2(), n);

        // Convert only the orbital state to Hamiltonian canonical
        // variables if required.
        if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
            model.velocityToHamiltonian();
        }

        // Set the initial deviation vector exactly as specified
        // in the input file.
        std::copy_n(init.getDy(), 4, model.getY() + 4);

        // -----------------------------------------------------------------
        // Initialize the selected chaos indicator.
        // -----------------------------------------------------------------

        double indicator_value = 0.0;

        if (model.getIndicator() == Model::IndicatorType::FLI) {
            indicator_value = chaos_indicator::computeFLI(model.getY(), 1.0);
        }

        bool valid = true;

        // -----------------------------------------------------------------
        // Integrate the current grid point.
        // -----------------------------------------------------------------

        while (model.getT() < tDimless - TIME_EPS) {
            // Check the numerical state periodically.
            if (step.n_tst % 10 == 0) {
                if (!checkFinite(model.getY(), model.getNVar())) {
                    valid = false;
                    break;
                }
            }

            // Force the last integration step to end exactly at
            // the final dimensionless integration time.
            limitStep(model.getT(), tDimless, step);

            // Integrate the orbit and variational equations.
            ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

            ++step.n_int;
            ++step.n_tst;

            // FLI is a running maximum and must therefore be updated
            // after every accepted integration step.
            if (model.getIndicator() == Model::IndicatorType::FLI) {
                indicator_value = chaos_indicator::computeFLI(model.getY(), indicator_value);
            }
        } /* while */

        // -----------------------------------------------------------------
        // Final indicator value.
        // -----------------------------------------------------------------

        if (valid) {
            // LCI only needs to be evaluated at the final integration time.
            if (model.getIndicator() == Model::IndicatorType::LCI) {
                // Convert the final dimensionless CRTBP time to physical time [day].
                const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);
                indicator_value = chaos_indicator::computeLCI(physicalTime, init.getT0(), model.getY(), init.getDy());
            }

        } else {
            indicator_value = std::numeric_limits<double>::quiet_NaN();
            std::cerr << "\nNon-finite state encountered at grid point:\n";
            grid.printCurrentPoint(std::cerr);
            std::cerr << "Proceeding to the next grid point.\n";
        }

        // -----------------------------------------------------------------
        // Write the current grid coordinates and the indicator value.
        // -----------------------------------------------------------------

        for (std::size_t axisIndex = 0; axisIndex < grid.getAxisCount(); ++axisIndex) {
            out << std::setw(W) << grid.getValue(axisIndex);
        }

        out << std::setw(W) << indicator_value << '\n';

        // Display the grid progress.
        grid.printProgress(std::cerr);
    } while (grid.next());

    std::cerr << '\n';
}

void run(const InitData &init, const CommandLineOptions &opt)
{
    const double mu    = init.getM2() / (init.getM1() + init.getM2());
    const double mu_12 = astro::sqr(astro::k) * (init.getM1() + init.getM2());
    const double mu_13 = astro::sqr(astro::k) * init.getM1();
    const double n     = std::sqrt(mu_12 / astro::cube(init.getA2()));

    CRTBP2D model(mu, init.getFormalism(), init.getIndicator());

    // Dimensionless output time interval.
    const double outputDtDimless = CRTBP2D::toDimlessTime(init.getOutputDt(), n);
    StepControl  step            = createStepControl();

    std::ofstream fout;
    std::ostream *out = openOutputStream(opt, init, fout);
    // Write program identification to the output.
    print::outputHeader(*out);

    switch (init.getRunMode()) {
        case RunMode::ORBIT:
            runOrbit(model, init, *out, step, mu_13, n, outputDtDimless);
            break;

        case RunMode::INDICATOR:
            runIndicator(model, init, *out, step, mu_13, n, outputDtDimless);
            break;

        case RunMode::GRID:
            // runGrid(     model, init, *out, step, mu_13, n);
            runGridGeneral(model, init, *out, step, mu_13, n);
            break;

        default:
            throw std::runtime_error("Unknown run mode.");
    }
}

int main(int argc, char *argv[])
{
    try {
        // Test GridIterator traversal on a simple two-dimensional (a,e) grid.
        // The test verifies the grid-point order, axis indexing, current values,
        // and the total number of visited grid points.
#if 0
        {
            const std::vector<GridAxis> axes = {
                {
                    OrbitalElement::ECCENTRICITY,
                    0.0,
                    0.2,
                    1
                },
                {
                    OrbitalElement::SEMIMAJOR_AXIS,
                    3.0,
                    5.0,
                    2
                }
            };

            testGridTraversal(axes);

            return 0;
        }
#endif

        // Test GridIterator traversal with a full-period mean-anomaly axis.
        // The test verifies that the periodic M grid excludes the duplicated
        // 360-degree endpoint and visits exactly 2 x 4 = 8 grid points.
#if 0
        const std::vector<GridAxis> axes = {
            {
                OrbitalElement::SEMIMAJOR_AXIS,
                3.0,
                4.0,
                1
            },
            {
                OrbitalElement::MEAN_ANOMALY,
                0.0,
                360.0,
                4
            }
        };

        testGridTraversal(axes);

        return 0;
#endif

        // Test GridIterator traversal on a three-dimensional (a,e,M) grid.
        // The test verifies multi-level index carry and periodic endpoint handling.
#if 0
        const std::vector<GridAxis> axes = {
            {
                OrbitalElement::SEMIMAJOR_AXIS,
                3.0,
                4.0,
                1
            },
            {
                OrbitalElement::ECCENTRICITY,
                0.0,
                0.2,
                2
            },
            {
                OrbitalElement::MEAN_ANOMALY,
                0.0,
                360.0,
                4
            }
        };
        testGridTraversal(axes);
        return 0;
#endif

#if 0
        // Test GridIterator traversal on a partial mean-anomaly interval.
        // The test verifies that both endpoints are included for a non-full-period
        // angular grid and that M advances from 0 to 180 degrees in 10-degree steps.
        const std::vector<GridAxis> axes = {{OrbitalElement::MEAN_ANOMALY, 0.0, 180.0, 18}};
        testGridTraversal(axes);
        return 0;
#endif

        // Test applying current GridIterator values to orbital elements.
        // The test verifies that a, e, and omega are updated according
        // to the current grid point and that angular values are converted
        // from degrees to radians.
#if 0
        {
            const std::vector<GridAxis> axes = {{OrbitalElement::SEMIMAJOR_AXIS, 3.0, 5.0, 2},
                                                {OrbitalElement::ECCENTRICITY, 0.0, 0.2, 2},
                                                {OrbitalElement::ARGUMENT_OF_PERICENTER, 0.0, 180.0, 2}};

            GridIterator grid(axes);

            std::cout << "----------------------------------------\n";
            std::cout << "GridIterator apply() test\n";
            std::cout << "----------------------------------------\n";

            std::size_t visited = 0;

            do {
                astro::OrbitalElements elements{};

                // Set fixed orbital elements to easily recognizable values.
                elements.i     = astro::toRad(5.0);
                elements.Omega = astro::toRad(15.0);
                elements.tau   = 123.0;

                // Apply the current grid-controlled orbital elements.
                grid.apply(elements);

                std::cout << "point " << visited << " : "
                          << "a=" << elements.a << "  e=" << elements.e << "  omega=" << elements.omega << " rad"
                          << "  i=" << elements.i << " rad"
                          << "  Omega=" << elements.Omega << " rad"
                          << "  tau=" << elements.tau << '\n';

                ++visited;

            } while (grid.next());

            std::cout << '\n';
            std::cout << "visited     : " << visited << '\n';
            std::cout << "expected    : " << grid.getTotalPointCount() << '\n';

            if (visited != grid.getTotalPointCount()) {
                throw std::runtime_error(
                    "GridIterator apply() test failed: "
                    "incorrect number of visited grid points.");
            }

            std::cout << "result      : OK\n";

            return 0;
        }
#endif

        CommandLineOptions opt;
        parseCommandLine(argc, argv, opt);

#if 0
        {
            const InitData init(opt.input_path);

            std::cout << "----------------------------------------\n";
            std::cout << "Automatic GRID output filename test\n";
            std::cout << "----------------------------------------\n";

            std::cout << buildGridOutputFileName(init) << '\n';

            return EXIT_SUCCESS;
        }
#endif

#if 0
        {
            const InitData init(opt.input_path);

            std::cout << "----------------------------------------\n";
            std::cout << "Automatic ORBIT output filename test\n";
            std::cout << "----------------------------------------\n";

            std::cout << buildTimeSeriesOutputFileName(init) << '\n';

            return EXIT_SUCCESS;
        }
#endif

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
