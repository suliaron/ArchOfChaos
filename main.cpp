#include <algorithm>  /**< std::copy, std::copy_n, std::remove_if. */
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
static constexpr const char *PROGRAM_VERSION = "1.0";

static constexpr double PI     = 3.1415926535897932;
static constexpr double TWO_PI = 2.0 * PI;

static constexpr double kG  = 1.720209895e-2;
static constexpr double kG2 = 2.9591220828559110e-4;

constexpr double TIME_EPS = 1.0e-15;

typedef double var_t;
/* prototype of the right-hand side of the diff. eq. system */
typedef void rhs_t(double x, const double *y, double *dy, void *par);

/**
 * @brief Command-line options.
 *
 * Stores the input and output file and directory names together with
 * flags controlling the display of help and version information.
 */
struct CommandLineOptions {
    /** Input file name. */
    std::string input_path;
    /** Input directory. */
    std::string input_dir = ".";
    /** Output file name. */
    std::string output_path;
    /** Output directory. */
    std::string output_dir = ".";
    /** Display help message. */
    bool show_help = false;
    /** Display program version. */
    bool show_version = false;
    /** Display verbose output. */
    bool verbose = false;

    /**
     * @brief Prints the command-line options.
     *
     * Prints all stored command-line options in a human-readable form.
     *
     * @param os Output stream.
     */
    void print(std::ostream& os = std::cout) const
    {
        os << "----------------------------------------\n";
        os << "Command-line options\n";
        os << "----------------------------------------\n";
        os << "Input file    : " << input_path << '\n'
           << "Input dir     : " << input_dir << '\n'
           << "Output file   : " << output_path << '\n'
           << "Output dir    : " << output_dir << '\n'
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
    /// Mass parameter.
    double mu = 0.0;

    /// Initial epoch [time unit].
    double t0 = 0.0;
    /// Length of integration [time unit].
    double T = 0.0;

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

    /// Constructs an empty initialization data structure.
    InitData() = default;

    /// Destroys the initialization data structure.
    ~InitData() = default;

    void print(std::ostream& os = std::cout) const
    {
        os << std::scientific << std::setprecision(6) << std::showpos;

        constexpr int W = 16;

        os << "----------------------------------------\n";
        os << "InitData contents\n";
        os << "----------------------------------------\n";

        os << "mu    = " << std::setw(W) << mu << "\n";

        os << "t0    = " << std::setw(W) << t0 << " [time unit]\n";
        os << "T     = " << std::setw(W) << T << " [time unit]\n";

        os << "a0    = " << std::setw(W) << a0 << " [distance unit]\n";
        os << "a1    = " << std::setw(W) << a1 << " [distance unit]\n";
        os << "Na    = " << std::setw(W) << Na << '\n';

        os << "e0    = " << std::setw(W) << e0 << '\n';
        os << "e1    = " << std::setw(W) << e1 << '\n';
        os << "Ne    = " << std::setw(W) << Ne << '\n';
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
    /**
     * @brief Parameters of the planar circular restricted three-body problem.
     */
    struct CRTBP2DParams
    {
        /**
         * @brief Gravitational mass parameter.
         */
        double mu;

        /**
         * @brief Constructs the parameter set.
         *
         * @param mu Gravitational mass parameter.
         */
        explicit CRTBP2DParams(double mu)
            : mu(mu)
        {}
    };

    /**
     * @brief Computes the right-hand side of the planar circular restricted three-body problem.
     *
     * Evaluates the first-order equations of motion in the rotating (synodic)
     * reference frame. The state vector is defined as
     * y = (x, y, vx, vy), and the output vector contains the corresponding
     * time derivatives
     * dydt = (dx/dt, dy/dt, dvx/dt, dvy/dt).
     *
     * @param t Current time (unused, as the equations are autonomous).
     * @param y Input state vector (x, y, vx, vy).
     * @param dydt Output time derivatives of the state vector.
     * @param par Pointer to a @c CRTBP2DParams structure containing the mass parameter.
     */
    void CRTBP2Dfun(double t, const double *y, double *dydt, void *par)
    {
        const auto  *p  = static_cast<const CRTBP2DParams *>(par);
        const double mu = p->mu;

        const double r1 = std::sqrt(sqr(y[0] + mu) + sqr(y[1]));
        const double r2 = std::sqrt(sqr(y[0] - 1.0 + mu) + sqr(y[1]));

        const double r1_3 = 1.0 / (r1 * r1 * r1);
        const double r2_3 = 1.0 / (r2 * r2 * r2);

        // Equations of motion in the rotating frame.
        dydt[0] = y[2];                                                                                 ///< dx/dt
        dydt[1] = y[3];                                                                                 ///< dy/dt
        dydt[2] = 2.0 * y[3] + y[0] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;  ///< dvx/dt
        dydt[3] = -2.0 * y[2] + y[1] * (1.0 - (1.0 - mu) * r1_3 - mu * r2_3);                           ///< dvy/dt
    }
}  // namespace model

namespace ode_integrator {
    void rkf54(var_t t, StepControl &step, var_t *y_in, var_t *y_out, int n_var, var_t relTol, var_t absTol, rhs_t *fun,
               void *par)
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
        var_t  temax = 0.0;

        fun(t0, y_in, dy.get(), par);
        do {
            temax = 0.0;
            for (int k = 1; k < 7; k++) {
                t = t0 + Ci[k] * step.h;

                for (int n = 0; n < n_var; n++) {
                    y[n] = y_in[n];

                    for (int l = 0; l < k; l++)
                        y[n] += step.h * Aij[k][l] * dy[l * n_var + n];
                }
                fun(t, y.get(), dy.get() + k * n_var, par);
            }

            for (int n = 0; n < n_var; ++n) {
                y_out[n] = y_in[n];

                for (int k = 0; k < 7; ++k) {
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
    }
}  // namespace ode_integrator

namespace {
    /**
     * @brief Prints the program version.
     *
     * Displays the program name together with its version number.
     */
    void print_version()
    {
        std::cout << PROGRAM_NAME << " version " << PROGRAM_VERSION << '\n';
    }

    /**
     * @brief Prints the command-line help.
     *
     * Displays the program usage together with the supported
     * command-line options.
     */
    void print_help()
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
     * @brief Parses the command-line arguments.
     *
     * Processes the supported command-line options and stores the results
     * in the supplied CommandLineOptions structure. The input and output
     * directory names are automatically extracted from the corresponding
     * file paths. If no directory is specified, the current working
     * directory (".") is used.
     *
     * Supported options:
     *   -i <path>    Input file.
     *   -o <path>    Output file.
     *   -h           Display help message.
     *   --help       Display help message.
     *   -v           Display program version.
     *   --version    Display program version.
     *
     * @param[in] argc Number of command-line arguments.
     * @param[in] argv Command-line argument vector.
     * @param[out] opt Parsed command-line options.
     *
     * @throws std::runtime_error
     *         If an unknown option is encountered or if a required
     *         argument is missing.
     */
    void parse_command_line(int argc, char *argv[], CommandLineOptions &opt)
    {
        for (int i = 1; i < argc; ++i) {
            const std::string arg(argv[i]);

            if (arg == "-i") {
                if (++i >= argc)
                    throw std::runtime_error("Missing argument after '-i'.");

                fs::path p(argv[i]);

                opt.input_path = p.filename().string();
                opt.input_dir  = p.has_parent_path() ? p.parent_path().string() : ".";
            } else if (arg == "-o") {
                if (++i >= argc)
                    throw std::runtime_error("Missing argument after '-o'.");

                fs::path p(argv[i]);

                opt.output_path = p.filename().string();
                opt.output_dir  = p.has_parent_path() ? p.parent_path().string() : ".";
            } else if (arg == "-h" || arg == "--help") {
                opt.show_help = true;
            } else if (arg == "-v" || arg == "--version") {
                opt.show_version = true;
            } else if (arg == "--verbose") {
                opt.verbose = true;
            }

            else {
                throw std::runtime_error("Unknown command-line option: " + arg);
            }
        }
    }

    /**
     * @brief Removes all whitespace characters from a string.
     *
     * @param text String to be modified.
     */
    void remove_spaces(std::string &text)
    {
        text.erase(std::remove_if(text.begin(), text.end(), ::isspace), text.end());
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
    void parse_line(const std::string &line, InitData &data)
    {
        std::string text = line;

        const std::size_t comment = text.find('#');
        if (comment != std::string::npos) {
            text.erase(comment);
        }

        remove_spaces(text);

        if (text.empty()) {
            return;
        }

        const std::size_t pos = text.find('=');

        if (pos == std::string::npos) {
            throw std::runtime_error("Missing '=' in initialization file.");
        }

        const std::string key   = text.substr(0, pos);
        const std::string value = text.substr(pos + 1);

        std::istringstream is(value);

        if (key == "mu") {
            is >> data.mu;
        } else if (key == "t0") {
            is >> data.t0;
        } else if (key == "T") {
            is >> data.T;
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
        } else {
            throw std::runtime_error("Unknown keyword: " + key);
        }

        if (!is || !is.eof()) {
            throw std::runtime_error("Invalid value for '" + key + "'.");
        }
    }

    /**
     * @brief Validates the initialization data.
     *
     * @param data Initialization data.
     *
     * @throws std::runtime_error If any parameter is invalid.
     */
    void validate_init_data(const InitData &data)
    {
        if (data.mu <= 0.0 || data.mu >= 0.5) {
            throw std::runtime_error("Invalid mass parameter.");
        }

        if (data.T <= 0.0) {
            throw std::runtime_error("Integration time must be positive.");
        }

        if (data.a1 < data.a0) {
            throw std::runtime_error("a1 must be greater than or equal to a0.");
        }

        if (data.Na == 0) {
            throw std::runtime_error("Na must be positive.");
        }

        if (data.e0 < 0.0 || data.e0 >= 1.0) {
            throw std::runtime_error("e0 must be in the range [0, 1).");
        }

        if (data.e1 < 0.0 || data.e1 >= 1.0) {
            throw std::runtime_error("e1 must be in the range [0, 1).");
        }

        if (data.e1 < data.e0) {
            throw std::runtime_error("e1 must be greater than or equal to e0.");
        }

        if (data.Ne == 0) {
            throw std::runtime_error("Ne must be positive.");
        }
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
    InitData read_init_data(const std::string &file_name)
    {
        std::ifstream file(file_name);

        if (!file) {
            throw std::runtime_error("Cannot open initialization file.");
        }

        InitData data;

        std::string line;

        while (std::getline(file, line)) {
            parse_line(line, data);
        }

        validate_init_data(data);

        return data;
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
    void print_input_data(const CommandLineOptions &options, const InitData &init, std::ostream& os = std::cout)
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
     * @brief Initializes the state vector at the L4 Lagrange point.
     *
     * Places the massless body at the triangular L4 equilibrium point
     * with zero velocity in the rotating reference frame.
     *
     * @param param CRTBP parameters.
     * @param y State vector [x, y, vx, vy].
     */
    void initialize_L4(const model::CRTBP2DParams& param, double* y)
    {
        y[0] = 0.5 - param.mu;
        y[1] = std::sqrt(3.0) / 2.0;
        y[2] = 0.0;
        y[3] = 0.0;
    }

    void getInitialCondition(double mu, double a, double e, double *y)
    {
        y[0] = a * (1.0 - e) - mu;                                      /// x_0
        y[1] = 0.0;                                                     /// y_0
        y[2] = 0.0;                                                     /// vx_0
        y[3] = std::sqrt((1 - mu) / a * (1.0 + e) / (1.0 - e)) - y[0];  /// vy_0
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

    /**
     * @brief Opens the output stream.
     *
     * If an output file is specified on the command line, the file is opened
     * and the returned stream points to it. Otherwise, the standard output
     * stream is returned.
     *
     * @param opt Command-line options.
     * @param fout Output file stream.
     *
     * @return Pointer to the selected output stream.
     *
     * @throws std::runtime_error If the output file cannot be opened.
     */
    std::ostream* open_output_stream(const CommandLineOptions& opt, std::ofstream& fout)
    {
        if (opt.output_path.empty()) {
            return &std::cout;
        }

        const std::filesystem::path output_file =
            std::filesystem::path(opt.output_dir) / opt.output_path;

        fout.open(output_file);

        if (!fout) {
            throw std::runtime_error("Cannot open output file.");
        }

        return &fout;
    }

    /**
     * @brief Prints the current integration state.
     *
     * Prints the current time and the state vector
     * (x, y, vx, vy) in a fixed-width formatted table.
     *
     * @param os Output stream.
     * @param t Current integration time.
     * @param y State vector (x, y, vx, vy).
     */
    void printState(std::ostream& os, double t, const double* y)
    {
        os << std::fixed
            << std::setprecision(6)
            << std::setw(10) << t
            << std::setw(10) << y[0]
            << std::setw(10) << y[1]
            << std::setw(10) << y[2]
            << std::setw(10) << y[3]
            << '\n';
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
}  // namespace

int main(int argc, char *argv[])
{
    CommandLineOptions opt;
    try {
        parse_command_line(argc, argv, opt);
        if (opt.show_version) {
            print_version();
        }
        if (opt.show_help) {
            print_help();
        }

        const InitData init = read_init_data(opt.input_path);

        model::CRTBP2DParams param(init.mu);

        GridIterator grid(init);
        if (opt.verbose) {
            print_input_data(opt, init);
        }

        std::unique_ptr<double[]> y_in =
            allocate_array<double>(4);  // Allocate an array of 4 doubles for initial conditions
        std::unique_ptr<double[]> y_out =
            allocate_array<double>(4);  // Allocate an array of 4 doubles for intermedient results

        /* If the -o option is specified then the output is written to a file */
        /** Output file stream. */
        std::ofstream fout;
        /** Output stream used by the program. */
        std::ostream* out = open_output_stream(opt, fout);

        int    n_var  = 4;  // Number of variables in the system (x, y, vx, vy)
        double relTol = 1.0e-6;
        double absTol = 1.0e-10;
        double t      = init.t0;
        do {
            double a = grid.a();
            double e = grid.e();
            // getInitialCondition(param.mu, a, e, y_in.get());  // Get initial conditions for the current (a,e)
            initialize_L4(param, y_in.get());  // Initialize at L4 point)
            y_in[0] += 3.0e-2;
            // printInitialCondition(grid, y_in.get());          // Print the initial conditions

            t = init.t0;
            StepControl step;  // Step control structure for adaptive integration
            step.h     = 0.1;
            step.h_max = 0.1;
            step.h_nxt = step.h;
            printState(*out, t, y_in.get());  // Print the current state (time and variables)
            do {
                limitStep(t, init.T, step);  // Limit the step size to not exceed the final time
                ode_integrator::rkf54(t, step, y_in.get(), y_out.get(), n_var, relTol, absTol, model::CRTBP2Dfun,
                                      (void *)&param);
                std::copy_n(y_out.get(), n_var, y_in.get());
                t += step.h_did;
                printState(*out, t, y_in.get());  // Print the current state (time and variables)
            } while (std::abs(init.T - t) > 1.0e-10);

            break;
        } while (grid.next());

        return EXIT_SUCCESS;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }

    return EXIT_FAILURE;
}
