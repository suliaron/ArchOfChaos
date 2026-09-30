#include "io.h"            // Input/output helpers and reproducibility metadata
#include "command_line.h"  // CommandLineOptions
#include "grid.h"          // GridAxis, OrbitalElement, and grid-element helpers
#include "init_data.h"     // InitData, RunMode, ProblemType
#include "math_utils.h"    // astro::toDeg
#include "model.h"         // Model and chaos-indicator/formalism helpers
#include "orbit.h"         // Orbital-element definitions
#include "version.h"       // Program identification and version information

#include <algorithm>  /**<  */
#include <chrono>     /**< std::chrono::system_clock */
#include <cstddef>    /**< std::size_t */
#include <ctime>      /**< std::time_t, std::tm */
#include <filesystem> /**< */
#include <fstream>    /**< std::ifstream */
#include <iomanip>    /**< std::put_time, std::setprecision, std::setw */
#include <iostream>   /**< */
#include <ostream>    /**< std::ostream */
#include <sstream>    /**< std::ostringstream */
#include <stdexcept>  /**< std::runtime_error */
#include <string>     /**< std::string */
#include <vector>     /**< std::vector */

#if defined(_WIN32)
#include <windows.h>  // GetComputerNameA
#else
#include <unistd.h>  // gethostname
#endif

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
     * @brief Returns a human-readable description of an orbital element.
     *
     * @param element Orbital element.
     *
     * @return Human-readable description.
     *
     * @throws std::runtime_error If the orbital element is unknown.
     */
    const char *orbitalElementDescription(OrbitalElement element)
    {
        switch (element) {
            case OrbitalElement::SEMIMAJOR_AXIS:
                return "Semimajor axis";

            case OrbitalElement::ECCENTRICITY:
                return "Eccentricity";

            case OrbitalElement::INCLINATION:
                return "Inclination";

            case OrbitalElement::ARGUMENT_OF_PERICENTER:
                return "Argument of pericenter";

            case OrbitalElement::LONGITUDE_OF_ASCENDING_NODE:
                return "Longitude of ascending node";

            case OrbitalElement::PERICENTER_TIME:
                return "Time of pericenter passage";

            case OrbitalElement::MEAN_ANOMALY:
                return "Mean anomaly";

            default:
                throw std::runtime_error("Unknown orbital element while building output schema.");
        }
    }

    /**
     * @brief Returns the current local date and time.
     *
     * The timestamp is formatted as
     *
     *     YYYY-MM-DD HH:MM:SS
     *
     * The implementation supports both Microsoft Visual C++ and GCC.
     *
     * @return Current local timestamp.
     *
     * @throws std::runtime_error If the local time cannot be determined.
     */
    std::string currentTimestamp()
    {
        const auto now = std::chrono::system_clock::now();

        const std::time_t time = std::chrono::system_clock::to_time_t(now);

        std::tm localTime{};

#if defined(_MSC_VER)

        if (localtime_s(&localTime, &time) != 0) {
            throw std::runtime_error("Cannot convert the current time to local time.");
        }

#elif defined(__unix__) || defined(__APPLE__)

        if (localtime_r(&time, &localTime) == nullptr) {
            throw std::runtime_error("Cannot convert the current time to local time.");
        }

#else

        const std::tm *tmp = std::localtime(&time);

        if (tmp == nullptr) {
            throw std::runtime_error("Cannot convert the current time to local time.");
        }

        localTime = *tmp;

#endif

        std::ostringstream os;

        os << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S");

        return os.str();
    }

    /**
     * @brief Returns the current build configuration.
     *
     * The result is inferred from the standard NDEBUG macro.
     *
     * @return "Release" if NDEBUG is defined, otherwise "Debug".
     */
    const char *buildConfiguration()
    {
#ifdef NDEBUG
        return "Release";
#else
        return "Debug";
#endif
    }

    /**
     * @brief Returns a description of the compiler used to build the program.
     *
     * @return Compiler name and version.
     */
    std::string compilerDescription()
    {
        std::ostringstream os;

#if defined(_MSC_VER)
        os << "Microsoft Visual C++ " << _MSC_VER;
#elif defined(__clang__)
        os << "Clang " << __clang_major__ << '.' << __clang_minor__ << '.' << __clang_patchlevel__;
#elif defined(__GNUC__)
        os << "GCC " << __GNUC__ << '.' << __GNUC_MINOR__ << '.' << __GNUC_PATCHLEVEL__;
#else
        os << "Unknown compiler";
#endif
        return os.str();
    }

    /**
     * @brief Returns the C++ language standard used for compilation.
     *
     * @return C++ standard name.
     */
    const char *cppStandardDescription()
    {
#if defined(_MSVC_LANG)
        constexpr long cppVersion = _MSVC_LANG;
#else
        constexpr long cppVersion = __cplusplus;
#endif

        if (cppVersion >= 202302L) {
            return "C++23";
        }
        if (cppVersion >= 202002L) {
            return "C++20";
        }
        if (cppVersion >= 201703L) {
            return "C++17";
        }
        if (cppVersion >= 201402L) {
            return "C++14";
        }
        if (cppVersion >= 201103L) {
            return "C++11";
        }
        return "pre-C++11";
    }

    /**
     * @brief Returns a short description of the target platform.
     *
     * @return Operating-system and architecture description.
     */
    const char *platformDescription()
    {
#if defined(_WIN32)

#if defined(_WIN64) || defined(__x86_64__) || defined(__amd64__)
        return "Windows x86-64";
#elif defined(_M_IX86) || defined(__i386__)
        return "Windows x86";
#elif defined(_M_ARM64) || defined(__aarch64__)
        return "Windows ARM64";
#else
        return "Windows";
#endif

#elif defined(__linux__)

#if defined(__x86_64__) || defined(__amd64__)
        return "Linux x86-64";
#elif defined(__aarch64__)
        return "Linux ARM64";
#elif defined(__i386__)
        return "Linux x86";
#else
        return "Linux";
#endif

#elif defined(__APPLE__)

#if defined(__aarch64__)
        return "macOS ARM64";
#elif defined(__x86_64__)
        return "macOS x86-64";
#else
        return "macOS";
#endif

#else
        return "Unknown platform";
#endif
    }

    /**
     * @brief Returns the host name of the computer running the program.
     *
     * Uses the native Windows API on Windows systems and gethostname()
     * on POSIX-compatible systems.
     *
     * @return Host name of the current computer, or "unknown" if it
     *         cannot be determined.
     */
    std::string hostName()
    {
        char name[256] = {};

#if defined(_WIN32)

        DWORD size = static_cast<DWORD>(sizeof(name));

        if (GetComputerNameA(name, &size) != 0) {
            return std::string(name, size);
        }

#else

        if (gethostname(name, sizeof(name)) == 0) {
            name[sizeof(name) - 1] = '\0';
            return std::string(name);
        }

#endif

        return "unknown";
    }
}  // namespace

namespace io {

    /**
     * @brief Builds an automatic output file name for ORBIT or INDICATOR mode.
     *
     * The file name contains the computation type, mathematical formalism,
     * integration-duration mode and value, output schedule, P1-P2
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
     *     LCI_NEWTONIAN_nP-100_log-0.0001-N9_a2-5.2026_
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
        }
        else {
            name << Model::indicatorTypeToString(init.getIndicator());
        }

        // ---------------------------------------------------------------------
        // Problem type.
        // ---------------------------------------------------------------------
        switch (init.getProblemType()) {
            case ProblemType::CRTBP2D:
                name << "_CRTBP2D";
                break;

            case ProblemType::CRTBP3D:
                name << "_CRTBP3D";
                break;

            default:
                throw std::runtime_error("Unknown problem type while building output file name.");
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
        }
        else if (init.usesOrbitalPeriods()) {
            name << "_nP-" << filenameNumber(init.getNPeriods());
        }
        else {
            throw std::runtime_error("Unknown integration-duration input mode.");
        }

        // ---------------------------------------------------------------------
        // Output schedule.
        // ---------------------------------------------------------------------
        if (init.getRunMode() == RunMode::ORBIT) {
            // Linear output interval for orbit integration.
            name << "_dt-" << filenameNumber(init.getOutputDt());
        }
        else {
            // Logarithmic output schedule for chaos indicators.
            name << "_log-" << filenameNumber(init.getOutputFirst()) << "-N" << init.getOutputPointsPerDecade();
        }

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
        }
        else if (init.usesFixedMeanAnomaly()) {
            name << "_M-" << filenameNumber(astro::toDeg(init.getMeanAnomaly()));
        }
        else {
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
        }
        else if (init.usesOrbitalPeriods()) {
            name << "_nP-" << filenameNumber(init.getNPeriods());
        }
        else {
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
            }
            else if (init.usesFixedMeanAnomaly()) {
                name << "_M-" << filenameNumber(astro::toDeg(init.getMeanAnomaly()));
            }
        }

        // ---------------------------------------------------------------------
        // File extension.
        // ---------------------------------------------------------------------
        name << ".txt";

        return name.str();
    }

    std::ostream *openOutputStream(const CommandLineOptions &opt, const InitData &init, std::ofstream &fout)
    {
        std::filesystem::path outputPath;

        // Use the explicitly specified output file if -o was given.
        if (!opt.output_path.empty()) {
            outputPath = opt.output_path;
        }
        else {
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

            // If -oDir was specified, place the automatically generated
            // output file in that directory.
            if (!opt.output_dir.empty()) {
                outputPath = std::filesystem::path(opt.output_dir) / outputPath;
            }
        }

        outputPath = std::filesystem::absolute(outputPath);

        fout.open(outputPath);

        if (!fout) {
            throw std::runtime_error("Cannot open output file: " + outputPath.string());
        }

        if (opt.verbose) {
            std::cout << "Output file    : " << outputPath.string() << '\n';
        }

        return &fout;
    }

    void configureNumericalOutput(std::ostream &out)
    {
        out << std::right << std::scientific << std::showpos << std::setprecision(DATA_PRECISION);
    }

    void writeOutputHeader(std::ostream &out, const std::filesystem::path &inputPath, const InitData &init)
    {
        writeProgramInformation(out, inputPath, init);
        out << "#\n";
        writeInputFileCopy(out, inputPath);
        out << "#\n";
        writeOutputStructure(out, init);
        out << "#\n";
        out << "# ==============================================================================\n";
        out << "# DATA\n";
        out << "# ==============================================================================\n";
    }

    void writeInputFileCopy(std::ostream &out, const std::filesystem::path &inputPath)
    {
        std::ifstream input(inputPath);

        if (!input) {
            throw std::runtime_error("Cannot open input file for header copy: " + inputPath.string());
        }

        out << "# ==============================================================================\n";
        out << "# BEGIN INPUT FILE\n";
        out << "# ==============================================================================\n";

        std::string line;

        while (std::getline(input, line)) {
            out << "# | " << line << '\n';
        }

        out << "# ==============================================================================\n";
        out << "# END INPUT FILE\n";
        out << "# ==============================================================================\n";
    }

    void writeProgramInformation(std::ostream &out, const std::filesystem::path &inputPath, const InitData &init)
    {
        out << "# ==============================================================================\n";
        out << "# " << program::name << '\n';
        out << "# ==============================================================================\n";
        out << "#\n";

        out << "# Program information\n";
        out << "# -------------------\n";
        out << "# program          : " << program::name << '\n';
        out << "# version          : " << program::version << '\n';
        out << "# author           : " << program::author << '\n';
        out << "# affiliation      : " << program::affiliation << '\n';
        out << "# description      : " << program::description << '\n';
        out << "# copyright        : " << program::copyright << '\n';
        out << "# run timestamp    : " << currentTimestamp() << '\n';
        out << "# build timestamp  : " << __DATE__ << ' ' << __TIME__ << '\n';
        out << "# build type       : " << buildConfiguration() << '\n';
        out << "# compiler         : " << compilerDescription() << '\n';
        out << "# C++ standard     : " << cppStandardDescription() << '\n';
        out << "# platform         : " << platformDescription() << '\n';
        out << "# host name        : " << hostName() << '\n';
        out << "# numerical method : RKF54 - adaptive Runge-Kutta-Fehlberg 5(4)\n";
        out << "# run mode         : " << runModeToString(init.getRunMode()) << '\n';
        out << "# indicator        : " << Model::indicatorTypeToString(init.getIndicator()) << '\n';
        out << "# formalism        : " << Model::formalismToString(init.getFormalism()) << '\n';
        // -------------------------------------------------------------------------
        // Output schedule.
        // -------------------------------------------------------------------------
        switch (init.getRunMode()) {
            case RunMode::ORBIT:
                out << "# output schedule  : linear\n";
                out << "# output interval  : " << init.getOutputDt() << " day\n";
                break;

            case RunMode::INDICATOR:
                out << "# output schedule  : logarithmic\n";
                out << "# first output     : " << init.getOutputFirst() << " day\n";
                out << "# points/decade    : " << init.getOutputPointsPerDecade() << '\n';
                break;

            case RunMode::GRID:
                out << "# output schedule  : final value only\n";
                break;

            default:
                break;
        }
        out << "# input file       : " << inputPath.filename().string() << '\n';
    }

    void writeOutputStructure(std::ostream &out, const InitData &init)
    {
        const std::vector<OutputColumn> columns = buildOutputSchema(init);

        out << "# Output data structure\n";
        out << "# ---------------------\n";
        out << "# number of columns : " << columns.size() << '\n';

        for (std::size_t i = 0; i < columns.size(); ++i) {
            const OutputColumn &column = columns[i];

            out << "# column " << (i + 1) << "          : " << column.name << " [" << column.unit << "]"
                << " - " << column.description << '\n';
        }
    }

    void writeTableHeader(std::ostream &out, const InitData &init)
    {
        const std::vector<OutputColumn> columns = buildOutputSchema(init);

        for (const OutputColumn &column : columns) {
            std::string label = column.name;

            if (!column.unit.empty()) {
                label += " [" + column.unit + "]";
            }
            out << std::left << std::setw(DATA_FIELD_WIDTH) << label << ' ';
        }
        out << '\n';

        // Restore right alignment for subsequent numerical output.
        out << std::right;
    }

    std::vector<OutputColumn> buildOutputSchema(const InitData &init)
    {
        std::vector<OutputColumn> columns;

        switch (init.getRunMode()) {
                // -----------------------------------------------------------------
                // ORBIT
                // -----------------------------------------------------------------
            case RunMode::ORBIT:
                columns.push_back({"t", "day", "Physical time"});
                columns.push_back({"x", "-", "Dimensionless rotating-frame x coordinate"});
                columns.push_back({"y", "-", "Dimensionless rotating-frame y coordinate"});
                if (init.getProblemType() == ProblemType::CRTBP3D) {
                    columns.push_back({"z", "-", "Dimensionless rotating-frame z coordinate"});
                }
                columns.push_back({"vx", "-", "Dimensionless rotating-frame x velocity"});
                columns.push_back({"vy", "-", "Dimensionless rotating-frame y velocity"});
                if (init.getProblemType() == ProblemType::CRTBP3D) {
                    columns.push_back({"vz", "-", "Dimensionless rotating-frame z velocity"});
                }
                columns.push_back({"CJ", "-", "Jacobi constant"});

                break;
                // -----------------------------------------------------------------
                // INDICATOR
                // -----------------------------------------------------------------
            case RunMode::INDICATOR:
                columns.push_back({"t", "day", "Physical time"});
                switch (init.getIndicator()) {
                    case Model::IndicatorType::FLI:
                        columns.push_back({"FLI", "-", "Fast Lyapunov Indicator"});
                        break;

                    case Model::IndicatorType::LCI:
                        columns.push_back({"LCI", "1/day", "Lyapunov Characteristic Indicator"});
                        break;

                    case Model::IndicatorType::RLI:
                        columns.push_back({"RLI", "-", "Relative Lyapunov Indicator"});
                        break;

                    case Model::IndicatorType::NONE:
                        throw std::runtime_error("INDICATOR output schema requires a chaos indicator.");

                    default:
                        throw std::runtime_error("Unknown chaos indicator while building output schema.");
                }

                break;

                // -----------------------------------------------------------------
                // GRID
                // -----------------------------------------------------------------
            case RunMode::GRID:
                for (const GridAxis &axis : init.getGridAxes()) {
                    columns.push_back({orbitalElementToString(axis.element), orbitalElementUnit(axis.element),
                                       orbitalElementDescription(axis.element)});
                }

                switch (init.getIndicator()) {
                    case Model::IndicatorType::FLI:
                        columns.push_back({"FLI", "-", "Fast Lyapunov Indicator"});
                        break;

                    case Model::IndicatorType::LCI:
                        columns.push_back({"LCI", "1/day", "Lyapunov Characteristic Indicator"});
                        break;

                    case Model::IndicatorType::RLI:
                        columns.push_back({"RLI", "-", "Relative Lyapunov Indicator"});
                        break;

                    case Model::IndicatorType::NONE:
                        throw std::runtime_error("GRID output schema requires a chaos indicator.");

                    default:
                        throw std::runtime_error("Unknown chaos indicator while building output schema.");
                }

                break;

            default:
                throw std::runtime_error("Unknown run mode while building output schema.");
        }

        return columns;
    }
}  // namespace io