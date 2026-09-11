#include "init_data.h"

#include "math_utils.h"  // astro::torad

#include <algorithm>  // std::remove_if, std::transform std::find_if_not, std::any_of
#include <cmath>      // std::abs, std::sqrt
#include <cctype>     // std::isspace, std::toupper
#include <fstream>    // std::ifstream
#include <iomanip>    // std::setw
#include <sstream>    // std::istringstream
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string

/**
 * @brief Returns the name of a run mode.
 *
 * @param mode Run mode.
 * @return Name of the run mode.
 */
const char *runModeToString(RunMode mode) noexcept
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

InitData::InitData(const std::string &file_name)
{
    std::ifstream file(file_name);
    if (!file) {
        throw std::runtime_error("Cannot open initialization file: " + file_name);
    }

    std::string line;
    while (std::getline(file, line)) {
        parseLine(line);
    }

    validate();
}

void InitData::removeSpaces(std::string &text)
{
    text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; }),
               text.end());
}

void InitData::Trim(std::string &text)
{
    const auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; });
    const auto last =
        std::find_if_not(text.rbegin(), text.rend(), [](unsigned char c) { return std::isspace(c) != 0; }).base();

    if (first >= last) {
        text.clear();
        return;
    }

    text = std::string(first, last);
}

RunMode InitData::parseRunMode(const std::string &text)
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

Model::IndicatorType InitData::parseIndicatorType(const std::string &text)
{
    std::string indicator(text);

    std::transform(indicator.begin(), indicator.end(), indicator.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (indicator == "NONE") {
        return Model::IndicatorType::NONE;
    }
    if (indicator == "FLI") {
        return Model::IndicatorType::FLI;
    }
    if (indicator == "LCI") {
        return Model::IndicatorType::LCI;
    }
    if (indicator == "RLI") {
        return Model::IndicatorType::RLI;
    }

    throw std::runtime_error("Unknown indicator type: " + text);
}

Model::Formalism InitData::parseFormalism(const std::string &text)
{
    std::string formalism(text);

    std::transform(formalism.begin(), formalism.end(), formalism.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (formalism == "NEWTONIAN") {
        return Model::Formalism::NEWTONIAN;
    }

    if (formalism == "HAMILTONIAN") {
        return Model::Formalism::HAMILTONIAN;
    }

    throw std::runtime_error("Unknown mathematical formalism: " + text);
}

void InitData::parseLine(const std::string &line)
{
    std::string text = line;

    // Remove comments.
    const std::size_t comment = text.find('#');

    if (comment != std::string::npos) {
        text.erase(comment);
    }

    // Remove leading and trailing whitespace.
    Trim(text);

    // Ignore empty lines.
    if (text.empty()) {
        return;
    }

    // Split the line into key and value.
    const std::size_t pos = text.find('=');

    if (pos == std::string::npos) {
        throw std::runtime_error("Missing '=' in initialization file: " + line);
    }

    std::string key   = text.substr(0, pos);
    std::string value = text.substr(pos + 1);

    // Remove whitespace surrounding the key and value.
    Trim(key);
    Trim(value);

    if (key.empty()) {
        throw std::runtime_error("Missing keyword in initialization file.");
    }

    if (value.empty()) {
        throw std::runtime_error("Missing value for '" + key + "'.");
    }

    // Parse enumeration values.
    if (key == "mode") {
        run_mode_ = parseRunMode(value);
        return;
    }

    if (key == "indicator") {
        indicator_ = parseIndicatorType(value);
        return;
    }

    if (key == "formalism") {
        formalism_ = parseFormalism(value);
        return;
    }

    if (key == "grid") {
        std::istringstream is(value);
        std::string        element_name;
        GridAxis           axis{};

        is >> element_name >> axis.min >> axis.max >> axis.nIntervals;

        if (!is) {
            throw std::runtime_error("Invalid grid definition: " + value);
        }
        is >> std::ws;
        if (!is.eof()) {
            throw std::runtime_error("Invalid grid definition: " + value);
        }

        axis.element = orbitalElementFromString(element_name);
        grid_axes_.push_back(axis);
        return;
    }

    // Parse numerical values.
    std::istringstream is(value);

    if (key == "m1") {
        is >> m1_;

    } else if (key == "m2") {
        is >> m2_;

    } else if (key == "a2") {
        is >> a2_;

    } else if (key == "t0") {
        is >> t0_;

    } else if (key == "T") {
        if (integration_duration_input_ == IntegrationDurationInput::ORBITAL_PERIODS) {
            throw std::runtime_error("T and nPeriods cannot be specified simultaneously.");
        }
        is >> T_;
        integration_duration_input_ = IntegrationDurationInput::PHYSICAL_TIME;
    } else if (key == "nPeriods") {
        if (integration_duration_input_ == IntegrationDurationInput::PHYSICAL_TIME) {
            throw std::runtime_error("T and nPeriods cannot be specified simultaneously.");
        }
        is >> n_periods_;
        integration_duration_input_ = IntegrationDurationInput::ORBITAL_PERIODS;
    } else if (key == "output_dt") {
        is >> output_dt_;
    } else if (key == "output_first") {
        is >> output_first_;
    } else if (key == "output_points_per_decade") {
        is >> output_points_per_decade_;
    } else if (key == "a") {
        registerFixedOrbitalElement(OrbitalElement::SEMIMAJOR_AXIS);
        is >> elements_.a;
    } else if (key == "e") {
        registerFixedOrbitalElement(OrbitalElement::ECCENTRICITY);
        is >> elements_.e;
    } else if (key == "i") {
        registerFixedOrbitalElement(OrbitalElement::INCLINATION);
        double value_deg = 0.0;
        is >> value_deg;
        elements_.i = astro::toRad(value_deg);
    } else if (key == "omega") {
        registerFixedOrbitalElement(OrbitalElement::ARGUMENT_OF_PERICENTER);
        double value_deg = 0.0;
        is >> value_deg;
        elements_.omega = astro::toRad(value_deg);
    } else if (key == "Omega") {
        registerFixedOrbitalElement(OrbitalElement::LONGITUDE_OF_ASCENDING_NODE);
        double value_deg = 0.0;
        is >> value_deg;
        elements_.Omega = astro::toRad(value_deg);
    } else if (key == "tau") {
        if (orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY) {
            throw std::runtime_error("Both tau and M are specified. Use only one.");
        }
        registerFixedOrbitalElement(OrbitalElement::PERICENTER_TIME);
        is >> elements_.tau;
        orbital_phase_input_ = OrbitalPhaseInput::TAU;
    } else if (key == "M") {
        if (orbital_phase_input_ == OrbitalPhaseInput::TAU) {
            throw std::runtime_error("Both tau and M are specified. Use only one.");
        }
        registerFixedOrbitalElement(OrbitalElement::MEAN_ANOMALY);
        double value_deg = 0.0;
        is >> value_deg;
        mean_anomaly_        = astro::toRad(value_deg);
        orbital_phase_input_ = OrbitalPhaseInput::MEAN_ANOMALY;
    } else if (key == "a0") {
        is >> a0_;
    } else if (key == "a1") {
        is >> a1_;
    } else if (key == "Na") {
        is >> Na_;
    } else if (key == "e0") {
        is >> e0_;
    } else if (key == "e1") {
        is >> e1_;
    } else if (key == "Ne") {
        is >> Ne_;
    } else if (key == "dy1") {
        is >> dy_[0];
    } else if (key == "dy2") {
        is >> dy_[1];
    } else if (key == "dy3") {
        is >> dy_[2];
    } else if (key == "dy4") {
        is >> dy_[3];
    } else if (key == "relTol") {
        is >> rel_tol_;
    } else if (key == "absTol") {
        is >> abs_tol_;
    } else {
        throw std::runtime_error("Unknown keyword: " + key);
    }

    // Check that the complete value was parsed successfully.
    if (!is) {
        throw std::runtime_error("Invalid value for '" + key + "'.");
    }

    is >> std::ws;

    if (!is.eof()) {
        throw std::runtime_error("Invalid value for '" + key + "'.");
    }
}

double InitData::calc_tau(double mu_grav, double a) const
{
    if (orbital_phase_input_ == OrbitalPhaseInput::TAU) {
        return elements_.tau;
    }

    if (orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY) {
        const double n = std::sqrt(mu_grav / astro::cube(a));
        return (t0_ - mean_anomaly_ / n);
    }

    throw std::runtime_error("Neither tau nor M has been specified.");
}

double InitData::calcIntegrationDuration(double mu_13, double a) const
{
    switch (integration_duration_input_) {
        case IntegrationDurationInput::PHYSICAL_TIME:
            return T_;

        case IntegrationDurationInput::ORBITAL_PERIODS: {
            const double period = 2.0 * astro::pi * std::sqrt(astro::cube(a) / mu_13);

            return n_periods_ * period;
        }

        case IntegrationDurationInput::NONE:
            throw std::runtime_error("Integration duration has not been specified.");
    }

    throw std::runtime_error("Unknown integration-duration input method.");
}

bool InitData::hasGridAxis(OrbitalElement element) const noexcept
{
    return std::any_of(grid_axes_.begin(), grid_axes_.end(),
                       [element](const GridAxis &axis) { return axis.element == element; });
}

bool InitData::hasFixedOrbitalElement(OrbitalElement element) const noexcept
{
    return fixed_elements_.find(element) != fixed_elements_.end();
}

void InitData::registerFixedOrbitalElement(OrbitalElement element)
{
    const bool inserted = fixed_elements_.insert(element).second;

    if (!inserted) {
        throw std::runtime_error("Duplicate fixed orbital element: " + std::string(orbitalElementToString(element)));
    }
}

void InitData::validateOrbitalElementSource(OrbitalElement element) const
{
    const bool fixed = hasFixedOrbitalElement(element);
    const bool grid  = hasGridAxis(element);

    if (fixed && grid) {
        throw std::runtime_error("Orbital element '" + std::string(orbitalElementToString(element)) +
                                 "' cannot be specified both as a fixed value and as a grid axis.");
    }

    if (!fixed && !grid) {
        throw std::runtime_error("Orbital element '" + std::string(orbitalElementToString(element)) +
                                 "' must be specified either as a fixed value or as a grid axis.");
    }
}

void InitData::validate() const
{
    constexpr double PLANAR_EPS = 1.0e-12;

    // Primary-system parameters.
    if (m1_ <= 0.0) {
        throw std::runtime_error("m1 must be greater than zero.");
    }
    if (m2_ <= 0.0) {
        throw std::runtime_error("m2 must be greater than zero.");
    }
    if (a2_ <= 0.0) {
        throw std::runtime_error("a2 must be greater than zero.");
    }
    // Integration duration.
    switch (integration_duration_input_) {
        case IntegrationDurationInput::PHYSICAL_TIME:
            if (T_ <= 0.0) {
                throw std::runtime_error("Integration duration T must be greater than zero.");
            }
            break;

        case IntegrationDurationInput::ORBITAL_PERIODS:
            if (n_periods_ <= 0.0) {
                throw std::runtime_error("nPeriods must be greater than zero.");
            }
            break;

        case IntegrationDurationInput::NONE:
            throw std::runtime_error("Either T or nPeriods must be specified.");
    }

    if (rel_tol_ <= 0.0) {
        throw std::runtime_error("Relative tolerance relTol must be greater than zero.");
    }

    if (abs_tol_ <= 0.0) {
        throw std::runtime_error("Absolute tolerance absTol must be greater than zero.");
    }

    switch (run_mode_) {
        case RunMode::ORBIT:
            if (indicator_ != Model::IndicatorType::NONE) {
                throw std::runtime_error("ORBIT mode requires indicator = NONE.");
            }

            if (orbital_phase_input_ == OrbitalPhaseInput::NONE) {
                throw std::runtime_error("ORBIT mode requires either tau or M.");
            }

            if (output_dt_ <= 0.0) {
                throw std::runtime_error("ORBIT mode requires output_dt > 0.");
            }

            if (elements_.a <= 0.0) {
                throw std::runtime_error("Semimajor axis must be greater than zero.");
            }

            if (elements_.e < 0.0 || elements_.e >= 1.0) {
                throw std::runtime_error("Eccentricity must satisfy 0 <= e < 1.");
            }

            if (elements_.i < 0.0 || elements_.i >= PLANAR_EPS) {
                throw std::runtime_error("CRTBP2D requires inclination 0 <= i < PLANAR_EPS.");
            }

            break;

        case RunMode::INDICATOR:
            // INDICATOR mode requires a chaos indicator.
            if (indicator_ == Model::IndicatorType::NONE) {
                throw std::runtime_error("INDICATOR mode requires an indicator.");
            }

            // The first logarithmic output time is an elapsed physical time
            // measured from the initial epoch and must therefore be positive.
            if (output_first_ <= 0.0) {
                throw std::runtime_error("INDICATOR mode requires output_first > 0.");
            }

            // At least one output point must be written in each time decade.
            if (output_points_per_decade_ == 0) {
                throw std::runtime_error("INDICATOR mode requires output_points_per_decade > 0.");
            }

            // Validate the initial osculating orbit.
            if (elements_.a <= 0.0) {
                throw std::runtime_error("Semimajor axis must be greater than zero.");
            }

            if (elements_.e < 0.0 || elements_.e >= 1.0) {
                throw std::runtime_error("Eccentricity must satisfy 0 <= e < 1.");
            }

            if (elements_.i < 0.0 || elements_.i >= PLANAR_EPS) {
                throw std::runtime_error("CRTBP2D requires inclination 0 <= i < PLANAR_EPS.");
            }

            break;

        case RunMode::GRID: {
            if (indicator_ == Model::IndicatorType::NONE) {
                throw std::runtime_error("GRID mode requires an indicator.");
            }

            // At least one grid axis must be specified.
            if (grid_axes_.empty()) {
                throw std::runtime_error("GRID mode requires at least one grid axis.");
            }

            // A physical orbit has six independent orbital elements.
            // The orbital phase may be represented either by tau or by M.
            if (grid_axes_.size() > 6) {
                throw std::runtime_error("GRID mode supports at most six independent grid axes.");
            }

            // Check for duplicate grid axes.
            for (std::size_t i = 0; i < grid_axes_.size(); ++i) {
                for (std::size_t j = i + 1; j < grid_axes_.size(); ++j) {
                    if (grid_axes_[i].element == grid_axes_[j].element) {
                        throw std::runtime_error("Duplicate grid axis: " +
                                                 std::string(orbitalElementToString(grid_axes_[i].element)));
                    }
                }
            }

            // Validate every grid axis.
            for (const GridAxis &axis : grid_axes_) {
                // Check that each geometric orbital element is specified exactly once:
                // either as a fixed input value or as a grid axis.
                validateOrbitalElementSource(OrbitalElement::SEMIMAJOR_AXIS);
                validateOrbitalElementSource(OrbitalElement::ECCENTRICITY);
                validateOrbitalElementSource(OrbitalElement::INCLINATION);
                validateOrbitalElementSource(OrbitalElement::ARGUMENT_OF_PERICENTER);
                validateOrbitalElementSource(OrbitalElement::LONGITUDE_OF_ASCENDING_NODE);

                if (axis.nIntervals == 0) {
                    throw std::runtime_error("Grid axis '" + std::string(orbitalElementToString(axis.element)) +
                                             "' requires nIntervals > 0.");
                }

                if (axis.max <= axis.min) {
                    throw std::runtime_error("Grid axis '" + std::string(orbitalElementToString(axis.element)) +
                                             "' requires max > min.");
                }

                switch (axis.element) {
                    case OrbitalElement::SEMIMAJOR_AXIS:
                        if (axis.min <= 0.0) {
                            throw std::runtime_error("Semimajor-axis grid requires a > 0.");
                        }
                        break;

                    case OrbitalElement::ECCENTRICITY:
                        if (axis.min < 0.0 || axis.max >= 1.0) {
                            throw std::runtime_error("Eccentricity grid requires 0 <= e < 1.");
                        }
                        break;

                    case OrbitalElement::INCLINATION:
                        // Grid angular values are still stored in degrees here.
                        if (axis.min < 0.0 || axis.max >= PLANAR_EPS) {
                            throw std::runtime_error("CRTBP2D currently requires inclination 0 <= i < PLANAR_EPS.");
                        }
                        break;

                    case OrbitalElement::ARGUMENT_OF_PERICENTER:
                        if (axis.min < 0.0 || axis.max > 360.0) {
                            throw std::runtime_error("Argument-of-pericenter grid requires 0 <= omega <= 360 deg.");
                        }
                        break;

                    case OrbitalElement::LONGITUDE_OF_ASCENDING_NODE:
                        if (axis.min < 0.0 || axis.max > 360.0) {
                            throw std::runtime_error(
                                "Longitude-of-ascending-node grid requires 0 <= Omega <= 360 deg.");
                        }
                        break;

                    case OrbitalElement::MEAN_ANOMALY:
                        if (axis.min < 0.0 || axis.max > 360.0) {
                            throw std::runtime_error("Mean-anomaly grid requires 0 <= M <= 360 deg.");
                        }
                        break;

                    case OrbitalElement::PERICENTER_TIME:
                        break;
                }
            }

            const bool hasA   = hasGridAxis(OrbitalElement::SEMIMAJOR_AXIS);
            const bool hasE   = hasGridAxis(OrbitalElement::ECCENTRICITY);
            const bool hasI   = hasGridAxis(OrbitalElement::INCLINATION);
            const bool hasTau = hasGridAxis(OrbitalElement::PERICENTER_TIME);
            const bool hasM   = hasGridAxis(OrbitalElement::MEAN_ANOMALY);

            // tau and M are alternative representations of the orbital phase.
            if (hasTau && hasM) {
                throw std::runtime_error("tau and M cannot both be used as grid axes.");
            }

            // The orbital phase cannot be both fixed and grid-controlled.
            if ((hasTau || hasM) && orbital_phase_input_ != OrbitalPhaseInput::NONE) {
                throw std::runtime_error(
                    "Orbital phase cannot be specified both as a fixed "
                    "input value and as a grid axis.");
            }

            // If the phase is not controlled by the grid, a fixed phase
            // must be specified.
            if (!hasTau && !hasM && orbital_phase_input_ == OrbitalPhaseInput::NONE) {
                throw std::runtime_error(
                    "GRID mode requires tau or M either as a fixed "
                    "orbital-phase input or as a grid axis.");
            }

            // Validate orbital elements that are not controlled by the grid.
            if (!hasA && elements_.a <= 0.0) {
                throw std::runtime_error("Semimajor axis must be greater than zero.");
            }

            if (!hasE && (elements_.e < 0.0 || elements_.e >= 1.0)) {
                throw std::runtime_error("Eccentricity must satisfy 0 <= e < 1.");
            }

            if (!hasI && (elements_.i < 0.0 || elements_.i >= PLANAR_EPS)) {
                throw std::runtime_error("CRTBP2D requires inclination 0 <= i < PLANAR_EPS.");
            }

            break;
        }
            /*
        case RunMode::GRID:
            if (indicator_ == Model::IndicatorType::NONE) {
                throw std::runtime_error("GRID mode requires an indicator.");
            }
            if (a0_ <= 0.0) {
                throw std::runtime_error("GRID mode requires a0 > 0.");
            }
            if (a1_ < a0_) {
                throw std::runtime_error("GRID mode requires a1 >= a0.");
            }
            if (Na_ == 0) {
                throw std::runtime_error("GRID mode requires Na > 0.");
            }
            if (e0_ < 0.0 || e1_ >= 1.0) {
                throw std::runtime_error("GRID mode requires 0 <= e0 <= e1 < 1.");
            }
            if (e1_ < e0_) {
                throw std::runtime_error("GRID mode requires e1 >= e0.");
            }
            if (Ne_ == 0) {
                throw std::runtime_error("GRID mode requires Ne > 0.");
            }
            break;
            */
    }
}

const char *InitData::orbitalPhaseInputToString(OrbitalPhaseInput input) noexcept
{
    switch (input) {
        case OrbitalPhaseInput::NONE:
            return "NONE";

        case OrbitalPhaseInput::TAU:
            return "TAU";

        case OrbitalPhaseInput::MEAN_ANOMALY:
            return "MEAN_ANOMALY";
    }

    return "UNKNOWN";
}

const char *InitData::integrationDurationInputToString(IntegrationDurationInput input) noexcept
{
    switch (input) {
        case IntegrationDurationInput::NONE:
            return "NONE";

        case IntegrationDurationInput::PHYSICAL_TIME:
            return "PHYSICAL_TIME";

        case IntegrationDurationInput::ORBITAL_PERIODS:
            return "ORBITAL_PERIODS";
    }

    return "UNKNOWN";
}

void InitData::print(std::ostream &os) const
{
    constexpr int W = 18;

    os << "----------------------------------------\n";
    os << "Initialization data\n";
    os << "----------------------------------------\n";

    os << "run mode      : " << runModeToString(run_mode_) << '\n';
    os << "indicator     : " << Model::indicatorTypeToString(indicator_) << '\n';
    os << "formalism     : " << Model::formalismToString(formalism_) << '\n';

    // Orbital-phase input.
    if (run_mode_ == RunMode::GRID) {
        if (hasGridAxis(OrbitalElement::MEAN_ANOMALY)) {
            os << "phase input   : GRID(M)\n";
        } else if (hasGridAxis(OrbitalElement::PERICENTER_TIME)) {
            os << "phase input   : GRID(tau)\n";
        } else {
            os << "phase input   : " << orbitalPhaseInputToString(orbital_phase_input_) << '\n';
        }
    } else {
        os << "phase input   : " << orbitalPhaseInputToString(orbital_phase_input_) << '\n';
    }

    // Numerical integration tolerances.
    os << "relTol        : " << std::setw(W) << rel_tol_ << '\n';
    os << "absTol        : " << std::setw(W) << abs_tol_ << '\n';

    // Primary-system parameters.
    os << "m1            : " << std::setw(W) << m1_ << " [M_sun]\n";
    os << "m2            : " << std::setw(W) << m2_ << " [M_sun]\n";
    os << "a2            : " << std::setw(W) << a2_ << " [AU]\n";

    // Initial physical epoch.
    os << "t0            : " << std::setw(W) << t0_ << " [day]\n";

    // Integration-duration input.
    os << "duration input: " << integrationDurationInputToString(integration_duration_input_) << '\n';

    switch (integration_duration_input_) {
        case IntegrationDurationInput::PHYSICAL_TIME:
            os << "T             : " << std::setw(W) << T_ << " [day]\n";
            break;

        case IntegrationDurationInput::ORBITAL_PERIODS:
            os << "nPeriods      : " << std::setw(W) << n_periods_ << '\n';
            break;

        case IntegrationDurationInput::NONE:
            break;
    }

    // ---------------------------------------------------------------------
    // ORBIT and INDICATOR modes
    // ---------------------------------------------------------------------

    if (run_mode_ == RunMode::ORBIT || run_mode_ == RunMode::INDICATOR) {
        os << "output_dt     : " << std::setw(W) << output_dt_ << " [day]\n";

        // Orbital elements.
        os << "a             : " << std::setw(W) << elements_.a << " [AU]\n";
        os << "e             : " << std::setw(W) << elements_.e << '\n';
        os << "i             : " << std::setw(W) << elements_.i << " [rad]\n";
        os << "omega         : " << std::setw(W) << elements_.omega << " [rad]\n";
        os << "Omega         : " << std::setw(W) << elements_.Omega << " [rad]\n";
        if (orbital_phase_input_ == OrbitalPhaseInput::TAU) {
            os << "tau           : " << std::setw(W) << elements_.tau << " [day]\n";
        } else if (orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY) {
            os << "M             : " << std::setw(W) << mean_anomaly_ << " [rad]\n";
        }
    }

    // ---------------------------------------------------------------------
    // GRID mode
    // ---------------------------------------------------------------------

    if (run_mode_ == RunMode::GRID) {
        // Fixed orbital elements not controlled by the grid.
        if (!hasGridAxis(OrbitalElement::SEMIMAJOR_AXIS)) {
            os << "a             : " << std::setw(W) << elements_.a << " [AU]\n";
        }
        if (!hasGridAxis(OrbitalElement::ECCENTRICITY)) {
            os << "e             : " << std::setw(W) << elements_.e << '\n';
        }
        if (!hasGridAxis(OrbitalElement::INCLINATION)) {
            os << "i             : " << std::setw(W) << elements_.i << " [rad]\n";
        }
        if (!hasGridAxis(OrbitalElement::ARGUMENT_OF_PERICENTER)) {
            os << "omega         : " << std::setw(W) << elements_.omega << " [rad]\n";
        }
        if (!hasGridAxis(OrbitalElement::LONGITUDE_OF_ASCENDING_NODE)) {
            os << "Omega         : " << std::setw(W) << elements_.Omega << " [rad]\n";
        }
        // Fixed orbital phase, if the phase is not controlled by the grid.
        if (!hasGridAxis(OrbitalElement::PERICENTER_TIME) && !hasGridAxis(OrbitalElement::MEAN_ANOMALY)) {
            if (orbital_phase_input_ == OrbitalPhaseInput::TAU) {
                os << "tau           : " << std::setw(W) << elements_.tau << " [day]\n";

            } else if (orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY) {
                os << "M             : " << std::setw(W) << mean_anomaly_ << " [rad]\n";
            }
        }

        // Grid-axis definitions.
        os << "grid axes     : " << grid_axes_.size() << '\n';

        os << std::left << std::setw(16) << " " << std::setw(10) << "element" << std::right << std::setw(W) << "min"
           << std::setw(W) << "max" << std::setw(W) << "intervals" << std::setw(10) << "unit" << '\n';

        for (const GridAxis &axis : grid_axes_) {
            os << std::left << std::setw(14) << "grid" << std::setw(10) << ":" << orbitalElementToString(axis.element)
               << std::right << std::setw(W) << axis.min << std::setw(W) << axis.max << std::setw(W) << axis.nIntervals
               << std::setw(10) << orbitalElementUnit(axis.element) << '\n';
        }
    }

    // Initial deviation vector.
    if (indicator_ != Model::IndicatorType::NONE) {
        os << "dy1           : " << std::setw(W) << dy_[0] << '\n';
        os << "dy2           : " << std::setw(W) << dy_[1] << '\n';
        os << "dy3           : " << std::setw(W) << dy_[2] << '\n';
        os << "dy4           : " << std::setw(W) << dy_[3] << '\n';
    }
}

// void InitData::print(std::ostream &os) const
//{
//     constexpr int W = 18;
//
//     os << "----------------------------------------\n";
//     os << "Initialization data\n";
//     os << "----------------------------------------\n";
//     os << "run mode      : " << runModeToString(run_mode_) << '\n';
//     os << "indicator     : " << Model::indicatorTypeToString(indicator_) << '\n';
//     os << "formalism     : " << Model::formalismToString(formalism_) << '\n';
//     os << "phase input   : " << orbitalPhaseInputToString(orbital_phase_input_) << '\n';
//     os << "relTol        : " << std::setw(W) << rel_tol_ << '\n';
//     os << "absTol        : " << std::setw(W) << abs_tol_ << '\n';
//     os << "m1            : " << std::setw(W) << m1_ << " [M_sun]\n";
//     os << "m2            : " << std::setw(W) << m2_ << " [M_sun]\n";
//     os << "a2            : " << std::setw(W) << a2_ << " [AU]\n";
//     os << "t0            : " << std::setw(W) << t0_ << " [day]\n";
//     os << "duration input: " << integrationDurationInputToString(integration_duration_input_) << '\n';
//
//     switch (integration_duration_input_) {
//         case IntegrationDurationInput::PHYSICAL_TIME:
//             os << "T             : " << std::setw(W) << T_ << " [day]\n";
//             break;
//
//         case IntegrationDurationInput::ORBITAL_PERIODS:
//             os << "nPeriods      : " << std::setw(W) << n_periods_ << '\n';
//             break;
//
//         case IntegrationDurationInput::NONE:
//             break;
//     }
//
//     if (run_mode_ == RunMode::ORBIT || run_mode_ == RunMode::INDICATOR) {
//         os << "output_dt     : " << std::setw(W) << output_dt_ << " [day]\n";
//
//         // Orbital elements.
//         os << "a             : " << std::setw(W) << elements_.a << " [AU]\n";
//         os << "e             : " << std::setw(W) << elements_.e << '\n';
//         os << "i             : " << std::setw(W) << elements_.i << " [rad]\n";
//         os << "omega         : " << std::setw(W) << elements_.omega << " [rad]\n";
//         os << "Omega         : " << std::setw(W) << elements_.Omega << " [rad]\n";
//
//         if (orbital_phase_input_ == OrbitalPhaseInput::TAU) {
//             os << "tau           : " << std::setw(W) << elements_.tau << " [day]\n";
//         } else if (orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY) {
//             os << "M             : " << std::setw(W) << mean_anomaly_ << " [rad]\n";
//         }
//     }
//
//     if (run_mode_ == RunMode::GRID) {
//         os << "grid axes     : " << grid_axes_.size() << '\n';
//
//         os << std::left << std::setw(16) << " " << std::setw(10) << "element" << std::right << std::setw(W) << "min"
//            << std::setw(W) << "max" << std::setw(W) << "intervals" << '\n';
//
//         for (const GridAxis &axis : grid_axes_) {
//             os << std::left << std::setw(16) << "grid" << std::setw(10) << orbitalElementToString(axis.element)
//                << std::right << std::setw(W) << axis.min << std::setw(W) << axis.max << std::setw(W) <<
//                axis.nIntervals
//                << '\n';
//         }
//     }
//
//     // if (run_mode_ == RunMode::GRID) {
//     //     // Grid parameters.
//     //     os << "a0            : " << std::setw(W) << a0_ << " [AU]\n";
//     //     os << "a1            : " << std::setw(W) << a1_ << " [AU]\n";
//     //     os << "Na            : " << std::setw(W) << Na_ << '\n';
//     //     os << "e0            : " << std::setw(W) << e0_ << '\n';
//     //     os << "e1            : " << std::setw(W) << e1_ << '\n';
//     //     os << "Ne            : " << std::setw(W) << Ne_ << '\n';
//
//     //    // Fixed orbital elements.
//     //    os << "i             : " << std::setw(W) << elements_.i << " [rad]\n";
//     //    os << "omega         : " << std::setw(W) << elements_.omega << " [rad]\n";
//     //    os << "Omega         : " << std::setw(W) << elements_.Omega << " [rad]\n";
//
//     //    if (orbital_phase_input_ == OrbitalPhaseInput::TAU) {
//     //        os << "tau           : " << std::setw(W) << elements_.tau << '\n';
//     //    } else if (orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY) {
//     //        os << "M             : " << std::setw(W) << mean_anomaly_ << " [rad]\n";
//     //    }
//     //}
//
//     if (indicator_ != Model::IndicatorType::NONE) {
//         os << "dy1           : " << std::setw(W) << dy_[0] << '\n';
//         os << "dy2           : " << std::setw(W) << dy_[1] << '\n';
//         os << "dy3           : " << std::setw(W) << dy_[2] << '\n';
//         os << "dy4           : " << std::setw(W) << dy_[3] << '\n';
//     }
// }
