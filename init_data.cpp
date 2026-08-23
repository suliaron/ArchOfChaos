#include "init_data.h"

#include "math_utils.h"  // astro::torad

#include <algorithm>  // std::remove_if, std::transform
#include <cmath>      // std::abs
#include <cctype>     // std::isspace, std::toupper
#include <fstream>    // std::ifstream
#include <iomanip>    // std::setw
#include <sstream>    // std::istringstream
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string

namespace {

    /**
     * @brief Returns the name of a run mode.
     *
     * @param mode Run mode.
     * @return Name of the run mode.
     */
    const char *RunModeToString(RunMode mode) noexcept
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
}  // namespace

InitData::InitData(const std::string &file_name)
{
    std::ifstream file(file_name);

    if (!file) {
        throw std::runtime_error("Cannot open initialization file: " + file_name);
    }

    std::string line;

    while (std::getline(file, line)) {
        ParseLine(line);
    }

    Validate();
}

void InitData::RemoveSpaces(std::string &text)
{
    text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; }),
               text.end());
}

RunMode InitData::ParseRunMode(const std::string &text)
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

Model::IndicatorType InitData::ParseIndicatorType(const std::string &text)
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

Model::Formalism InitData::ParseFormalism(const std::string &text)
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

void InitData::ParseLine(const std::string &line)
{
    std::string text = line;

    // Remove comments.
    const std::size_t comment = text.find('#');

    if (comment != std::string::npos) {
        text.erase(comment);
    }

    // Remove whitespace.
    RemoveSpaces(text);

    // Ignore empty lines.
    if (text.empty()) {
        return;
    }

    // Split the line into key and value.
    const std::size_t pos = text.find('=');

    if (pos == std::string::npos) {
        throw std::runtime_error("Missing '=' in initialization file: " + line);
    }

    const std::string key   = text.substr(0, pos);
    const std::string value = text.substr(pos + 1);

    if (key.empty()) {
        throw std::runtime_error("Missing keyword in initialization file.");
    }

    if (value.empty()) {
        throw std::runtime_error("Missing value for '" + key + "'.");
    }

    // Parse enumeration values.
    if (key == "mode") {
        run_mode_ = ParseRunMode(value);
        return;
    }

    if (key == "indicator") {
        indicator_ = ParseIndicatorType(value);
        return;
    }

    if (key == "formalism") {
        formalism_ = ParseFormalism(value);
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
        is >> T_;

    } else if (key == "output_dt") {
        is >> output_dt_;

    } else if (key == "a") {
        is >> elements_.a;

    } else if (key == "e") {
        is >> elements_.e;

    } else if (key == "i") {
        double value_deg = 0.0;
        is >> value_deg;
        elements_.i = astro::toRad(value_deg);

    } else if (key == "omega") {
        double value_deg = 0.0;
        is >> value_deg;
        elements_.omega = astro::toRad(value_deg);

    } else if (key == "Omega") {
        double value_deg = 0.0;
        is >> value_deg;
        elements_.Omega = astro::toRad(value_deg);

    } else if (key == "tau") {
        is >> elements_.tau;

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

void InitData::Validate() const
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
    // Integration interval.
    if (T_ <= t0_) {
        throw std::runtime_error("T must be greater than t0.");
    }

    switch (run_mode_) {
        case RunMode::ORBIT:
            if (indicator_ != Model::IndicatorType::NONE) {
                throw std::runtime_error("ORBIT mode requires indicator = NONE.");
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
            if (std::abs(elements_.i) > PLANAR_EPS) {
                throw std::runtime_error("CRTBP2D requires inclination i = 0.");
            }
            break;

        case RunMode::INDICATOR:
            if (indicator_ == Model::IndicatorType::NONE) {
                throw std::runtime_error("INDICATOR mode requires an indicator.");
            }
            if (output_dt_ <= 0.0) {
                throw std::runtime_error("INDICATOR mode requires output_dt > 0.");
            }
            if (elements_.a <= 0.0) {
                throw std::runtime_error("Semimajor axis must be greater than zero.");
            }
            if (elements_.e < 0.0 || elements_.e >= 1.0) {
                throw std::runtime_error("Eccentricity must satisfy 0 <= e < 1.");
            }
            if (std::abs(elements_.i) > PLANAR_EPS) {
                throw std::runtime_error("CRTBP2D requires inclination i = 0.");
            }
            break;

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
    }
}

void InitData::Print(std::ostream &os) const
{
    constexpr int W = 18;

    os << "----------------------------------------\n";
    os << "Initialization data\n";
    os << "----------------------------------------\n";
    os << "run mode      : " << RunModeToString(run_mode_) << '\n';
    os << "indicator     : " << Model::IndicatorTypeToString(indicator_) << '\n';
    os << "formalism     : " << Model::FormalismToString(formalism_) << '\n';
    os << "m1            : " << std::setw(W) << m1_ << " [M_sun]\n";
    os << "m2            : " << std::setw(W) << m2_ << " [M_sun]\n";
    os << "a2            : " << std::setw(W) << a2_ << " [AU]\n";
    os << "t0            : " << std::setw(W) << t0_ << '\n';
    os << "T             : " << std::setw(W) << T_ << '\n';

    if (run_mode_ == RunMode::ORBIT || run_mode_ == RunMode::INDICATOR) {
        os << "output_dt     : " << std::setw(W) << output_dt_ << '\n';
        os << "a             : " << std::setw(W) << elements_.a << " [AU]\n";
        os << "e             : " << std::setw(W) << elements_.e << '\n';
        os << "i             : " << std::setw(W) << elements_.i << " [rad]\n";
        os << "omega         : " << std::setw(W) << elements_.omega << " [rad]\n";
        os << "Omega         : " << std::setw(W) << elements_.Omega << " [rad]\n";
        os << "tau           : " << std::setw(W) << elements_.tau << " [JD]\n";
    }

    if (run_mode_ == RunMode::GRID) {
        os << "a0            : " << std::setw(W) << a0_ << '\n';
        os << "a1            : " << std::setw(W) << a1_ << '\n';
        os << "Na            : " << std::setw(W) << Na_ << '\n';
        os << "e0            : " << std::setw(W) << e0_ << '\n';
        os << "e1            : " << std::setw(W) << e1_ << '\n';
        os << "Ne            : " << std::setw(W) << Ne_ << '\n';
    }

    if (indicator_ != Model::IndicatorType::NONE) {
        os << "dy1           : " << std::setw(W) << dy_[0] << '\n';
        os << "dy2           : " << std::setw(W) << dy_[1] << '\n';
        os << "dy3           : " << std::setw(W) << dy_[2] << '\n';
        os << "dy4           : " << std::setw(W) << dy_[3] << '\n';
    }
}