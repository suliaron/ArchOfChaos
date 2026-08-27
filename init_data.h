#pragma once

#include "astro_types.h"  // astro::OrbitalElements
#include "model.h"        // Model::Formalism, Model::IndicatorType

#include <cstddef>   // std::size_t
#include <cstdint>   // std::uint32_t
#include <iostream>  // std::cout
#include <ostream>   // std::ostream
#include <string>    // std::string

/**
 * @brief Specifies the overall computation mode.
 */
enum class RunMode {
    ORBIT,     /**< Integrate a single orbit and save the state evolution. */
    INDICATOR, /**< Integrate a single orbit and save an indicator evolution. */
    GRID       /**< Integrate orbits over a parameter grid. */
};

/**
 * @brief Stores, parses, and validates program initialization data.
 *
 * Reads the initialization parameters from an input file, converts textual
 * values to the corresponding internal types, validates the resulting data,
 * and provides read-only access to the initialized parameters.
 *
 * Physical quantities in the input file use the following units:
 *
 * - masses: solar mass,
 * - distances and semimajor axes: AU,
 * - physical times: day,
 * - angular orbital elements: degree in the input file and radian internally.
 *
 * The initial epoch @c t0, integration duration @c T, output interval
 * @c output_dt, and pericenter passage time @c tau are specified in days.
 *
 * The CRTBP integration itself uses normalized dimensionless time; conversion
 * from physical time to normalized CRTBP time is performed outside this class.
 */
class InitData {
   public:
    /**
     * @brief Constructs initialization data from an input file.
     *
     * Reads, parses, and validates the specified initialization file.
     *
     * @param file_name Initialization file name or path.
     *
     * @throws std::runtime_error If the file cannot be opened or contains
     *         invalid initialization data.
     */
    explicit InitData(const std::string &file_name);

    /**
     * @brief Prints the initialization data.
     *
     * @param os Output stream. Defaults to std::cout.
     */
    void print(std::ostream &os = std::cout) const;

    /** @return Selected run mode. */
    RunMode getRunMode() const noexcept
    {
        return run_mode_;
    }

    /** @return Selected chaos indicator. */
    Model::IndicatorType getIndicator() const noexcept
    {
        return indicator_;
    }

    /** @return Selected CRTBP mathematical formulation. */
    Model::Formalism getFormalism() const noexcept
    {
        return formalism_;
    }

    /** @return Mass of the primary body [solar mass]. */
    double getM1() const noexcept
    {
        return m1_;
    }

    /** @return Mass of the secondary body [solar mass]. */
    double getM2() const noexcept
    {
        return m2_;
    }

    /** @return Constant distance between the two primary bodies [AU]. */
    double getA2() const noexcept
    {
        return a2_;
    }

    /**
     * @brief Returns the initial physical epoch.
     *
     * @return Initial epoch [day].
     */
    double getT0() const noexcept
    {
        return t0_;
    }

    /**
     * @brief Returns the physical integration duration.
     *
     * @return Integration duration [day].
     */
    double getT() const noexcept
    {
        return T_;
    }

    /**
     * @brief Returns the physical output time interval.
     *
     * @return Output time interval [day].
     */
    double getOutputDt() const noexcept
    {
        return output_dt_;
    }

    /** @return Initial osculating orbital elements. */
    const astro::OrbitalElements &getElements() const noexcept
    {
        return elements_;
    }

    /** @return Minimum semimajor axis of the grid. */
    double getA0() const noexcept
    {
        return a0_;
    }

    /** @return Maximum semimajor axis of the grid. */
    double getA1() const noexcept
    {
        return a1_;
    }

    /** @return Number of semimajor-axis grid intervals. */
    std::uint32_t getNa() const noexcept
    {
        return Na_;
    }

    /** @return Minimum eccentricity of the grid. */
    double getE0() const noexcept
    {
        return e0_;
    }

    /** @return Maximum eccentricity of the grid. */
    double getE1() const noexcept
    {
        return e1_;
    }

    /** @return Number of eccentricity grid intervals. */
    std::uint32_t getNe() const noexcept
    {
        return Ne_;
    }

    /**
     * @brief Returns the initial deviation vector.
     *
     * @return Pointer to the first element of the four-component
     *         deviation vector.
     */
    const double *getDy() const noexcept
    {
        return dy_;
    }

    /**
     * @brief Calculates the pericenter passage time.
     *
     * If the pericenter passage time was specified directly in the input file,
     * the stored value is returned. If the mean anomaly was specified instead,
     * the pericenter passage time is calculated from
     *
     *     tau = t0 - M / n,
     *
     * where
     *
     *     n = sqrt(mu_grav / a^3).
     *
     * The gravitational parameter is expressed in AU^3/day^2, the semimajor
     * axis in AU, and the resulting pericenter passage time in days.
     *
     * @param mu_grav Gravitational parameter of the Keplerian orbit [AU^3/day^2].
     * @param a Semimajor axis [AU].
     *
     * @return Pericenter passage time [day].
     *
     * @throws std::runtime_error If neither tau nor M was specified.
     */
    double calc_tau(double mu_grav, double a) const;

   private:
    /**
     * @brief Specifies how the initial orbital phase was provided.
     */
    enum class OrbitalPhaseInput {
        NONE,        /**< No orbital phase has been specified. */
        TAU,         /**< Pericenter passage time was specified. */
        MEAN_ANOMALY /**< Mean anomaly was specified. */
    };

    /**
     * @brief Returns the name of an orbital-phase input type.
     *
     * @param input Orbital-phase input type.
     * @return Name of the orbital-phase input type.
     */
    static const char *orbitalPhaseInputToString(OrbitalPhaseInput input) noexcept;

    /**
     * @brief Parses a single line of the initialization file.
     *
     * Empty lines and comments are ignored.
     *
     * @param line Input line.
     */
    void parseLine(const std::string &line);

    /**
     * @brief Validates the initialization data.
     *
     * @throws std::runtime_error If the stored input parameters are
     *         inconsistent or invalid.
     */
    void validate() const;

    /**
     * @brief Removes whitespace characters from a string.
     *
     * @param text String to modify.
     */
    static void removeSpaces(std::string &text);

    /**
     * @brief Converts text to a run mode.
     *
     * @param text Run mode as text.
     * @return Corresponding run mode.
     *
     * @throws std::runtime_error If the value is unknown.
     */
    static RunMode parseRunMode(const std::string &text);

    /**
     * @brief Converts text to an indicator type.
     *
     * @param text Indicator type as text.
     * @return Corresponding indicator type.
     *
     * @throws std::runtime_error If the value is unknown.
     */
    static Model::IndicatorType parseIndicatorType(const std::string &text);

    /**
     * @brief Converts text to a CRTBP mathematical formulation.
     *
     * @param text Formalism as text.
     * @return Corresponding CRTBP formulation.
     *
     * @throws std::runtime_error If the value is unknown.
     */
    static Model::Formalism parseFormalism(const std::string &text);

    RunMode              run_mode_  = RunMode::ORBIT;
    Model::Formalism     formalism_ = Model::Formalism::NEWTONIAN;
    Model::IndicatorType indicator_ = Model::IndicatorType::NONE;

    double m1_ = 0.0;
    double m2_ = 0.0;
    double a2_ = 0.0;

    /// Initial physical epoch [day].
    double t0_ = 0.0;
    /// Physical integration duration [day].
    double T_ = 0.0;
    /// Physical output time interval [day].
    double output_dt_ = 0.0;

    /// Specifies how the initial orbital phase was provided.
    OrbitalPhaseInput      orbital_phase_input_ = OrbitalPhaseInput::NONE;
    astro::OrbitalElements elements_{};
    /// Initial mean anomaly [rad].
    double mean_anomaly_ = 0.0;

    double        a0_ = 0.0;
    double        a1_ = 0.0;
    std::uint32_t Na_ = 0;

    double        e0_ = 0.0;
    double        e1_ = 0.0;
    std::uint32_t Ne_ = 0;

    double dy_[4] = {1.0, 0.0, 0.0, 0.0};
};