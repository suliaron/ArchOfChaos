#pragma once

#include "astro_types.h"  // astro::OrbitalElements
#include "grid.h"         // GridAxis, OrbitalElement, and grid-related declarations.
#include "model.h"        // Model::Formalism, Model::IndicatorType

#include <cstddef>   // std::size_t
#include <cstdint>   // std::uint32_t
#include <iostream>  // std::cout
#include <ostream>   // std::ostream
#include <set>       // std::set container for explicitly specified fixed orbital elements.
#include <string>    // std::string
#include <vector>    // std::vector container for storing grid-axis definitions.

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

    /**
     * @brief Checks whether the integration duration is specified in days.
     *
     * @return True if T is used to specify the integration duration.
     */
    bool usesPhysicalIntegrationTime() const noexcept
    {
        return integration_duration_input_ == IntegrationDurationInput::PHYSICAL_TIME;
    }

    /**
     * @brief Checks whether the integration duration is specified by
     *        the number of P3 orbital periods.
     *
     * @return True if nPeriods is used to specify the integration duration.
     */
    bool usesOrbitalPeriods() const noexcept
    {
        return integration_duration_input_ == IntegrationDurationInput::ORBITAL_PERIODS;
    }

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
     * @brief Returns the requested number of initial P3 orbital periods.
     *
     * @return Number of orbital periods.
     */
    double getNPeriods() const noexcept
    {
        return n_periods_;
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

    /**
     * @brief Returns the fixed mean anomaly.
     *
     * @return Mean anomaly [rad].
     */
    double getMeanAnomaly() const noexcept
    {
        return mean_anomaly_;
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
     * @brief Checks whether the fixed orbital phase is specified
     *        by the time of pericenter passage.
     *
     * @return True if tau is used as the fixed orbital-phase input.
     */
    bool usesFixedTau() const noexcept
    {
        return orbital_phase_input_ == OrbitalPhaseInput::TAU;
    }

    /**
     * @brief Checks whether the fixed orbital phase is specified
     *        by the mean anomaly.
     *
     * @return True if M is used as the fixed orbital-phase input.
     */
    bool usesFixedMeanAnomaly() const noexcept
    {
        return orbital_phase_input_ == OrbitalPhaseInput::MEAN_ANOMALY;
    }

    /**
     * @brief Returns the relative integration tolerance.
     *
     * @return Relative tolerance of the adaptive numerical integrator.
     */
    double getRelTol() const noexcept
    {
        return rel_tol_;
    }

    /**
     * @brief Returns the absolute integration tolerance.
     *
     * @return Absolute tolerance of the adaptive numerical integrator.
     */
    double getAbsTol() const noexcept
    {
        return abs_tol_;
    }

    /**
     * @brief Returns the orbital-element grid definitions.
     *
     * @return Constant reference to the grid-axis definitions.
     */
    const std::vector<GridAxis> &getGridAxes() const noexcept
    {
        return grid_axes_;
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

    /**
     * @brief Calculates the physical integration duration for a given orbit.
     *
     * If the integration duration is specified directly by T, the stored
     * physical duration is returned.
     *
     * If the integration duration is specified by nPeriods, the Keplerian
     * orbital period of P3 is calculated from the semimajor axis @p a supplied
     * for the current orbit:
     *
     *     P3(a) = 2*pi*sqrt(a^3 / mu_13).
     *
     * The integration duration is then set to
     *
     *     T(a) = nPeriods * P3(a).
     *
     * Thus, in GRID mode, every grid point is integrated for the same number
     * of its own Keplerian orbital periods, with the period recalculated from
     * the current semimajor-axis value of that grid point.
     *
     * @param mu_13 Gravitational parameter of the P1-P3 two-body problem
     *              [AU^3/day^2].
     * @param a Semimajor axis of the current P3 orbit [AU].
     *
     * @return Physical integration duration for the current orbit [day].
     *
     * @throws std::runtime_error If the integration-duration input method
     *         has not been specified.
     */
    double calcIntegrationDuration(double mu_13, double a) const;

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
     * @brief Specifies how the integration duration is defined.
     *
     * The integration duration can be specified either directly as a
     * physical time interval in days or as a number of initial Keplerian
     * orbital periods of the massless body P3.
     */
    enum class IntegrationDurationInput {
        NONE,           /**< Integration duration has not been specified. */
        PHYSICAL_TIME,  /**< Integration duration is specified by T [day]. */
        ORBITAL_PERIODS /**< Integration duration is specified by nPeriods. */
    };
    /**
     * @brief Returns the name of an orbital-phase input type.
     *
     * @param input Orbital-phase input type.
     * @return Name of the orbital-phase input type.
     */
    static const char *orbitalPhaseInputToString(OrbitalPhaseInput input) noexcept;

    /**
     * @brief Returns the name of an integration-duration input method.
     *
     * @param input Integration-duration input method.
     * @return Name of the integration-duration input method.
     */
    static const char *integrationDurationInputToString(IntegrationDurationInput input) noexcept;

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
     * @brief Checks whether an orbital element is controlled by a grid axis.
     *
     * @param element Orbital element to search for.
     *
     * @return True if the orbital element is present among the grid axes,
     *         false otherwise.
     */
    bool hasGridAxis(OrbitalElement element) const noexcept;

    /**
     * @brief Checks whether an orbital element was explicitly specified
     *        as a fixed input value.
     *
     * @param element Orbital element to search for.
     *
     * @return True if the orbital element was specified as a fixed input
     *         value, false otherwise.
     */
    bool hasFixedOrbitalElement(OrbitalElement element) const noexcept;

    /**
     * @brief Registers an orbital element as an explicitly specified
     *        fixed input value.
     *
     * @param element Orbital element to register.
     *
     * @throws std::runtime_error If the orbital element has already been
     *         specified as a fixed input value.
     */
    void registerFixedOrbitalElement(OrbitalElement element);

    /**
     * @brief Validates how an orbital element is specified in GRID mode.
     *
     * The orbital element must be specified exactly once: either as a fixed
     * input value or as a grid axis, but not both.
     *
     * @param element Orbital element to validate.
     *
     * @throws std::runtime_error If the orbital element is specified both
     *         as a fixed value and as a grid axis, or by neither method.
     */
    void validateOrbitalElementSource(OrbitalElement element) const;

    /**
     * @brief Removes whitespace characters from a string.
     *
     * @param text String to modify.
     */
    static void removeSpaces(std::string &text);

    /**
     * @brief Removes leading and trailing whitespace from a string.
     *
     * Internal whitespace characters are preserved.
     *
     * @param text String to modify.
     */
    static void Trim(std::string &text);

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
    /// Number of initial Keplerian orbital periods of P3.
    double n_periods_ = 0.0;
    /// Physical output time interval [day].
    double output_dt_ = 0.0;

    /// Method used to specify the integration duration.
    IntegrationDurationInput integration_duration_input_ = IntegrationDurationInput::NONE;
    /// Specifies how the initial orbital phase was provided.
    OrbitalPhaseInput      orbital_phase_input_ = OrbitalPhaseInput::NONE;
    astro::OrbitalElements elements_{};
    /**
     * @brief Orbital elements explicitly specified as fixed input values.
     *
     * An orbital element stored here was explicitly provided in the input
     * file rather than obtained from a grid axis.
     */
    std::set<OrbitalElement> fixed_elements_;

    /// Initial mean anomaly [rad].
    double mean_anomaly_ = 0.0;

    /**
     * @brief Definitions of the orbital-element grid axes.
     */
    std::vector<GridAxis> grid_axes_;

    double        a0_ = 0.0;
    double        a1_ = 0.0;
    std::uint32_t Na_ = 0;

    double        e0_ = 0.0;
    double        e1_ = 0.0;
    std::uint32_t Ne_ = 0;

    double dy_[4] = {1.0, 0.0, 0.0, 0.0};

    /// Relative tolerance of the numerical integrator.
    double rel_tol_ = 1.0e-6;

    /// Absolute tolerance of the numerical integrator.
    double abs_tol_ = 1.0e-10;
};