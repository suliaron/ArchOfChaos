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
    void Print(std::ostream &os = std::cout) const;

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

    /** @return Initial integration time. */
    double getT0() const noexcept
    {
        return t0_;
    }

    /** @return Final integration time. */
    double getT() const noexcept
    {
        return T_;
    }

    /** @return Output time interval. */
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

   private:
    /**
     * @brief Parses a single line of the initialization file.
     *
     * Empty lines and comments are ignored.
     *
     * @param line Input line.
     */
    void ParseLine(const std::string &line);

    /**
     * @brief Validates the initialization data.
     *
     * @throws std::runtime_error If the stored input parameters are
     *         inconsistent or invalid.
     */
    void Validate() const;

    /**
     * @brief Removes whitespace characters from a string.
     *
     * @param text String to modify.
     */
    static void RemoveSpaces(std::string &text);

    /**
     * @brief Converts text to a run mode.
     *
     * @param text Run mode as text.
     * @return Corresponding run mode.
     *
     * @throws std::runtime_error If the value is unknown.
     */
    static RunMode ParseRunMode(const std::string &text);

    /**
     * @brief Converts text to an indicator type.
     *
     * @param text Indicator type as text.
     * @return Corresponding indicator type.
     *
     * @throws std::runtime_error If the value is unknown.
     */
    static Model::IndicatorType ParseIndicatorType(const std::string &text);

    /**
     * @brief Converts text to a CRTBP mathematical formulation.
     *
     * @param text Formalism as text.
     * @return Corresponding CRTBP formulation.
     *
     * @throws std::runtime_error If the value is unknown.
     */
    static Model::Formalism ParseFormalism(const std::string &text);

    RunMode              run_mode_  = RunMode::ORBIT;
    Model::Formalism     formalism_ = Model::Formalism::NEWTONIAN;
    Model::IndicatorType indicator_ = Model::IndicatorType::NONE;

    double m1_ = 0.0;
    double m2_ = 0.0;
    double a2_ = 0.0;

    double t0_        = 0.0;
    double T_         = 0.0;
    double output_dt_ = 0.0;

    astro::OrbitalElements elements_{};

    double        a0_ = 0.0;
    double        a1_ = 0.0;
    std::uint32_t Na_ = 0;

    double        e0_ = 0.0;
    double        e1_ = 0.0;
    std::uint32_t Ne_ = 0;

    double dy_[4] = {1.0, 0.0, 0.0, 0.0};
};