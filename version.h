#pragma once

/**
 * @file version.h
 * @brief Program identification and version information.
 *
 * Contains compile-time constants describing the Arch of Chaos program,
 * including its name, version number, author, affiliation, and copyright
 * information.
 */

namespace program {

    /**
     * @brief Full program name.
     */
    inline constexpr char name[] = "Arch of Chaos";

    /**
     * @brief Command-line executable name.
     */
    inline constexpr char executable[] = "archofchaos";

    /**
     * @brief Major version number.
     */
    inline constexpr int versionMajor = 1;

    /**
     * @brief Minor version number.
     */
    inline constexpr int versionMinor = 5;

    /**
     * @brief Patch version number.
     */
    inline constexpr int versionPatch = 0;

    /**
     * @brief Complete program version string.
     *
     * The version follows the MAJOR.MINOR.PATCH convention.
     */
    inline constexpr char version[] = "1.5.0";

    /**
     * @brief Program author.
     */
    inline constexpr char author[] = "Dr. Áron Süli";

    /**
     * @brief Author affiliation.
     */
    inline constexpr char affiliation[] = "Eötvös Loránd University (ELTE)";

    /**
     * @brief Short description of the program.
     */
    inline constexpr char description[] =
        "Numerical investigation of dynamical structures "
        "in the circular restricted three-body problem.";

    /**
     * @brief Copyright notice.
     */
    inline constexpr char copyright[] = "Copyright (c) 2026 Dr. Áron Süli";

}  // namespace program