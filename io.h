#pragma once

#include <filesystem>  // std::filesystem::path
#include <iosfwd>      // std::ostream
#include <string>      // std::string
#include <vector>      // std::vector

class InitData;  // Forward declaration used by the I/O interface.

namespace io {

    /**
     * @brief Width of each numerical output column.
     */
    inline constexpr int DATA_FIELD_WIDTH = 18;

    /**
     * @brief Number of digits written after the decimal point in scientific notation.
     */
    inline constexpr int DATA_PRECISION = 10;

    /**
     * @brief Describes one column of a numerical output table.
     *
     * Stores the machine-readable column name, its physical unit, and a
     * human-readable description. The same structure is used both for the
     * reproducibility header and for generating the actual table header.
     */
    struct OutputColumn {
        std::string name;         ///< Column name.
        std::string unit;         ///< Physical unit, or "-" if dimensionless.
        std::string description;  ///< Human-readable description.
    };

    /**
     * @brief Configures an output stream for numerical data.
     *
     * Numerical values are written in scientific notation, right-aligned,
     * with an explicit plus sign for positive values, using the common
     * output precision defined by DATA_PRECISION.
     *
     * The function centralizes the numerical formatting used by ORBIT,
     * INDICATOR, and GRID output.
     *
     * @param out Output stream.
     */
    void configureNumericalOutput(std::ostream &out);

    /**
     * @brief Writes the complete reproducibility header.
     *
     * The header contains:
     *  - program and run metadata;
     *  - an exact commented copy of the input file;
     *  - a detailed description of the numerical output columns.
     *
     * The same header is written regardless of whether the numerical output
     * is directed to a file or to the standard output stream.
     *
     * @param out Output stream.
     * @param inputPath Path of the initialization file used for the run.
     * @param init Initialization data.
     */
    void writeOutputHeader(std::ostream &out, const std::filesystem::path &inputPath, const InitData &init);

    /**
     * @brief Writes a commented copy of the original input file.
     *
     * Each original input line is written unchanged after the "# | " prefix,
     * so that the complete input configuration is embedded in the numerical
     * output while remaining a comment for data-processing tools.
     *
     * @param out Output stream.
     * @param inputPath Path of the original initialization file.
     *
     * @throws std::runtime_error If the input file cannot be opened.
     */
    void writeInputFileCopy(std::ostream &out, const std::filesystem::path &inputPath);

    /**
     * @brief Writes program and run metadata to the output stream.
     *
     * The metadata identifies the software version, build environment,
     * numerical integration method, computation mode, chaos indicator,
     * mathematical formalism, and input file associated with the run.
     *
     * This information forms part of the reproducibility header written
     * before the numerical output data.
     *
     * @param out Output stream.
     * @param inputPath Path of the initialization file used for the run.
     * @param init Initialization data.
     */
    void writeProgramInformation(std::ostream &out, const std::filesystem::path &inputPath, const InitData &init);

    /**
     * @brief Writes a detailed description of the numerical output columns.
     *
     * The structure is generated from the same OutputColumn definitions that
     * are later used for the actual table header, ensuring consistency between
     * the documented and written data layout.
     *
     * @param out Output stream.
     * @param init Initialization data.
     */
    void writeOutputStructure(std::ostream &out, const InitData &init);

    /**
     * @brief Writes the numerical table header.
     *
     * The column labels are generated from the same output schema that is used
     * in the reproducibility header. The common output column width defined by
     * DATA_FIELD_WIDTH is used for every column.
     *
     * @param out Output stream.
     * @param init Initialization data.
     */
    void writeTableHeader(std::ostream &out, const InitData &init);

    /**
     * @brief Builds the description of the numerical output columns.
     *
     * The returned schema depends on the selected run mode, chaos indicator,
     * mathematical formalism, and, in GRID mode, the selected grid axes.
     *
     * @param init Initialization data.
     *
     * @return Ordered list of output-column definitions.
     */
    std::vector<OutputColumn> buildOutputSchema(const InitData &init);

}  // namespace io