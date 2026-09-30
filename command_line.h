#pragma once

#include <iostream>
#include <ostream>
#include <string>

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
    void print(std::ostream& os = std::cout) const;
};

/**
 * @brief Prints program version information.
 */
void printVersion();

/**
 * @brief Prints the command-line help.
 */
void printHelp();

/**
 * @brief Parses command-line arguments.
 *
 * Processes command-line options related to input and output files and
 * directories, as well as general program options.
 *
 * Input and output files can be specified in two ways:
 *
 * @code
 * -i <file>
 * -o <file>
 * @endcode
 *
 * where <file> may contain either a file name only or a complete or relative
 * path, or by specifying the directory separately:
 *
 * @code
 * -iDir <directory> -i <file>
 * -oDir <directory> -o <file>
 * @endcode
 *
 * If -iDir is specified, the argument of -i must contain a file name only and
 * must not contain a path separator. The same restriction applies to -o when
 * -oDir is specified.
 *
 * Directory options are processed in a first pass so that the order of
 * -iDir/-oDir relative to -i/-o does not matter.
 *
 * The resulting input and output paths are converted to absolute paths.
 *
 * Supported options are:
 *
 * @code
 * -i <file>            Input file.
 * -iDir <directory>    Input directory.
 * -idir <directory>    Alias for -iDir.
 * -o <file>            Output file.
 * -oDir <directory>    Output directory.
 * -odir <directory>    Alias for -oDir.
 * -h, --help           Display help information.
 * -v, --version        Display version information.
 * --verbose            Enable verbose output.
 * @endcode
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param opt Structure receiving the parsed command-line options.
 *
 * @throws std::runtime_error If an option requiring an argument is missing
 *         its argument, if both a directory option and a path-containing
 *         file argument are specified, or if an unknown command-line option
 *         is encountered.
 */
void parseCommandLine(int argc, char* argv[], CommandLineOptions& opt);
