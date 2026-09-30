#include "command_line.h"
#include "version.h"

#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {
    /**
     * @brief Checks whether a string contains a path separator.
     *
     * On Windows, both '/' and '\\' are treated as path separators.
     * On other platforms, only '/' is considered a path separator.
     *
     * @param str The string to examine.
     * @return true if the string contains a path separator, false otherwise.
     */
    bool containsPathSeparator(const std::string &str)
    {
#ifdef _WIN32
        return str.find_first_of("/\\") != std::string::npos;
#else
        return str.find('/') != std::string::npos;
#endif
    }

    /**
     * @brief Checks whether a string begins or ends with whitespace.
     *
     * @param str String to examine.
     * @return true if the first or last character is whitespace, otherwise false.
     */
    bool hasLeadingOrTrailingWhitespace(const std::string &str)
    {
        if (str.empty()) {
            return false;
        }

        return std::isspace(static_cast<unsigned char>(str.front())) ||
               std::isspace(static_cast<unsigned char>(str.back()));
    }

}  // namespace

void printVersion()
{
    std::cout << program::name << '\n'
              << "Version     : " << program::version << '\n'
              << "Author      : " << program::author << '\n'
              << "Affiliation : " << program::affiliation << '\n'
              << "Description : " << program::description << '\n'
              << program::copyright << '\n';
}

/**
 * @brief Prints the command-line help.
 *
 * Displays the program name and version, command-line syntax,
 * supported options, and usage examples.
 */
void printHelp()
{
    std::cout << program::name << " " << program::version << "\n";
    std::cout << "========================================\n\n";
    std::cout << program::description << "\n\n";
    std::cout << "Usage:\n";
    std::cout << "  " << program::executable << " -i <input file> -o <output file>\n";
    std::cout << "  " << program::executable << " -iDir <input directory> -i <input file>"
              << " -oDir <output directory> -o <output file>\n\n";
    std::cout << "Options:\n";
    std::cout << "  -i <file>          Input file. The argument may contain a path unless\n";
    std::cout << "                     -iDir is specified.\n";
    std::cout << "  -iDir <directory>  Input directory.\n";
    std::cout << "  -idir <directory>  Alias for -iDir.\n";
    std::cout << "  -o <file>          Output file. The argument may contain a path unless\n";
    std::cout << "                     -oDir is specified.\n";
    std::cout << "  -oDir <directory>  Output directory.\n";
    std::cout << "  -odir <directory>  Alias for -oDir.\n";
    std::cout << "  -h, --help         Display this help message.\n";
    std::cout << "  -v, --version      Display program version information.\n";
    std::cout << "  --verbose          Display detailed input information.\n\n";

    std::cout << "Notes:\n";
    std::cout << "  If -iDir is specified, -i must contain a file name only.\n";
    std::cout << "  If -oDir is specified, -o must contain a file name only.\n";
    std::cout << "  Paths containing spaces must be enclosed in quotation marks.\n";
    std::cout << "  A quoted path must not end with a backslash ('\\').\n\n";

    std::cout << "Examples:\n";
    std::cout << "  " << program::executable << " -i input.txt -o output.txt\n";
    std::cout << "  " << program::executable << " -i data/input.txt -o results/output.txt\n";
    std::cout << "  " << program::executable << " -iDir data -i input.txt"
              << " -oDir results -o output.txt\n";
    std::cout << "  " << program::executable << " -iDir \"D:\\OneDrive - elte.hu\\Work\\Input\""
              << " -i input.txt"
              << " -oDir \"D:\\OneDrive - elte.hu\\Work\\Output\""
              << " -o output.txt\n";
}

void CommandLineOptions::print(std::ostream &os) const
{
    os << "----------------------------------------\n";
    os << "Command-line options\n";
    os << "----------------------------------------\n";

    os << "Input file    : " << input_file << '\n'
       << "Input dir     : " << input_dir << '\n'
       << "Input path    : " << input_path << '\n'
       << "Output file   : " << output_file << '\n'
       << "Output dir    : " << output_dir << '\n'
       << "Output path   : " << output_path << '\n'
       << "Show help     : " << std::boolalpha << show_help << '\n'
       << "Show version  : " << std::boolalpha << show_version << '\n'
       << "Verbose       : " << std::boolalpha << verbose << '\n';
}

void parseCommandLine(int argc, char *argv[], CommandLineOptions &opt)
{
    bool inputDirSpecified  = false;
    bool outputDirSpecified = false;

    //
    // First pass: process directory options.
    //
    for (int i = 1; i < argc; ++i) {
        const std::string key(argv[i]);

        if (key == "-iDir" || key == "-idir") {
            if (++i >= argc) {
                throw std::runtime_error("Missing argument after '-iDir'.");
            }

            const std::string value(argv[i]);

            if (!value.empty() && value[0] == '-') {
                throw std::runtime_error("Missing argument after '-iDir'.");
            }

            if (hasLeadingOrTrailingWhitespace(value)) {
                throw std::runtime_error("Input directory must not begin or end with whitespace: '" + value + "'");
            }

            const fs::path p(value);

            opt.input_dir     = fs::absolute(p).string();
            inputDirSpecified = true;
        }
        else if (key == "-oDir" || key == "-odir") {
            if (++i >= argc) {
                throw std::runtime_error("Missing argument after '-oDir'.");
            }

            const std::string value(argv[i]);

            if (!value.empty() && value[0] == '-') {
                throw std::runtime_error("Missing argument after '-oDir'.");
            }

            if (hasLeadingOrTrailingWhitespace(value)) {
                throw std::runtime_error("Output directory must not begin or end with whitespace: '" + value + "'");
            }

            const fs::path p(value);

            opt.output_dir     = fs::absolute(p).string();
            outputDirSpecified = true;
        }
    }

    //
    // Second pass: process file and general options.
    //
    for (int i = 1; i < argc; ++i) {
        const std::string key(argv[i]);

        //
        // Directory options were already processed in the first pass.
        //
        if (key == "-iDir" || key == "-idir" || key == "-oDir" || key == "-odir") {
            ++i;
            continue;
        }

        if (key == "-i") {
            if (++i >= argc) {
                throw std::runtime_error("Missing argument after '-i'.");
            }

            const std::string value(argv[i]);

            if (!value.empty() && value[0] == '-') {
                throw std::runtime_error("Missing argument after '-i'.");
            }

            if (inputDirSpecified) {
                if (containsPathSeparator(value)) {
                    throw std::runtime_error("-i and -iDir cannot contain directory at the same time.");
                }
                opt.input_file = value;
                opt.input_path = (fs::path(opt.input_dir) / opt.input_file).string();
            }
            else {
                const fs::path p = fs::absolute(fs::path(value));

                opt.input_file = p.filename().string();
                opt.input_dir  = p.parent_path().string();
                opt.input_path = p.string();
            }
        }
        else if (key == "-o") {
            if (++i >= argc) {
                throw std::runtime_error("Missing argument after '-o'.");
            }

            const std::string value(argv[i]);

            if (!value.empty() && value[0] == '-') {
                throw std::runtime_error("Missing argument after '-o'.");
            }

            if (outputDirSpecified) {
                if (containsPathSeparator(value)) {
                    throw std::runtime_error("-o and -oDir cannot contain directory at the same time.");
                }

                opt.output_file = value;
                opt.output_path = (fs::path(opt.output_dir) / opt.output_file).string();
            }
            else {
                const fs::path p = fs::absolute(fs::path(value));

                opt.output_file = p.filename().string();
                opt.output_dir  = p.parent_path().string();
                opt.output_path = p.string();
            }
        }
        else if (key == "-h" || key == "--help") {
            opt.show_help = true;
        }
        else if (key == "-v" || key == "--version") {
            opt.show_version = true;
        }
        else if (key == "--verbose") {
            opt.verbose = true;
        }
        else {
            throw std::runtime_error("Unknown command-line option: " + key);
        }
    }
}