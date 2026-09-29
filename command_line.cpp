#include "command_line.h"

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
}  // namespace

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

            opt.input_dir     = value;
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

            opt.output_dir     = value;
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
                opt.input_path = fs::absolute(fs::path(opt.input_dir) / opt.input_file).string();
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
                opt.output_path = fs::absolute(fs::path(opt.output_dir) / opt.output_file).string();
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