#pragma once
#include <cxxopts.hpp>
#include "common.hpp"
#include "spec_options.hpp"

struct CLIOptions {
    bool encode = false;
    bool decode = false;
    std::string input;
    std::string output;
    SpectrogramOptions spectrogramOptions;
};

class CLI {
public:
    static CLIOptions parse(int argc, char** argv);
};
