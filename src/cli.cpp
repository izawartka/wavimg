#include "cli.hpp"

CLIOptions CLI::parse(int argc, char** argv) {
    cxxopts::Options options("wavimg", "Convert WAV <-> PNG Spectrogram");
    options.positional_help("<input> <output>");
    options.add_options()
        ("c,encode", "Convert WAV to PNG Spectrogram", cxxopts::value<bool>())
        ("d,decode", "Convert PNG Spectrogram to WAV", cxxopts::value<bool>())
        ("n_fft", "FFT size", cxxopts::value<size_t>()->default_value("2048"))
        ("hop_length", "Hop length", cxxopts::value<size_t>()->default_value("512"))
        ("win_length", "Window length (0 - same as n_fft)", cxxopts::value<size_t>()->default_value("0"))
        ("min_db", "Minimum dB", cxxopts::value<float>()->default_value("-80.0"))
        ("max_db", "Maximum dB", cxxopts::value<float>()->default_value("0.0"))
        ("gi", "Griffin-Lim iterations", cxxopts::value<size_t>()->default_value("60"))
        ("sr", "Sample rate", cxxopts::value<size_t>()->default_value("44100"))
        ("h,help", "Print help")
        ("input", "Input file", cxxopts::value<std::string>())
        ("output", "Output file", cxxopts::value<std::string>());

    options.parse_positional({ "input", "output" });

    auto result = options.parse(argc, argv);

    bool show_help = result.count("help") > 0;
    bool valid_mode = result.count("encode") ^ result.count("decode"); // XOR
    bool has_input_output = result.count("input") > 0 && result.count("output") > 0;

    if (show_help || !valid_mode || !has_input_output) {
        std::cout << options.help() << std::endl;
        exit(1);
    }

    CLIOptions opts;
    opts.encode = result["encode"].as<bool>();
    opts.decode = result["decode"].as<bool>();
    opts.input = result["input"].as<std::string>();
    opts.output = result["output"].as<std::string>();

    opts.spectrogramOptions.n_fft = result["n_fft"].as<size_t>();
    opts.spectrogramOptions.hop_length = result["hop_length"].as<size_t>();
    opts.spectrogramOptions.win_length = result["win_length"].as<size_t>();
    opts.spectrogramOptions.min_db = result["min_db"].as<float>();
    opts.spectrogramOptions.max_db = result["max_db"].as<float>();
    opts.spectrogramOptions.griffin_iter = result["gi"].as<size_t>();
    opts.spectrogramOptions.sample_rate = result["sr"].as<size_t>();
    if (opts.spectrogramOptions.win_length == 0) {
        opts.spectrogramOptions.win_length = opts.spectrogramOptions.n_fft;
    }

    return opts;
}
