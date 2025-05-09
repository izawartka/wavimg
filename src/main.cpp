#include "common.hpp"
#include "cli.hpp"
#include "file_audio.hpp"
#include "file_image.hpp"
#include "spectrogram.hpp"

int main(int argc, char** argv) {
    try {
        auto opts = CLI::parse(argc, argv);

        if (opts.encode) {
            FileAudio audio = FileAudio::load(opts.input);
            FileImage img;
            Spectrogram::encode(audio, img, opts.spectrogramOptions);

            img.save(opts.output);
        }
        else {
            FileImage img = FileImage::load(opts.input);
            FileAudio audio;
            Spectrogram::decode(img, audio, opts.spectrogramOptions);

            audio.save(opts.output);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Done!" << std::endl;

    return 0;
}
