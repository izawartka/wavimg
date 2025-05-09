#include "common.hpp"
#include "spec_types.hpp"
#include "spec_options.hpp"
#include "file_image.hpp"
#include "file_audio.hpp"

class Spectrogram {
public:
    static void encode(
        FileAudio& audio,
        FileImage& img,
        const SpectrogramOptions& opt
    );

    static void decode(
        FileImage& img,
        FileAudio& audio,
        const SpectrogramOptions& opt
    );

private:
    static inline size_t getWinLength(const SpectrogramOptions& opt);

    static void stftMatrixToDB(
        SpectrogramMatrix& stft_matrix,
        float min_db,
        float max_db
    );

    static void normalizeSamples(
        std::vector<double>& samples
    );

    static void validateImageSize(
        const FileImage& img,
        size_t n_fft,
        size_t hop_length
    );

    static void saveToImage(
        const SpectrogramMatrix& mag_db,
        SpectrogramOptions opt,
        FileImage& img_out
    );

    static SpectrogramMatrix loadFromImage(
        FileImage& img,
        float min_db,
        float max_db
    );

    static std::vector<double> addPadding(
        const std::vector<double>& samples,
        size_t n_fft
    );

    static std::vector<double> removePadding(
        const std::vector<double>& samples,
        size_t n_fft
    );
};
