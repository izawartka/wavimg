#include "spectrogram.hpp"
#include "spec_math.hpp"

void Spectrogram::encode(FileAudio& audio, FileImage& img, const SpectrogramOptions& opt)
{
    auto samples = addPadding(audio.samples, opt.n_fft);
    auto out_matrix = SpecMath::stft(samples, opt.n_fft, opt.hop_length, getWinLength(opt));
    stftMatrixToDB(out_matrix, opt.min_db, opt.max_db);

    saveToImage(out_matrix, opt, img);
}

void Spectrogram::decode(FileImage& img, FileAudio& audio, const SpectrogramOptions& opt)
{
    validateImageSize(img, opt.n_fft, opt.hop_length);
    size_t win_length = getWinLength(opt);

    SpectrogramMatrix matrix = loadFromImage(img, opt.min_db, opt.max_db);
    auto [recon_real, recon_imag] = SpecMath::griffinLim(matrix, opt.n_fft, opt.hop_length, win_length, opt.griffin_iter);
    size_t length = (recon_real[0].size() - 1) * opt.hop_length + win_length;
    auto samples = SpecMath::istftComplex({ recon_real, recon_imag }, opt.hop_length, win_length, length);
    samples = removePadding(samples, opt.n_fft);
    normalizeSamples(samples);

    audio.sampleRate = opt.sample_rate;
    audio.samples = std::move(samples);
}

inline size_t Spectrogram::getWinLength(const SpectrogramOptions& opt)
{
    return opt.win_length == 0 || opt.win_length > opt.n_fft ? opt.n_fft : opt.win_length;
}

void Spectrogram::stftMatrixToDB(SpectrogramMatrix& matrix, float min_db, float max_db)
{
    double ref = 0.0;
    for (auto& row : matrix) {
        for (double& v : row) {
            v = std::abs(v);
            ref = std::max(ref, v);
        }
    }

    for (auto& row : matrix) {
        for (auto& v : row) {
            v = 20 * std::log10(std::max(v, 1e-10) / ref);
            v = std::min(std::max(v, (double)min_db), (double)max_db);
        }
    }
}

void Spectrogram::normalizeSamples(std::vector<double>& samples)
{
    double max_val = 0.0;
    for (double v : samples) max_val = std::max(max_val, std::abs(v));
    for (auto& v : samples) v = v / max_val * 0.95;
}

void Spectrogram::validateImageSize(const FileImage& img, size_t n_fft, size_t hop_length)
{
    if (img.height != n_fft / 2 + 1) {
        size_t calculated_n_fft = (img.height - 1) * 2;
        throw std::invalid_argument("Image size does not match the expected size for the given FFT parameters. Try using --n_fft=" + std::to_string(calculated_n_fft) + ".");
    }
}

void Spectrogram::saveToImage(
    const SpectrogramMatrix& matrix,
    SpectrogramOptions opt,
    FileImage& img_out
) {
    size_t w = matrix[0].size(), h = matrix.size();
    img_out.data.resize(w * h * 4);
    double range = opt.max_db - opt.min_db;
    for (size_t x = 0; x < w; ++x)
        for (size_t y = 0; y < h; ++y) {
            float norm = (matrix[h - 1 - y][x] - opt.min_db) / range;
            uint8_t val = static_cast<uint8_t>(norm * 255);
            size_t idx = (y * w + x) * 4;
            img_out.data[idx] = val;
            img_out.data[idx + 1] = val;
            img_out.data[idx + 2] = val;
            img_out.data[idx + 3] = 255;
        }
    img_out.width = w;
    img_out.height = h;
}

SpectrogramMatrix Spectrogram::loadFromImage(
    FileImage& img,
    float min_db,
    float max_db
) {
    size_t w = img.width, h = img.height;
    SpectrogramMatrix matrix(h, std::vector<double>(w));
    double range = max_db - min_db;
    for (size_t x = 0; x < w; ++x)
        for (size_t y = 0; y < h; ++y) {
            uint8_t byte_val = img.data[(y * w + x) * 4];
            double db_val = min_db + (byte_val / 255.0) * range;
            double mag_val = std::pow(10.0, db_val / 20.0);
            matrix[h - 1 - y][x] = mag_val;
        }

    return matrix;
}

std::vector<double> Spectrogram::addPadding(const std::vector<double>& samples, size_t n_fft)
{
    std::vector<double> padded_samples(n_fft + samples.size(), 0.0);
    std::copy(samples.begin(), samples.end(), padded_samples.begin() + n_fft / 2);

    return padded_samples;
}

std::vector<double> Spectrogram::removePadding(const std::vector<double>& samples, size_t n_fft)
{
    size_t start = n_fft / 2;
    size_t end = samples.size() - n_fft / 2;

    return std::vector<double>(samples.begin() + start, samples.begin() + end);
}
