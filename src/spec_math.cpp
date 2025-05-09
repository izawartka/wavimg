#include "spec_math.hpp"

void SpecMath::fft(
    std::vector<double>& real, 
    std::vector<double>& imag, 
    bool inverse
) {
    const size_t n = real.size();
    if (n != imag.size()) {
        throw std::invalid_argument("FFT error: real and imag parts must be the same length");
    }

    // Check power of two
    if ((n & (n - 1)) != 0) {
        throw std::invalid_argument("FFT error: size must be a power of two (got " + std::to_string(n) + ")");
    }

    // Bit-reversal permutation
    unsigned int j = 0;
    for (unsigned int i = 1; i < n; ++i) {
        unsigned int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j |= bit;
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
    }

    // FFT
    for (size_t len = 2; len <= n; len <<= 1) {
        double angle = 2 * M_PI / len * (inverse ? -1 : 1);
        double wlen_real = std::cos(angle);
        double wlen_imag = std::sin(angle);
        for (size_t i = 0; i < n; i += len) {
            double w_real = 1.0;
            double w_imag = 0.0;
            for (size_t k = 0; k < len / 2; ++k) {
                size_t u = i + k;
                size_t v = i + k + len / 2;
                double u_real = real[u];
                double u_imag = imag[u];
                double v_real = real[v] * w_real - imag[v] * w_imag;
                double v_imag = real[v] * w_imag + imag[v] * w_real;

                real[u] = u_real + v_real;
                imag[u] = u_imag + v_imag;
                real[v] = u_real - v_real;
                imag[v] = u_imag - v_imag;

                double next_w_real = w_real * wlen_real - w_imag * wlen_imag;
                double next_w_imag = w_real * wlen_imag + w_imag * wlen_real;
                w_real = next_w_real;
                w_imag = next_w_imag;
            }
        }
    }
    if (inverse) {
        for (size_t i = 0; i < n; ++i) {
            real[i] /= n;
            imag[i] /= n;
        }
    }
}

SpectrogramMatrix SpecMath::stft(
    const std::vector<double>& samples,
    size_t n_fft,
    size_t hop_length,
    size_t win_length
) {
    auto window = hannWindow(win_length);
    size_t n_frames = (samples.size() - win_length) / hop_length + 1;
    SpectrogramMatrix matrix(n_fft / 2 + 1, std::vector<double>(n_frames));
    for (size_t frame = 0; frame < n_frames; ++frame) {
        std::vector<double> real(n_fft, 0.0), imag(n_fft, 0.0);
        for (size_t n = 0; n < win_length; ++n) {
            size_t idx = frame * hop_length + n;
            if (idx < samples.size())
                real[n] = samples[idx] * window[n];
            else
                real[n] = 0.0;
        }
        SpecMath::fft(real, imag, false);
        for (size_t k = 0; k <= n_fft / 2; ++k)
            matrix[k][frame] = std::hypot(real[k], imag[k]);
    }
    return matrix;
}

ComplexSpectrogramMatrix SpecMath::stftComplex(
    const std::vector<double>& samples,
    size_t n_fft,
    size_t hop_length,
    size_t win_length
) {
    auto window = hannWindow(win_length);
    size_t n_frames = (samples.size() - win_length) / hop_length + 1;
    SpectrogramMatrix real_spec(n_fft / 2 + 1, std::vector<double>(n_frames));
    SpectrogramMatrix imag_spec(n_fft / 2 + 1, std::vector<double>(n_frames));

    for (size_t frame = 0; frame < n_frames; ++frame) {
        std::vector<double> real(n_fft, 0.0), imag(n_fft, 0.0);
        for (size_t n = 0; n < win_length; ++n) {
            size_t idx = frame * hop_length + n;
            if (idx < samples.size())
                real[n] = samples[idx] * window[n];
            else
                real[n] = 0.0;
        }
        SpecMath::fft(real, imag, false);
        for (size_t k = 0; k <= n_fft / 2; ++k) {
            real_spec[k][frame] = real[k];
            imag_spec[k][frame] = imag[k];
        }
    }
    return { real_spec, imag_spec };
}

std::vector<double> SpecMath::istft(
    const SpectrogramMatrix& spec,
    size_t hop_length,
    size_t win_length,
    size_t length
) {
    size_t n_fft = (spec.size() - 1) * 2;
    auto window = hannWindow(win_length);
    std::vector<double> signal(length, 0.0);
    std::vector<double> sum_window(length, 0.0);
    size_t n_frames = spec[0].size();
    for (size_t frame = 0; frame < n_frames; ++frame) {
        std::vector<double> real(n_fft), imag(n_fft);
        // create full spectrum with phase = 0
        for (size_t k = 0; k <= n_fft / 2; ++k) {
            real[k] = spec[k][frame];
            imag[k] = 0;

            if (k != 0 && k != n_fft / 2) {
                real[n_fft - k] = spec[k][frame];
                imag[n_fft - k] = 0;
            }
        }
        SpecMath::fft(real, imag, true);
        for (size_t n = 0; n < win_length; ++n) {
            signal[frame * hop_length + n] += real[n] * window[n];
            sum_window[frame * hop_length + n] += window[n] * window[n];
        }
    }
    // normalize by window overlap
    for (size_t i = 0; i < length; ++i)
        if (sum_window[i] > 0)
            signal[i] /= sum_window[i];

    return signal;
}

std::vector<double> SpecMath::istftComplex(
    const ComplexSpectrogramMatrix& spec,
    size_t hop_length,
    size_t win_length,
    size_t length
) {
    size_t n_fft = (spec.first.size() - 1) * 2;
    size_t n_frames = spec.first[0].size();
    auto window = hannWindow(win_length);
    std::vector<double> y(length, 0.0);
    std::vector<double> window_sum(length, 0.0);

    std::vector<double> real(n_fft), imag(n_fft);
    std::vector<double> frame(n_fft);

    for (size_t frame_idx = 0; frame_idx < n_frames; ++frame_idx) {
        for (size_t k = 0; k <= n_fft / 2; ++k) {
            real[k] = spec.first[k][frame_idx];
            imag[k] = spec.second[k][frame_idx];

            if (k != 0 && k != n_fft / 2) {
                real[n_fft - k] = spec.first[k][frame_idx];
                imag[n_fft - k] = -spec.second[k][frame_idx];
            }
        }

        SpecMath::fft(real, imag, true); // Inverse FFT

        size_t start = frame_idx * hop_length;
        for (size_t n = 0; n < win_length; ++n) {
            if (start + n < length) {
                y[start + n] += real[n] * window[n];
                window_sum[start + n] += window[n] * window[n];
            }
        }
    }

    // Normalize to compensate for windowing
    for (size_t i = 0; i < length; ++i)
        if (window_sum[i] > 1e-8)
            y[i] /= window_sum[i];

    return y;
}

std::vector<double> SpecMath::hannWindow(size_t win_length) {
    std::vector<double> win(win_length);
    for (size_t n = 0; n < win_length; ++n)
        win[n] = 0.5 * (1 - std::cos(2 * M_PI * n / (win_length - 1)));
    return win;
}

ComplexSpectrogramMatrix SpecMath::griffinLim(
    const SpectrogramMatrix& spec,
    size_t n_fft,
    size_t hop_length,
    size_t win_length,
    size_t iterations
) {
    size_t n_frames = spec[0].size();

    // Initialize random phase
    SpectrogramMatrix angle(spec.size(), std::vector<double>(n_frames));
    for (size_t k = 0; k < spec.size(); ++k)
        for (size_t t = 0; t < n_frames; ++t)
            angle[k][t] = 2 * M_PI * ((double)rand() / RAND_MAX);

    SpectrogramMatrix real_spec(spec.size(), std::vector<double>(n_frames));
    SpectrogramMatrix ispec_spec(spec.size(), std::vector<double>(n_frames));

    for (size_t iter = 0; iter < iterations; ++iter) {
        // Construct complex spectrogram with current phase
        for (size_t k = 0; k < spec.size(); ++k)
            for (size_t t = 0; t < n_frames; ++t) {
                real_spec[k][t] = spec[k][t] * std::cos(angle[k][t]);
                ispec_spec[k][t] = spec[k][t] * std::sin(angle[k][t]);
            }

        // Inverse STFT
        size_t length = (n_frames - 1) * hop_length + win_length;
        auto y = SpecMath::istftComplex({ real_spec, ispec_spec }, hop_length, win_length, length);

        // STFT
        auto [stft_real, stft_imag] = SpecMath::stftComplex(y, n_fft, hop_length, win_length);

        // Update phase
        for (size_t k = 0; k < spec.size(); ++k)
            for (size_t t = 0; t < n_frames; ++t)
                angle[k][t] = std::atan2(stft_imag[k][t], stft_real[k][t]);
    }

    // Final complex spectrogram
    for (size_t k = 0; k < spec.size(); ++k)
        for (size_t t = 0; t < n_frames; ++t) {
            real_spec[k][t] = spec[k][t] * std::cos(angle[k][t]);
            ispec_spec[k][t] = spec[k][t] * std::sin(angle[k][t]);
        }

    return { real_spec, ispec_spec };
}