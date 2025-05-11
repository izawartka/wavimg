#include "spec_math.hpp"

void SpecMath::fft(
    std::vector<double>& real,
    std::vector<double>& imag,
    bool inverse
) {
    const size_t n = real.size();
    if (n != imag.size() || (n & (n - 1)))
        throw std::invalid_argument("FFT size must be power-of-two and real/imag lengths match");

    // Bit-reversal permutation
    unsigned j = 0;
    for (unsigned i = 1; i < n; ++i) {
        unsigned bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j |= bit;
        if (i < j) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
    }

    const auto& [W_real, W_imag] = twiddles(n, inverse);
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len >> 1;
        size_t stride = n / len;
        for (size_t i = 0; i < n; i += len) {
            for (size_t k = 0; k < half; ++k) {
				size_t w_idx = k * stride;
				double w_r = W_real[w_idx];
				double w_i = W_imag[w_idx];
                size_t u = i + k;
                size_t v = u + half;
                double u_r = real[u];
                double u_i = imag[u];
				double v_r = real[v] * w_r - imag[v] * w_i;
				double v_i = real[v] * w_i + imag[v] * w_r;
                real[u] = u_r + v_r;
                imag[u] = u_i + v_i;
                real[v] = u_r - v_r;
                imag[v] = u_i - v_i;
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
            if (idx < samples.size()) {
                real[n] = samples[idx] * window[n];
            } else {
                real[n] = 0.0;
            }
        }
        SpecMath::fft(real, imag, false);
        for (size_t k = 0; k <= n_fft / 2; ++k) {
            matrix[k][frame] = std::hypot(real[k], imag[k]);
        }
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
            if (idx < samples.size()) {
                real[n] = samples[idx] * window[n];
            } else {
                real[n] = 0.0;
            }
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
    for (size_t i = 0; i < length; ++i) {
        if (sum_window[i] > 0) {
            signal[i] /= sum_window[i];
        }
    }

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
    for (size_t i = 0; i < length; ++i) {
        if (window_sum[i] > 1e-8) {
            y[i] /= window_sum[i];
        }
    }

    return y;
}

std::map<size_t, std::vector<double>> SpecMath::m_hannWindowCache;

std::vector<double> SpecMath::hannWindow(size_t win_length) {
	if (m_hannWindowCache.find(win_length) != m_hannWindowCache.end()) {
		return m_hannWindowCache[win_length];
	}

	std::vector<double> window(win_length);
	for (size_t n = 0; n < win_length; ++n) {
		window[n] = 0.5 * (1 - std::cos(2 * M_PI * n / (win_length - 1)));
	}
	m_hannWindowCache[win_length] = window;
	return window;
}

std::map<size_t, std::pair<std::vector<double>, std::vector<double>>> SpecMath::m_twiddlesCache;

std::pair<std::vector<double>, std::vector<double>> SpecMath::twiddles(size_t n_fft, bool inverse)
{
	if (m_twiddlesCache.find(n_fft) != m_twiddlesCache.end()) {
		return m_twiddlesCache[n_fft];
	}

	std::vector<double> W_real(n_fft), W_imag(n_fft);
	for (size_t k = 0; k < n_fft; ++k) {
		double angle = 2 * M_PI * k / n_fft;
		if (inverse) angle = -angle;
		W_real[k] = std::cos(angle);
		W_imag[k] = std::sin(angle);
	}
	m_twiddlesCache[n_fft] = { W_real, W_imag };
	return { W_real, W_imag };
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
    for (size_t k = 0; k < spec.size(); ++k) {
        for (size_t t = 0; t < n_frames; ++t) {
            angle[k][t] = 2 * M_PI * ((double)rand() / RAND_MAX);
        }
    }

    SpectrogramMatrix real_spec(spec.size(), std::vector<double>(n_frames));
    SpectrogramMatrix ispec_spec(spec.size(), std::vector<double>(n_frames));

    for (size_t iter = 0; iter < iterations; ++iter) {
        // Construct complex spectrogram with current phase
        for (size_t k = 0; k < spec.size(); ++k) {
            for (size_t t = 0; t < n_frames; ++t) {
                real_spec[k][t] = spec[k][t] * std::cos(angle[k][t]);
                ispec_spec[k][t] = spec[k][t] * std::sin(angle[k][t]);
            }
        }

        // Inverse STFT
        size_t length = (n_frames - 1) * hop_length + win_length;
        auto y = SpecMath::istftComplex({ real_spec, ispec_spec }, hop_length, win_length, length);

        // STFT
        auto [stft_real, stft_imag] = SpecMath::stftComplex(y, n_fft, hop_length, win_length);

        // Update phase
        for (size_t k = 0; k < spec.size(); ++k) {
            for (size_t t = 0; t < n_frames; ++t) {
                angle[k][t] = std::atan2(stft_imag[k][t], stft_real[k][t]);
            }
        }
    }

    // Final complex spectrogram
    for (size_t k = 0; k < spec.size(); ++k) {
        for (size_t t = 0; t < n_frames; ++t) {
            real_spec[k][t] = spec[k][t] * std::cos(angle[k][t]);
            ispec_spec[k][t] = spec[k][t] * std::sin(angle[k][t]);
        }
    }

    return { real_spec, ispec_spec };
}