#pragma once
#include "common.hpp"
#include "spec_types.hpp"

class SpecMath {
public:
    static void fft(
        std::vector<double>& real,
        std::vector<double>& imag,
        bool inverse
    );

    static SpectrogramMatrix stft(
        const std::vector<double>& samples,
        size_t n_fft,
        size_t hop_length,
        size_t win_length
    );

    static ComplexSpectrogramMatrix stftComplex(
        const std::vector<double>& samples,
        size_t n_fft,
        size_t hop_length,
        size_t win_length
    );

    static std::vector<double> istft(
        const SpectrogramMatrix& spec,
        size_t hop_length,
        size_t win_length,
        size_t length
    );

    static std::vector<double> istftComplex(
        const ComplexSpectrogramMatrix& spec,
        size_t hop_length,
        size_t win_length,
        size_t length
    );

    static std::map<size_t, std::vector<double>> m_hannWindowCache;
    static std::vector<double> hannWindow(size_t win_length);

    static std::map<size_t, std::pair<std::vector<double>, std::vector<double>>> m_twiddlesCache;
	static std::pair<std::vector<double>, std::vector<double>> twiddles(size_t n_fft, bool inverse);

    static ComplexSpectrogramMatrix griffinLim(
        const SpectrogramMatrix& spec,
        size_t n_fft,
        size_t hop_length,
        size_t win_length,
        size_t iterations
    );
};
