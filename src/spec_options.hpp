#pragma once
#include "common.hpp"

// Parameters for spectrogram generation and inversion
struct SpectrogramOptions {
    size_t n_fft = 2048;
    size_t hop_length = 512;
    size_t win_length = 0;       // 0 means same as n_fft
    float min_db = -80.0f;
    float max_db = 0.0f;
    size_t griffin_iter = 60;
    size_t sample_rate = 44100;
};
