#include "file_audio.hpp"
#include <AudioFile.h>

FileAudio FileAudio::load(const std::string& path) {
    FileAudio ah;
    AudioFile<double> af;

    if (!af.load(path)) throw std::runtime_error("Failed to load audio: " + path);

    ah.sampleRate = af.getSampleRate();
    size_t sampleCount = af.getNumSamplesPerChannel();

    ah.samples.resize(sampleCount);
    // only mono (first channel)
    for (size_t i = 0; i < sampleCount; ++i) {
        ah.samples[i] = af.samples[0][i];
    }

    return ah;
}

void FileAudio::save(const std::string& path) const {
    AudioFile<double> af;
    af.setSampleRate(sampleRate);
    af.setNumChannels(1); // mono
    af.setAudioBufferSize(1, samples.size());
    af.samples[0] = samples;

    if (!af.save(path)) throw std::runtime_error("Failed to save audio: " + path);
}
