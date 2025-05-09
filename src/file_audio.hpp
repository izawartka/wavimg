#pragma once
#include "common.hpp"

class FileAudio {
public:
    uint32_t sampleRate = 0;
    std::vector<double> samples;

    static FileAudio load(const std::string& path);
    void save(const std::string& path) const;
};
