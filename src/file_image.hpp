#pragma once
#include "common.hpp"

class FileImage {
public:
    unsigned int width = 0, height = 0;
    std::vector<uint8_t> data; // RGBA format (8 bits per channel)

    void save(const std::string& path) const;
    static FileImage load(const std::string& path);
};
