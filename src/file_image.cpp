#include "file_image.hpp"
#include "lodepng.h"

void FileImage::save(const std::string& path) const {
    if (width == 0 || height == 0 || data.size() != width * height * 4) {
        throw std::runtime_error("FileImage is empty or has invalid dimensions");
    }

    unsigned int err = lodepng::encode(
        path.c_str(),
        data,
        width,
        height
    );

    if (err) {
        throw std::runtime_error("PNG encode error: " + std::to_string(err));
    }
}

FileImage FileImage::load(const std::string& path) {
    FileImage img;
    unsigned int err = lodepng::decode(
        img.data,
        img.width,
        img.height,
        path.c_str()
    );

    if (err) {
        throw std::runtime_error("PNG decode error: " + std::to_string(err));
    }

    if (img.width == 0 || img.height == 0) {
        throw std::runtime_error("FileImage is empty or has invalid dimensions");
    }

    return img;
}
