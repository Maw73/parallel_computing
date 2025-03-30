#include <iostream>
#include <vector>
#include <fstream>

#include "pixel.hpp"

struct Image {
    // Pixel** data = nullptr;
    std::vector<Pixel> data;
    
    int height;
    int width;
    
    Image (int height, int width, const Pixel& val = {0, 0, 0}) 
        : height(height), width(width), data(height*width, val )
    { }

    ~Image() {
        clear();
        height = 0;
        width = 0;
    }

    void clear() {
        data.clear();
    }

    Pixel* operator[](unsigned int row) {
        return &data[row*width];
    }

    void save_to_ppm(std::string filename) {
        std::ofstream ofs(filename, std::ofstream::out);
        ofs << "P3" << std::endl;
        ofs << width << " " << height << std::endl;
        ofs << 255 << std::endl;

        for (int row = 0; row < height; row++)
        {
            for (int col = 0; col < width; col++)
            {
                ofs << " " << (unsigned int)this->operator[](row)[col][0] << " " << (unsigned int)this->operator[](row)[col][1] << " " << (unsigned int)this->operator[](row)[col][2] << "\n";
            }
        }
        ofs.close();
    }

    std::string to_string() {
        std::stringstream ss;
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                for (int channel = 0; channel < 3; channel++)
                    ss << this->operator[](row)[col][channel] << " ";
            }
            ss << std::endl;
        }
        return ss.str();
    }
};