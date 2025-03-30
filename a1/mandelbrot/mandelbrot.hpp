
#include <vector>
#include <complex>
#include <numeric>

#include "pixel.hpp"
#include "image.hpp"

class Mandelbrot {
    Image image;

    int max_iterations = 2048;

public:
    Mandelbrot(int rows, int cols, int max_iterations): image(rows, cols, {255, 255, 255}), max_iterations(max_iterations) { }

    void compute(int num_threads = 1) {
        worker(num_threads);
    }

    void worker(int num_threads, int thread_id=0)
    {
        
        unsigned char color = 0; // (254*(thread_id+1))/num_threads % 254; // use for your parallel code

        for (int y = 0; y < image.height; ++y) 
        {
            for (int x = 0; x < image.width; ++x) 
            {
                double dx = ((double)x / image.width - 0.75) * 2.0;
                double dy = ((double)y / image.height - 0.5) * 2.0;

                std::complex<double> c(dx, dy);

                if (check_pixel(c)) { 
                    image[y][x] = {color, color, color}; // {0, 0, 0} - black for a pixel inside of the mandelbrot set
                }
            }
        }
    }

    // Test if point c belongs to the Mandelbrot set
    bool check_pixel(std::complex<double> c)
    {
        std::complex<double> z(0, 0);
        for (int i=0; i<max_iterations; ++i) {
            z = z * z + c;
            if ( std::abs(z) > 4 ) return false;
        }

        return true;
    };

    void save_to_ppm(std::string filename){
        image.save_to_ppm(filename);
    }
};