
#include <vector>
#include <complex>
#include <numeric>
#include <thread>

#include "pixel.hpp"
#include "image.hpp"

using namespace std;

class Mandelbrot {
    Image image;

    int max_iterations = 2048;

//1st thread 0-127 x 0:95
//2nd thread 128-255 x 96:191
//3rd thread 256-383 x 192:287
//4th thread 384-511 x 288:383

public:
    Mandelbrot(int rows, int cols, int max_iterations): image(rows, cols, {255, 255, 255}), max_iterations(max_iterations) { }

    void compute(int num_threads) {
        vector<thread> threads;

        for (int i=0; i<num_threads; i++)
        {
            threads.push_back(thread([this, num_threads, i]() { worker(num_threads, i); }));
        }

        for (int i=0; i<num_threads; i++)
        {
            threads[i].join();
        }
        
    }

    void worker(int num_threads, int thread_id)
    {
        double h = image.height/num_threads;
        double from = thread_id*h;
        double to = (thread_id+1)*h;

        unsigned char color = (254*(thread_id+1))/num_threads % 254; // use for your parallel code

        //for (int y = 0; y < image.height; ++y)
        for (int y = from; y < to; ++y) 
        {
            for (int x = 0; x < image.width; ++x) 
            {
                double dx = ((double)x / image.width - 0.75) * 2.0;
                double dy = ((double)y / image.height - 0.5) * 2.0;

                complex<double> c(dx, dy);

                if (check_pixel(c)) { 
                    image[y][x] = {color, color, color}; // {0, 0, 0} - black for a pixel inside of the mandelbrot set
                }
            }
        }
    }

    // Test if point c belongs to the Mandelbrot set
    bool check_pixel(complex<double> c)
    {
        std::complex<double> z(0, 0);
        for (int i=0; i<max_iterations; ++i) {
            z = z * z + c;
            if ( abs(z) > 4 ) return false;
        }

        return true;
    };

    void save_to_ppm(string filename){
        image.save_to_ppm(filename);
    }
};