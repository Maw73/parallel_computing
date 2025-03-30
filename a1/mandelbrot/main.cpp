
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <tuple>
#include <time.h>
#include <cmath>
#include <complex>
#include <chrono>

#include "helper.hpp"
#include "mandelbrot.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // arguments 
    int num_threads = 1;
    int max_iterations = 2048;
    
    // height and width of the output image
    int width = 512, height = 384;

    helper::parse_args(argc, argv, num_threads, height, width, max_iterations);

    std::cout << "Executing with " << num_threads << " threads\n";
    std::cout << "Mandelbrot image [" << width << "x" << height << "], max iterations " << max_iterations << "\n";

    // Generate Mandelbrot set in this image
    Mandelbrot mb(height, width, max_iterations); 

    auto t1 = chrono::high_resolution_clock::now();

    mb.compute(num_threads);

    auto t2 = chrono::high_resolution_clock::now();

    // save image
    mb.save_to_ppm("mandelbrot.ppm");

    cout << chrono::duration<double>(t2 - t1).count() << endl;

    return 0;
}