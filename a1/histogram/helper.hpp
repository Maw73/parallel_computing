#include <iostream>
#include <string>
#include <random>


struct generator {
    // let's generate random number like this for consistent tests
    generator(const int& max) : dist(0, max-1) {	}
    int operator()() {
        return dist(engine);
    }
private:
    int max;
    std::minstd_rand engine; 
    std::uniform_int_distribution<int> dist;
};

// Parsing arguments
void parse_args(int argc, char **argv, int& num_threads, int &num_bins, int &sample_size, int& print_level) {
    std::string usage("Usage: --num-threads <integer> --num-bins <integer> --samlpe-size <integer> --print-level <integer>");

    for (int i = 1; i < argc; ++i) {
        if ( std::string(argv[i]).compare("--num-threads") == 0 ) {
            num_threads=std::stoi(argv[++i]);
        } else if ( std::string(argv[i]).compare("--num-bins") == 0 ) {
            num_bins=std::stoi(argv[++i]);
        } else if ( std::string(argv[i]).compare("--sample-size") == 0 ) {
            sample_size=std::stoi(argv[++i]);
        } else if ( std::string(argv[i]).compare("--print-level") == 0 ) {
            print_level=std::stoi(argv[++i]);
        } else if (  std::string(argv[i]).compare("--help") == 0 ) {
            std::cout << usage << std::endl;
            exit(-1);
        }
    }
};