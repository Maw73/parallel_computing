#include <iostream>

namespace helper {
    /**
     * Parsing arguments
     * 
    */
    void parse_args(int argc, char **argv, int& num_threads, int& height, int& width, int& max_iterations) {
        std::string usage("Usage: --num-threads <integer> --height <integer> --width <integer>  --max-iterations <integer>");

        for (int i = 1; i < argc; ++i) {
            if ( std::string(argv[i]).compare("--num-threads") == 0 ) {
                num_threads=std::stoi(argv[++i]);
            } else if ( std::string(argv[i]).compare("--height") == 0 ) {
                height=std::stoi(argv[++i]);
            } else if ( std::string(argv[i]).compare("--width") == 0 ) {
                width=std::stoi(argv[++i]);
            } else if ( std::string(argv[i]).compare("--max-iterations") == 0 ) {
                max_iterations=std::stoi(argv[++i]);            
            } else if (  std::string(argv[i]).compare("--help") == 0 ) {
                std::cout << usage << std::endl;
                exit(-1);
            }
        }
    };
}