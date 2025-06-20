#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <omp.h>

struct generator {
private:
	int bins;

public:
    generator(const int& max) : bins(max) {	}
    int operator()(int n) {
		int num_factors = 0;

		// Start with 2 and increment by 1 to check each number
		for (int p = 2; p <= n; ++p) {
			// While p divides n, print p and divide n
			while (n % p == 0) {
				n = n / p;
				num_factors++;
			}
		}
		return (num_factors < bins - 1) ? num_factors : bins - 1;
	}
};

struct histogram {
	int bins, *data;

	histogram(int count) : bins(count) {
		// allocate memory for histogram
		data = (int*) malloc(sizeof(int) * count);

		// initialize histogram with 0's
		for (int b = 0; b < count; b++) {
			data[b] = 0;
		}
	}

	~histogram() {free(data); }

	void populate(int sample_size) {

		#pragma omp parallel
		{
			// initialize prime factors generator
			generator number_generator(bins);

			//creating an arrays of bins for each thread
			int* local_data = new int[bins];

			//initialiting local_data with 0
			for (int i = 0; i < bins; i++) local_data[i] = 0;
			
			#pragma omp for
			for (int i = 2; i < sample_size; i++) {
				// count number of prime factors for integer i
				int number_of_primes = number_generator(i);

				// each thread will update the corresponding bin of their array
				local_data[number_of_primes]++;
			}

			// loop for adding the data from each local thread to the final data array
			
			#pragma omp critical
			for (int i=0; i < bins; i++){
				// omp critical because some bins have the risk of being updated simultaneously 
				 
				data[i] += local_data[i];
			}
		}
	}

	void print() {
		int total = 0;
		for (int b = 0; b < bins; ++b) {
			total += data[b];
			std::cout << b << ":" << data[b] << std::endl;
		}
		std::cout << "total: " << total << std::endl;
	}
};

int main(int argc, char **argv)
{
	int num_bins = 10;
	int sample_ceiling = 50000;
	
	std::cout << "Bins: " << num_bins << ", sample ceiling: " << sample_ceiling  << std::endl;

	// initialize and empty histogram with 'num_bins' bins
	histogram h(num_bins);

	auto t1 = std::chrono::high_resolution_clock::now();

	// populate the histogram that was just created
	h.populate(sample_ceiling);

	auto t2 = std::chrono::high_resolution_clock::now();

	// print the contents of the histogram and the time it took populate it
	h.print();

	std::cout << "\ntime elapsed: " << std::chrono::duration<double>(t2 - t1).count() << " seconds." << std::endl;
}