#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <numeric>
#include <thread>

#include "helper.hpp"

using namespace std;

struct histogram {
	std::vector<int> data;
	
	histogram(int count) : data(count) { }

	void add(int i) {
		++data[i];
	}

	int& get(int i)	{
		return data[i];
	}

	void print_total(std::ostream& str) {
		str << "total:" << accumulate(data.begin(), data.end(), 0) << "\n";
	}

	void print_bins(std::ostream& str) {
		for (size_t i = 0; i < data.size(); ++i) str << i << ":" << data[i] << "\n";
	}
};

void worker(int sample_count, histogram& h, int num_bins) 
{
	long count = 0.0;

	generator gen(num_bins);

	while (sample_count--) {
		int next = gen();
		h.add(next);
		count++;
	}
}

int main(int argc, char **argv)
{
	int num_bins = 10;
	int sample_count = 30000000;

	int num_threads = std::thread::hardware_concurrency();
	int print_level = 3; // 0: exec info + histogram total, 1: + histogram bins, 2: +exec time, 3: +bin info
	parse_args(argc, argv, num_threads, num_bins, sample_count, print_level);
	
	histogram h(num_bins);

	auto t1 = chrono::high_resolution_clock::now();

	worker(sample_count, h, num_bins);

	auto t2 = chrono::high_resolution_clock::now();

	if ( print_level >= 0 ) cout << "Bins: " << num_bins << ", sample size: " << sample_count << ", threads: " << num_threads << endl;
	if ( print_level >= 3 ) h.print_bins(cout);
	if ( print_level >= 1 ) h.print_total(cout);
	if ( print_level >= 2 || print_level == -1 ) cout << chrono::duration<double>(t2 - t1).count() << endl;
}