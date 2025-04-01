#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <numeric>
#include <thread>

#include "helper.hpp"

using namespace std;

struct histogram {
	vector<int> data;
	
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

void append(const histogram& other) {
		if (data.size() != other.data.size()) {
				throw std::invalid_argument("Histograms must have the same size to append.");
		}
		for (size_t i = 0; i < data.size(); ++i) {
				data[i] += other.data[i];
		}
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
	vector<thread> threads;

	int num_threads = std::thread::hardware_concurrency();
	int print_level = 3; // 0: exec info + histogram total, 1: + histogram bins, 2: +exec time, 3: +bin info
	int sample_per_thread = sample_count/num_threads;
	int remaining_samples = sample_count - sample_per_thread*(num_threads-1);
	
	parse_args(argc, argv, num_threads, num_bins, sample_count, print_level);
	
	histogram h(num_bins);
  vector<histogram> hists(num_threads, histogram(num_bins));

	auto t1 = chrono::high_resolution_clock::now();

	for (int i=0; i<num_threads; i++){
		threads.push_back(thread(worker, sample_per_thread, std::ref(hists[i]), num_bins));
	}

	for(auto& thread : threads) {
		thread.join();
	}

	for (const auto& h_th : hists) {
			h.append(h_th);
	}
	// worker(sample_count, h, num_bins);

	auto t2 = chrono::high_resolution_clock::now();

	if ( print_level >= 0 ) cout << "Bins: " << num_bins << ", sample size: " << sample_count << ", threads: " << num_threads << endl;
	if ( print_level >= 3 ) h.print_bins(cout);
	if ( print_level >= 1 ) h.print_total(cout);
	if ( print_level >= 2 || print_level == -1 ) cout << chrono::duration<double>(t2 - t1).count() << endl;
}