#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <numeric>
#include <thread>
#include <mutex>

#include "helper.hpp"

using namespace std;
std::mutex m;

struct histogram {
	vector<int> data;
	
	histogram(int count) : data(count) { }

	void add(int i) {
		// locking and unlocking the resource that is common to all threads
		m.lock();
		++data[i];
		m.unlock();
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
    vector<thread> threads;

    int num_threads = std::thread::hardware_concurrency();
    int print_level = 3; // 0: exec info + histogram total, 1: + histogram bins, 2: +exec time, 3: +bin info

    parse_args(argc, argv, num_threads, num_bins, sample_count, print_level);

	// creating a variable to hold the number of samples per thread
    int sample_per_thread = sample_count / num_threads;
	// and the number of remaining samples for the last thread
    int remaining_samples = sample_count - sample_per_thread * (num_threads - 1);
	
	histogram h(num_bins);

	auto t1 = chrono::high_resolution_clock::now();

	// creating the first n-1 threads
	for (int i=0; i<num_threads-1; i++){
		threads.push_back(thread(worker, sample_per_thread, std::ref(h), num_bins));
	}

	// creating the last thread
	threads.push_back(thread(worker, remaining_samples, std::ref(h), num_bins));

	// joining the threads to ensure they finish before the main thread
	for(auto& thread : threads) {
		thread.join();
	}
	// worker(sample_count, h, num_bins);

	auto t2 = chrono::high_resolution_clock::now();

	if ( print_level >= 0 ) cout << "Bins: " << num_bins << ", sample size: " << sample_count << ", threads: " << num_threads << endl;
	if ( print_level >= 3 ) h.print_bins(cout);
	if ( print_level >= 1 ) h.print_total(cout);
	if ( print_level >= 2 || print_level == -1 ) cout << chrono::duration<double>(t2 - t1).count() << endl;
}