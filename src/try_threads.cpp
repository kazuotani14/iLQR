#include <iostream>
#include <thread>
#include <chrono>
#include <functional>
#include <memory>
#include <vector>
#include <algorithm>
#include <numeric>
#include "eigen/Eigen/Core"

Eigen::MatrixXd mat(1000, 1000);

void fill_range(int start, int end)
{
    for(int i=start; i<end; i++)
    {
        for(int j=0; j<mat.cols(); j++)
        {
            if((i+j)%2 != 0)
                mat(i,j) = 1.0;
        }
    }
}

void fill_all()
{
    for(int i=0; i<mat.rows(); i++)
    {
        for(int j=0; j<mat.cols(); j++)
        {
            if((i+j)%2 == 0)
                mat(i,j) = 0.0;
            if((i+j)%2 != 0)
                mat(i,j) = 1.0;
        }
    }
}

long int fill_with_threads(int n_threads)
{
    mat.setZero();
    auto start = std::chrono::system_clock::now();

    std::vector<std::thread> threads(n_threads);

    for(int i=0; i<n_threads; i++)
    {
        int start = (mat.rows()/n_threads)*i;
        int end = (mat.rows()/n_threads)*(i+1);
        // std::cout << start << " " << end << std::endl;
        threads[i] = std::thread(std::bind(fill_range, start, end));
    }

    for(std::thread& t : threads) t.join();

    auto now = std::chrono::system_clock::now();
    long int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    // std::cout << n_threads << " threads took: " << elapsed/1000. << " seconds." << std::endl;
    return elapsed;
}

void run_and_average(int n_threads)
{
    int n_trials = 100;
    std::vector<long int> times(n_trials);
    for(int i=0; i<n_trials; i++)
    {
        times[i] = fill_with_threads(n_threads);
    }
    double average = std::accumulate(times.begin(), times.end(), 0.0)/times.size();
    std::cout << n_threads << " threads took: " << average/1000. << " seconds on average" << std::endl;
}


//g++ -o try_threads try_threads.cpp -std=c++11
int main()
{
    // 25 threads seems to do the best - why?
    std::vector<double> n_threads = {1, 2, 5, 10, 25, 50, 100, 500, 1000};
    for(const auto& num : n_threads)
    {
        // fill_with_threads(num);
        run_and_average(num);
    }

    return 0;
}
