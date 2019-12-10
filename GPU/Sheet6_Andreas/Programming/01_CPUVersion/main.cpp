#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <vector>
#include <numeric>

constexpr size_t ITERATIONS = 100;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Pass DIM as param\n";
        return 1;
    }
    int vectorSize = atoi(argv[1]);
    std::vector<float> vec(vectorSize);
    std::iota(vec.begin(), vec.end(), 1);

    volatile float result = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t count = 0; count < ITERATIONS; count++)
    {
        result = 0.0;
        for (const auto &it : vec)
            result += it;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> time = end - start;

    // Calculate the bandwith = Bytes / elapsed time
    double bytes = sizeof(float) * vec.size();
    double timeInSeconds = (time.count() / ITERATIONS) / 1e6;
    double bandwith = bytes / timeInSeconds;
    bandwith /= 1e9; // Giga

    std::cout << "CPU Version: Iterations, Problem Size, Runtime[microseconds], Bandwith[GB/s], Result\n";
    std::cout << ITERATIONS << ", " << vec.size() << ", " << time.count() / ITERATIONS << ", " << bandwith << ", " << result << std::endl;

    return 0;
}