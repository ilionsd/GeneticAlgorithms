#include <array>
#include <iostream>
#include <random>

#include "genetic_algorithm/space.h"
#include "genetic_algorithm/cmn-ga.h"

template<typename CharT, typename T>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const std::valarray<T>& val);

auto main() -> int {

    auto f1 = [](const std::valarray<double>& x) -> double {
        auto val = std::sin(5.1 * M_PI * x[0] + 0.5);
        return std::pow(val, 6);
    };
    auto f2 = [f1](const std::valarray<double>& x) -> double {
        auto val = 4.0 * std::log(2.0) * std::pow(x[0] - 0.0667, 2) / 0.64;
        return std::exp(-val) * f1(x);
    };
    std::array<std::valarray<int>, 25> A;
    for (std::size_t i = 0; i < 5; ++i) {
        for (std::size_t j = 0; j < 5; ++j) {
            auto x = static_cast<int>(i) - 2;
            auto y = static_cast<int>(j) - 2;
            A[i * 5 + j] = std::move(std::valarray<int>({x, y}) * 16);
        }
    }
    auto f3 = [A](const std::valarray<double>& x) -> double {
        auto val = 0.0;
        for (std::size_t k = 0; k < A.size(); ++k) {
            auto denominator = static_cast<double>(k);
            for (std::size_t dim = 0; dim < x.size(); ++dim) {
                denominator += std::pow(x[dim] - A[k][dim], 6);
            }
            val += 1 / denominator;
        }
        return val;
    };
    std::array<std::valarray<double>, 5> B = {
        std::valarray<double>{-20, -20, 0.4, 0.02},
        std::valarray<double>{-5, -25, 0.2, 0.5},
        std::valarray<double>{0, 30, 0.7, 0.01},
        std::valarray<double>{30, 0, 1.0, 2.0},
        std::valarray<double>{30, -30, 0.05, 0.1}
    };
    auto f4 = [B](const std::valarray<double>& x) -> double {
        auto val = 0.0;
        auto dims = x.size();
        for (std::size_t k = 0; k < B.size(); ++k) {
            auto denominator = 1.0;
            for (std::size_t dim = 0; dim < x.size(); ++dim) {
                denominator += std::pow(x[dim] - B[k][dim], 2);
            }
            val += B[k][dims] / (denominator * B[k][dims + 1]);
        }
        return val;
    };

    auto params = ::genetic_algorithm::cmn_parameters {};
    auto cmnGA = params.make_optimizer();

    auto resolution = std::vector<std::uint32_t>{ 1000, 1000 };
    auto bounds = std::vector<std::pair<const double, const double>>{
        {-1, 1},
        {-1, 1}
    };
    auto space = ::genetic_algorithm::make_space(bounds, resolution);
    auto fitness = std::function<double(const std::valarray<double>&)> (
        [](const std::valarray<double>& val) {
            return 0.0;
        }
    );

    auto onNextGeneration = [](std::uint64_t gen, std::size_t size, const std::vector<std::valarray<double>>& population, const std::vector<double>& fitness) {
        if (gen % 10 == 0) {
            std::cout << "\nGeneration " << gen << ", population size " << size;
            for (std::size_t k = 0; k < size; ++k) {
                std::cout << "\n" << population[k] << " -> " << fitness[k];
            }
            std::cout << std::endl;
        }
    };

    std::random_device rd;
    std::mt19937_64 gen(rd());
    auto [real, vals] = cmnGA.optimize(gen, space, fitness, onNextGeneration);
    std::cout << "\nFinal population size: " << real.size();
    for (std::size_t k = 0; k < real.size(); ++k) {
        std::cout << "\n" << real[k] << " -> " << vals[k];
    }
    std::cout << std::endl;
    return 0;
}

template<typename CharT, typename T>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const std::valarray<T>& val) {
    os << "[ ";
    for (std::size_t k = 0; k < val.size(); ++k) {
        os << val[k] << " ";
    }
    return os << "]";
}
