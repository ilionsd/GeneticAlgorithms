#include <array>
#include <charconv>
#include <iostream>
#include <random>
#include <optional>
#include <unordered_map>

#include <boost/program_options.hpp>

#include "genetic_algorithm/common/space.h"
#include "genetic_algorithm/cmn-ga.h"
#include "f1.h"
#include "f2.h"
#include "f3.h"
#include "f4.h"


template<typename CharT, typename T>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const std::valarray<T>& val);

template<typename T>
auto select_target(int target) {
    switch (target) {
        case 1: {
            auto [fitness, bounds] = f1<T>();
            return std::make_pair(std::function<T(const std::valarray<T>&)>(fitness), bounds);
        }
        case 2: {
            auto [fitness, bounds] = f2<T>();
            return std::make_pair(std::function<T(const std::valarray<T>&)>(fitness), bounds);
        }
        case 3: {
            auto [fitness, bounds] = f3<T>();
            return std::make_pair(std::function<T(const std::valarray<T>&)>(fitness), bounds);
        }
        case 4: {
            auto [fitness, bounds] = f4<T>();
            return std::make_pair(std::function<T(const std::valarray<T>&)>(fitness), bounds);
        }
    default:
        return std::make_pair(std::function<T(const std::valarray<T>&)>(), std::vector<std::pair<const T, const T>>());
    }
}

template<typename CharT>
auto parse_function(const std::basic_string<CharT> &str) {
    if (str.find("-f") != 0) {
        return make_pair(std::basic_string<CharT>(), std::basic_string<CharT>());
    }
    else {
        return make_pair(std::basic_string<CharT>("function"), std::basic_string<CharT>(str.substr(2)));
    }
}

template<typename CharT>
std::optional<int> maybe_int(const std::basic_string<CharT> &str) {
    int val;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), val);
    if (ec == std::errc{}) {
        return std::make_optional<int>();
    }
    return std::make_optional<int>(val);
}

auto main(int argc, char *argv[]) -> int {

    std::optional<int> target;

    try {
        namespace po = boost::program_options;
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("method", po::value<std::string>(), "optimization method")
        ;
        po::positional_options_description positional;
        positional.add("method", 1);

        auto parsed = po::command_line_parser(argc, argv)
            .options(desc)
            .positional(positional)
            .extra_parser(parse_function<char>)
            .run();
        po::variables_map vm;
        po::store(parsed, vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
        if (vm.count("function")) {
            auto str = vm["function"].as<std::string>();
            target = maybe_int(str);
            if (target) {
                std::cout << "Invalid target: "  << str << "\n";
                return 1;
            }
        }

    }
    catch (std::exception& e) {

    }


    auto [fitness, bounds] = select_target<double>(0);

    auto params = ::genetic_algorithm::cmn_parameters {};
    params.initialPopulation = 100;
    params.maxPopulation = 1000;
    params.crossoversPerGeneration = 11;
    params.mutationsPerGeneration = 7;
    params.nMin = 10;
    params.nCrowd = 20;
    params.nTry = 20;
    params.targetClosestAvg = 2;

    auto printEvery10 = [](std::uint64_t gen, std::size_t size, const std::vector<std::valarray<double>>& population, const std::vector<double>& fitness) {
        if (gen % 10 == 0) {
            std::cout << "\nGeneration " << gen << ", population size " << size;
            for (std::size_t k = 0; k < size; ++k) {
                std::cout << "\n" << population[k] << " -> " << fitness[k];
            }
            std::cout << std::endl;
        }
    };
    auto printPopulationSize = [](std::uint64_t gen, std::size_t size, const std::vector<std::valarray<double>>& population, const std::vector<double>& fitness) {
        std::cout << "\nGeneration " << gen << ", population size " << size;
    };
    auto cmnGA = params.make_optimizer(printPopulationSize);

    auto resolution = std::vector<std::uint32_t>{ 10000, 10000 };
    auto space = ::genetic_algorithm::common::make_space(bounds, resolution);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    auto [result, generation, real, vals] = cmnGA.optimize(gen, space, fitness);
    
    switch (result) {
        case ::genetic_algorithm::cmn_result::CONVERGED:
            std::cout << "\nOptimization completed on generation " << generation 
                      << " by converging";
            break;
        case ::genetic_algorithm::cmn_result::REACHED_GENERATION_LIMIT:
            std::cout << "\nOptimization completed on generation " << generation 
                      << " by reaching generation limit";
            break;
        case ::genetic_algorithm::cmn_result::REACHED_POPULATION_LIMIT:
            std::cout << "\nOptimization completed on generation " << generation 
                      << " by reaching population limit";
            break;
        default:
            std::cout << "Optimization finished on generation " << generation 
                      << " with unexpected result: " << static_cast<int>(result);
            break;
    }
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
