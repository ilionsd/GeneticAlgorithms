#include <array>
#include <future>
#include <iostream>
#include <numbers>
#include <random>
#include <optional>
#include <thread>
#include <variant>
#include <unordered_set>
#include <unordered_map>

#include <boost/program_options.hpp>

#include <json/json.h>

#include "genetic_algorithm/space.h"
#include "genetic_algorithm/algorithms.h"
#include "f1.h"
#include "f2.h"
#include "f3.h"
#include "f4.h"
#include "utility.h"


using target_type = double;
using fitness_function = std::function<target_type(const std::valarray<target_type> &)>;
using bounds_type = std::vector<std::pair<const target_type, const target_type>>;

using namespace ::genetic_algorithm;
using variants_of_parameters = std::variant<
    dummy::parameters, 
    cmn::parameters
>;


template<typename CharT, typename T>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const std::valarray<T>& val);


template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;


auto parse_function(const std::string &str) -> std::pair<std::string, std::string>;


auto print_population_size(
    std::uint64_t gen, 
    std::size_t size, 
    const std::vector<std::valarray<target_type>> &population, 
    const std::vector<target_type> &fitness
) -> void;


auto print_every10(
    std::uint64_t gen, 
    std::size_t size, 
    const std::vector<std::valarray<target_type>> &population, 
    const std::vector<target_type> &fitness
) -> void;


template<typename T1, typename T2, class Params, class... Ts>
auto optimize(
    std::function<T1(const std::valarray<T1> &)> &&fitness, 
    const space<T1, T2> &space, 
    const Params &params, 
    const Ts... args
) -> std::promise<result<T1, Params>>;


template<typename T, class Params>
auto handle_result(result<T, Params> &&result);


auto main(int argc, char *argv[]) -> int {

    std::unordered_map<int, std::pair<fitness_function, bounds_type>> availableTargets {
        { 1, f1<target_type>() },
        { 2, f2<target_type>() },
        { 3, f3<target_type>() },
        { 4, f4<target_type>() }
    };
    
    std::unordered_map<std::string, variants_of_parameters> defaultParams {
        { "dummy", dummy::parameters {
            .maxPopulation = 1000000
        } },
        { "cmn", cmn::parameters {
            .initialPopulation = 100,
            .maxPopulation = 1000,
            .crossoversPerGeneration = 11,
            .mutationsPerGeneration = 7,
            .nMin = 10,
            .nCrowd = 20,
            .nTry = 20,
            .targetClosestAvg = 2
        } }
    };
    
    std::optional<std::string> oMethod;
    std::optional<int> oTarget;
    std::optional<double> oPrecision;
    std::optional<std::uint64_t> oGenerationLimit;
    try {
        namespace po = boost::program_options;
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("method", po::value<std::string>(), "optimization method")
            ("precision,p", po::value<double>(), "precision, mesh size")
            ("limit,l", po::value<std::uint64_t>(), "generation limit, max generations count")
        ;
        po::positional_options_description positional;
        positional.add("method", 1);

        auto parsed = po::command_line_parser(argc, argv)
            .options(desc)
            .positional(positional)
            .extra_parser(parse_function)
            .run();
        po::variables_map vm;
        po::store(parsed, vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        if (vm.count("method")) {
            auto str = vm["method"].as<std::string>();
            oMethod = str;
        }
        if (vm.count("function")) {
            auto str = vm["function"].as<std::string>();
            oTarget = optional_parse<int>(str);
            if (!oTarget) {
                throw std::logic_error("\nInvalid target: " + str);
            }
        }
        if (vm.count("precision")) {
            auto str = vm["precision"].as<std::string>();
            oPrecision = optional_parse<double>(str);
        }
        if (vm.count("limit")) {
            auto str = vm["limit"].as<std::string>();
            oGenerationLimit = optional_parse<std::uint64_t>(str);
        }
    }
    catch (std::exception& e) {
        std::cerr << "\nAn error occured: " << e.what() << std::endl;
        return 1;
    }

    auto target = oTarget.value();
    auto selectedTargetEntry = availableTargets.find(target);
    if (availableTargets.cend() == selectedTargetEntry) {
        std::cerr << "\nInvalid target: " << target << "\nAvailable targets: ";
        for (const auto &[key, value]: availableTargets) {
            std::cerr << key << ' ';
        }
        std::cerr << std::endl;
        return 1;
    }
    auto [fitness, bounds] = selectedTargetEntry->second;

    auto method = oMethod.value();
    auto defaultParamsEntry = defaultParams.find(method);
    if (defaultParams.cend() == defaultParamsEntry) {
        std::cerr << "\nUnsupported method: " << method << "\nSupported methods: ";
        for (const auto &[key, value]: defaultParams) {
            std::cerr << key << ' ';
        }
        std::cerr << std::endl;
        return 1;
    }
    auto vParams = defaultParamsEntry->second;

    Json::Value configuration;
    Json::CharReaderBuilder builder;
    //builder["collectComments"] = true;
    JSONCPP_STRING errs;
    bool success = Json::parseFromStream(builder, std::cin, &configuration, &errs);
    if (success) {
        Json::Value methodConfig = configuration[method];
        std::visit([&methodConfig](auto &params) {
            common::maybe_set<std::decay_t<decltype(params)>> {
                params
            }.from_json(methodConfig);
        }, vParams);
    }
    else {
        std::visit([](auto &params) {
            std::cin >> params;
        }, vParams);
    }

    auto generationLimit = oGenerationLimit.value_or(1000);
    auto vOptimizer = std::visit(overloaded {
        [](const dummy::parameters &params) {
            return params.make_optimizer(print_population_size);
        },
        [generationLimit](const cmn::parameters &params) {
            return params.make_optimizer(generationLimit, print_population_size);
        }
    }, vParams);

    auto precision = oPrecision.value_or(1.0) / std::numbers::sqrt2 * 2.0;
    auto resolution = std::vector<std::uint32_t>(bounds.size());
    for (std::size_t k = 0; k < bounds.size(); ++k) {
        auto count = std::abs(bounds[k].first - bounds[k].second) / precision;
        resolution[k] = static_cast<std::uint32_t>(std::ceil(count)) + 1;
    }

    auto space = make_space(bounds, resolution);

    std::random_device rd;
    std::mt19937_64 gen { rd() };
    auto vResult = std::visit([&gen, &space, fitness](const auto &optimizer) { 
        return optimizer.optimize(gen, space, fitness);
    }, vOptimizer);

    std::visit([](const auto &result) {
        handle_result(result);
    }, vResult);
    
    return 0;
}


template<typename CharT, typename T>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const std::valarray<T>& val) {
    os << '[' << ' ';
    for (std::size_t k = 0; k < val.size(); ++k) {
        os << val[k] << ' ';
    }
    return os << ']';
}


auto parse_function(const std::string &str) -> std::pair<std::string, std::string> {
    if (str.find("-f") != 0) {
        return make_pair(std::string(), std::string());
    }
    else {
        return make_pair(std::string("function"), std::string(str.substr(2)));
    }
}


auto print_population_size(
    std::uint64_t gen, 
    std::size_t size, 
    const std::vector<std::valarray<target_type>> &population, 
    const std::vector<target_type> &fitness
) -> void {
    std::cout << "\nGeneration " << gen << ", population size " << size;
}


auto print_every10(
    std::uint64_t gen, 
    std::size_t size, 
    const std::vector<std::valarray<target_type>> &population, 
    const std::vector<target_type> &fitness
) -> void {
    if (gen % 10 == 0) {
        std::cout << "\nGeneration " << gen << ", population size " << size;
        for (std::size_t k = 0; k < size; ++k) {
            std::cout << "\n" << population[k] << " -> " << fitness[k];
        }
        std::cout << std::endl;
    }
}


template<typename T1, typename T2, class Params, class... Ts>
auto optimize(
    std::function<T1(const std::valarray<T1> &)> &&fitness,
    const space<T1, T2> &space, 
    const Params &params, 
    const Ts... args
) -> std::thread {
    auto optimizer = params.make_optimizer(args...);
    auto thread = std::thread { 
        [&optimizer](const auto &space, auto &&fitness) {
            std::random_device rd;
            std::mt19937_64 gen { rd() };
            auto result = optimizer.optimize(gen, space, fitness);
            futureResult.set_value_at_thread_exit(result);
        }, space, fitness
    };

    return thread;
}


template<typename T>
auto handle_result(result<T, dummy::parameters> &&result) {
    std::cout << "\nDummy" << std::endl;
}

template<typename T>
auto handle_result(result<T, cmn::parameters> &&result) {
    switch (result.status) {
    case status::IN_PROGRESS:
        std::cout << "\nOptimization is ongoing: generation " << result.generation;
        break;
    case status::CONVERGED:
        std::cout << "\nOptimization completed on generation " << result.generation
                << " by converging";
        break;
    case status::REACHED_GENERATION_LIMIT:
        std::cout << "\nOptimization completed on generation " << result.generation
                << " by reaching generation limit";
        break;
    case status::REACHED_POPULATION_LIMIT:
        std::cout << "\nOptimization completed on generation " << result.generation
                << " by reaching population limit";
        break;
    default:
        std::cout << "Optimization finished on generation " << result.generation 
                << " with unexpected result: " << static_cast<int>(result.status);
        break;
    }
    std::cout << "\nFinal population size: " << result.population.size();
    for (std::size_t k = 0; k < result.population.size(); ++k) {
        std::cout << "\n" << result.population[k] << " -> " << result.values[k];
    }
    std::cout << std::endl;
}
