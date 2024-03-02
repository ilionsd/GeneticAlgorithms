#include <array>
#include <iostream>
#include <random>
#include <optional>
#include <variant>
#include <unordered_set>
#include <unordered_map>

#include <boost/program_options.hpp>

#include <json/json.h>

#include "genetic_algorithm/common/space.h"
#include "genetic_algorithm/dummy.h"
#include "genetic_algorithm/cmn-ga.h"
#include "f1.h"
#include "f2.h"
#include "f3.h"
#include "f4.h"
#include "utility.h"


template<typename CharT, typename T>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const std::valarray<T>& val);

template<typename CharT>
auto parse_function(const std::basic_string<CharT> &str) {
    if (str.find("-f") != 0) {
        return make_pair(std::basic_string<CharT>(), std::basic_string<CharT>());
    }
    else {
        return make_pair(std::basic_string<CharT>("function"), std::basic_string<CharT>(str.substr(2)));
    }
}

template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

using target_type = double;
using fitness_function = std::function<target_type(const std::valarray<target_type> &)>;
using bounds_type = std::vector<std::pair<const target_type, const target_type>>;

using variants_of_parameters = std::variant<::genetic_algorithm::dummy::parameters, ::genetic_algorithm::cmn::parameters>;


auto main(int argc, char *argv[]) -> int {

    std::unordered_map<int, std::pair<fitness_function, bounds_type>> availableTargets {
        { 1, f1<target_type>() },
        { 2, f2<target_type>() },
        { 3, f3<target_type>() },
        { 4, f4<target_type>() }
    };
    
    std::unordered_map<std::string, variants_of_parameters> defaultParams {
        { "dummy", ::genetic_algorithm::dummy::parameters {
            .maxPopulation = 1000000
        } },
        { "cmn", ::genetic_algorithm::cmn::parameters {
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
    std::optional<std::string> method;
    std::optional<int> target;
    std::optional<double> precision;
    std::optional<std::uint64_t> generationLimit;

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
            .extra_parser(parse_function<char>)
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
            method = str;
        }
        if (vm.count("function")) {
            auto str = vm["function"].as<std::string>();
            target = optional_parse<int>(str);
            if (target) {
                std::cerr << "\nInvalid target: "  << str << std::endl;
                return 1;
            }
        }
        if (vm.count("precision")) {
            auto str = vm["precision"].as<std::string>();
            precision = optional_parse<double>(str);
        }
        if (vm.count("limit")) {
            auto str = vm["limit"].as<std::string>();
            generationLimit = optional_parse<std::uint64_t>(str);
        }
    }
    catch (std::exception& e) {
        std::cerr << "\nAn error occured: " << e.what() << std::endl;
        //std::cout << desc << std::endl;
        return 1;
    }

    auto selectedTargetEntry = availableTargets.find(target.value());
    if (availableTargets.cend() == selectedTargetEntry) {
        std::cerr << "\nInvalid target: " << method.value() << "\nAvailable targets: ";
        for (const auto &[key, value]: availableTargets) {
            std::cerr << key << ' ';
        }
        std::cerr << std::endl;
        return 1;
    }
    auto [fitness, bounds] = selectedTargetEntry->second;

    auto defaultParamsEntry = defaultParams.find(method.value());
    if (defaultParams.cend() == defaultParamsEntry) {
        std::cerr << "\nUnsupported method: " << method.value() << "\nSupported methods: ";
        for (const auto &[key, value]: defaultParams) {
            std::cerr << key << ' ';
        }
        std::cerr << std::endl;
        return 1;
    }
    auto paramsVariant = defaultParamsEntry->second;

    Json::Value configuration;
    Json::CharReaderBuilder builder;
    //builder["collectComments"] = true;
    JSONCPP_STRING errs;
    bool success = Json::parseFromStream(builder, std::cin, &configuration, &errs);
    if (success) {
        Json::Value methodConfig = configuration[method.value()];
        std::visit([&methodConfig](auto &params) {
            ::genetic_algorithm::common::maybe_set<std::decay_t<decltype(params)>> {
                params
            }.from_json(methodConfig);
        }, paramsVariant);
    }
    else {
        std::visit([](auto &params) {
            std::cin >> params;
        }, paramsVariant);
    }

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

    auto optimizerVariant = std::visit(overloaded {
        [](const ::genetic_algorithm::dummy::parameters &params) {
            return params.make_optimizer();
        },
        [generationLimit, printPopulationSize](const ::genetic_algorithm::cmn::parameters &params) {
            return params.make_optimizer(generationLimit.value(), printPopulationSize);
        }
    }, paramsVariant);

    auto resolution = std::vector<std::uint32_t>(bounds.size());
    for (std::size_t k = 0; k < bounds.size(); ++k) {
        auto count = std::abs(bounds[k].first - bounds[k].second) / precision.value();
        resolution[k] = static_cast<std::uint32_t>(std::ceil(count)) + 1;
    }

    auto space = ::genetic_algorithm::common::make_space(bounds, resolution);

    std::random_device rd;
    std::mt19937_64 gen { rd() };
    auto resultVariant = std::visit([&gen, &space, fitness](const auto &optimizer) { 
        return optimizer.optimize(gen, space, fitness);
    }, optimizerVariant);

    std::visit(overloaded {
        [](const ::genetic_algorithm::dummy::result<target_type> &result) {
            std::cout << "\nDummy" << std::endl;
        },
        [](const ::genetic_algorithm::cmn::result<target_type> &result) {
            switch (result.status) {
            case ::genetic_algorithm::cmn::status::CONVERGED:
                std::cout << "\nOptimization completed on generation " << result.generation
                        << " by converging";
                break;
            case ::genetic_algorithm::cmn::status::REACHED_GENERATION_LIMIT:
                std::cout << "\nOptimization completed on generation " << result.generation
                        << " by reaching generation limit";
                break;
            case ::genetic_algorithm::cmn::status::REACHED_POPULATION_LIMIT:
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
    }, resultVariant);
    
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
