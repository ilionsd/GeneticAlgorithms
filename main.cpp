#include <iostream>
#include <random>

#include "genetic_algorithm/space.h"
#include "genetic_algorithm/cmn-ga.h"

auto main() -> int {
    std::random_device rd;
    auto params = ::genetic_algorithm::cmn_parameters {};
    auto cmnGA = params.make_optimizer<std::mt19937>(rd, std::cout);

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

    auto [real, vals] = cmnGA.optimize(space, fitness);
    
    return 0;
}
