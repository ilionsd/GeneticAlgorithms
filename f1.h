#ifndef F1_H_
#define F1_H_

#include <numbers>
#include <valarray>

template<class T>
auto f1() {
    auto bounds = std::vector<std::pair<const T, const T>>{
        {0, 1}
    };
    auto fitness = [](const std::valarray<T>& x) -> T {
        auto phi = static_cast<T>(5.1 * std::numbers::pi_v<T>) * x[0] + static_cast<T>(0.5);
        return std::pow(std::sin(phi), static_cast<T>(6));
    };
    return std::make_pair(fitness, bounds);
}

#endif //F1_H_
