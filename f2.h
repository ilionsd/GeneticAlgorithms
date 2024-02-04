#ifndef F2_H_
#define F2_H_

#include <valarray>

#include "f1.h"

template<class T>
auto f2() {
    auto [f, bounds] = f1<T>();
    auto fitness = [f](const std::valarray<T>& x) -> T {
        auto val = static_cast<T>(4.0 * std::log(2.0) / 0.64) * std::pow(x[0] - static_cast<T>(0.0667), 2) ;
        return std::exp(-val) * f(x);
    };
    return std::make_pair(fitness, bounds);
}

#endif //F2_H_
