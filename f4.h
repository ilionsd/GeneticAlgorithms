#ifndef F4_H_
#define F4_H_

#include <array>
#include <valarray>

template<typename T>
auto f4() {
    auto bounds = std::vector<std::pair<const T, const T>>{
        {-100, 100},
        {-100, 100}
    };

    std::array<std::valarray<T>, 5> B = {
        std::valarray<T>{-20, -20, 0.4, 0.02},
        std::valarray<T>{-5, -25, 0.2, 0.5},
        std::valarray<T>{0, 30, 0.7, 0.01},
        std::valarray<T>{30, 0, 1.0, 2.0},
        std::valarray<T>{30, -30, 0.05, 0.1}
    };
    auto f4 = [B](const std::valarray<T>& x) -> T {
        T val = 0.0;
        auto dims = x.size();
        for (std::size_t k = 0; k < B.size(); ++k) {
            T denominator = 1.0;
            for (std::size_t dim = 0; dim < x.size(); ++dim) {
                denominator += std::pow(x[dim] - B[k][dim], 2);
            }
            val += B[k][dims] / (denominator * B[k][dims + 1]);
        }
        return val;
    };

    return std::make_pair(f4, bounds);
}

#endif //F4_H_
