#ifndef F3_H_
#define F3_H_

#include <array>
#include <valarray>

template<typename T>
auto f3() {
    auto bounds = std::vector<std::pair<const T, const T>>{
        {-100, 100},
        {-100, 100}
    };

    std::array<std::valarray<int>, 25> A;
    for (std::size_t i = 0; i < 5; ++i) {
        for (std::size_t j = 0; j < 5; ++j) {
            auto x = static_cast<int>(i) - 2;
            auto y = static_cast<int>(j) - 2;
            A[i * 5 + j] = std::move(std::valarray<int>({x, y}) * 16);
        }
    }
    auto f3 = [A](const std::valarray<T>& x) -> T {
        T val = 0.0;
        for (std::size_t k = 0; k < A.size(); ++k) {
            auto denominator = static_cast<T>(k);
            for (std::size_t dim = 0; dim < x.size(); ++dim) {
                denominator += std::pow(x[dim] - A[k][dim], 6);
            }
            val += 1 / denominator;
        }
        return val;
    };

    return std::make_pair(f3, bounds);
}

#endif //F3_H_
