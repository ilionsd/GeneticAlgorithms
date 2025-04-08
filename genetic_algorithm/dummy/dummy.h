#ifndef GENETIC_ALGORITHM__DUMMY__DUMMY_H_
#define GENETIC_ALGORITHM__DUMMY__DUMMY_H_

#include <cstddef>
#include <cstdint>
#include <valarray>
#include <vector>

#include <json/json.h>

#include "../optimizer.h"
#include "../maybe_set.h"
#include "../space.h"

namespace genetic_algorithm {

namespace dummy {

struct parameters {

    template<class OnNextGeneration>
    auto make_optimizer(OnNextGeneration onNextGeneration) const {
        return optimizer<OnNextGeneration, parameters> { maxPopulation };
    }

    std::uint64_t maxPopulation;
};

template<typename CharT>
std::basic_istream<CharT>& operator>>(std::basic_istream<CharT> &is, parameters &params) {
    is >> params.maxPopulation;
    return is;
}

template<typename CharT>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, parameters &params) {
    os << params.maxPopulation;
    return os;
}

} //-- namespace dummy --


template<typename T>
struct result<T, dummy::parameters> {
    status status;
    std::uint64_t size;
    std::vector<std::valarray<T>> population;
    std::vector<T> values;
};


template<typename CharT, typename T>
inline std::basic_ostream<CharT>& operator<< (std::basic_ostream<CharT>& os, const result<T, dummy::parameters>& res) {
    os << "\nDummy" << std::endl;
    return os;
}


template<class OnNextGeneration>
class optimizer<OnNextGeneration, dummy::parameters> {
public:
    optimizer(const std::size_t maxPopulation)
    : mMaxPopulation(maxPopulation)
    {}
    ~optimizer() = default;

    template<class Gen, typename T1, typename T2>
    auto optimize(
        Gen& generator, 
        const space<T1, T2>& space, 
        const std::function<T1(const std::valarray<T1>&)> fitness
    ) const -> const result<T1, dummy::parameters> {
        std::vector<std::valarray<T1>> population = space.rand_n_real(generator, mMaxPopulation);
        std::vector<T1> evaluations(mMaxPopulation);
        for (std::size_t k = 0; k < mMaxPopulation; ++k) {
            evaluations[k] = fitness(population[k]);
        }
        return { status::CONVERGED, mMaxPopulation, population, evaluations };
    }

private:
    const std::size_t mMaxPopulation;
};


template<>
struct maybe_set<dummy::parameters> {

    dummy::parameters& params;

    maybe_set<dummy::parameters>& from_json(const Json::Value& json) {
        maybe_set<std::uint64_t> {
            params.maxPopulation
        }.from_json(json["max-population"]);
        return *this;
    }
};

}   // namespace genetic_algorithm

#endif // GENETIC_ALGORITHM__DUMMY__DUMMY_H_
