#ifndef GENETIC_ALGORITHM__DUMMY_H_
#define GENETIC_ALGORITHM__DUMMY_H_

#include <cstddef>
#include <cstdint>
#include <valarray>
#include <vector>

#include <json/json.h>

#include "common/optimizer.h"
#include "common/maybe_set.h"
#include "common/space.h"

namespace genetic_algorithm::dummy {

enum class status: int {
    // UNEXPECTED
    UNABLE_TO_ADD_OFFSPRINGS = -2,
    IN_PROGRESS = -1,
    // SUCCESS
    CONVERGED = 0,
    REACHED_GENERATION_LIMIT = 1,
    REACHED_POPULATION_LIMIT = 2,
};

template<typename T>
struct result {
    status status;
    std::vector<std::valarray<T>> population;
    std::vector<T> values;
};

class optimizer: public common::optimizer<optimizer> {
public:
    optimizer(const std::size_t maxPopulation)
    : mMaxPopulation(maxPopulation)
    {}
    ~optimizer() = default;

    template<class Gen, typename T1, typename T2>
    auto optimize(
        Gen& generator, 
        const common::space<T1, T2>& space, 
        const std::function<T1(const std::valarray<T1>&)> fitness
    ) const -> result<T1> {
        std::vector<std::valarray<T1>> population = space.rand_n_real(generator, mMaxPopulation);
        std::vector<T1> evaluations(mMaxPopulation);
        for (std::size_t k = 0; k < mMaxPopulation; ++k) {
            evaluations[k] = fitness(population[k]);
        }
        return result<T1> { status::CONVERGED, population, evaluations };
    }

private:
    const std::size_t mMaxPopulation;
};

struct parameters {
    auto make_optimizer() const {
        return optimizer { maxPopulation };
    }

    std::size_t maxPopulation;
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

}   // namespace genetic_algorithm::dummy

namespace genetic_algorithm::common {

template<>
struct maybe_set<dummy::parameters> {
    
    dummy::parameters &params;

    maybe_set<dummy::parameters> &from_json(const Json::Value &json) {
        maybe_set<std::size_t> { 
            params.maxPopulation 
        }.from_json(json["max-population"]);
        return *this;
    }
};

}   //-- namespace genetic_algorithm::common --

#endif // GENETIC_ALGORITHM__DUMMY_H_