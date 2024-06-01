#ifndef GENETIC_ALGORITHM__OPTIMIZER_H_
#define GENETIC_ALGORITHM__OPTIMIZER_H_

#include <functional>
#include <valarray>

#include "space.h"


namespace genetic_algorithm {

enum class status: int {
    NOT_IMPLEMENTED = -100,
    // UNEXPECTED
    UNABLE_TO_ADD_OFFSPRINGS = -2,
    IN_PROGRESS = -1,
    // SUCCESS
    CONVERGED = 0,
    REACHED_GENERATION_LIMIT = 1,
    REACHED_POPULATION_LIMIT = 2,
};

template<typename T, class Params>
struct result {
    const status status;
};

template<class OnNextGeneration, class Params>
struct optimizer {

    template<class Gen, typename T1, typename T2>
    const result<T1, Params> optimize(
        Gen& generator,
        const space<T1, T2>& space, 
        const std::function<T1(const std::valarray<T1>&)> fitness
    ) const {
        return { status::NOT_IMPLEMENTED };
    }

};



template<class Params, class... Ts>
auto make_optimizer(const Params &params, Ts... args) {
    return params.make_optimizer(args...);
}

struct max_optimizer {};
struct min_optimizer {};

} // namespace genetic_algorithm

#endif //GENETIC_ALGORITHM__OPTIMIZER_H_
