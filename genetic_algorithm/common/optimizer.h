#ifndef GENETIC_ALGORITHM__COMMON__OPTIMIZER_H_
#define GENETIC_ALGORITHM__COMMON__OPTIMIZER_H_

#include <functional>
#include <valarray>

#include "space.h"


namespace genetic_algorithm::common {

template<class Derived>
struct optimizer {

    template<class Gen, typename T1, typename T2>
    auto optimize(
        Gen& generator, 
        const space<T1, T2>& space, 
        const std::function<T1(const std::valarray<T1>&)> fitness
    ) const {
        return static_cast<Derived *>(this)->optimize(generator, space, fitness);
    }

};

struct max_optimizer {};
struct min_optimizer {};

} // namespace genetic_algorithm::common

#endif //GENETIC_ALGORITHM__COMMON__OPTIMIZER_H_
