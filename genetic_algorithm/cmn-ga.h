#ifndef GENETIC_ALGORITHM__CMN_GA_H_
#define GENETIC_ALGORITHM__CMN_GA_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <iostream>
#include <iterator>
#include <functional>
#include <random>
#include <tuple>
#include <valarray>
#include <vector>

#include "space.h"

namespace genetic_algorithm {

template<class Generator>
class cmn_optimizer {
public:
    cmn_optimizer(
        const std::size_t initialPopulation,
        const std::size_t maxPopulation,
	    const std::size_t crossoversPerGeneration,
	    const std::size_t mutationsPerGeneration,
	    const int nMin,
	    const int nCrowd,
	    const int nTry,
	    const double targetClosestAvg,
        std::random_device& rd) 
    : mInitialPopulation(initialPopulation)
    , mMaxPopulation(maxPopulation)
    , mCrossoversPerGeneration(crossoversPerGeneration)
    , mMutationsPerGeneration(mutationsPerGeneration)
    , mMin(nMin)
    , mCrowd(nCrowd)
    , mTry(nTry)
    , mGenerator(rd)
    {}
    ~cmn_optimizer() = default;

    template<typename T1, typename T2, class OnNextGeneration>
    auto optimize(const space<T1, T2>& space, const std::function<T1(const std::valarray<T1>&)> fitness, OnNextGeneration onNextGeneration) {
        // used to allocate memory
        const std::size_t increaseCapacity = (mInitialPopulation > mCrossoversPerGeneration + mMutationsPerGeneration) 
            ? mInitialPopulation
            : mCrossoversPerGeneration + mMutationsPerGeneration;
        const std::size_t populationCapacity = ((mInitialPopulation > mMaxPopulation)
            ? mInitialPopulation
            : mMaxPopulation) 
            + mCrossoversPerGeneration + mMutationsPerGeneration;
        
        std::size_t memoryUsage = 0;
        // Population mesh+real
        memoryUsage += sizeof(std::vector<std::valarray<T1>>) 
            + populationCapacity * (sizeof(std::valarray<T1>) + space.dimensions * sizeof(T1));
        memoryUsage += sizeof(std::vector<std::valarray<T2>>) 
            + populationCapacity * (sizeof(std::valarray<T2>) + space.dimensions * sizeof(T2));
        // Population increase mesh+real
        memoryUsage += sizeof(std::vector<std::valarray<T1>>) 
            + increaseCapacity * (sizeof(std::valarray<T1>) + space.dimensions * sizeof(T1));
        memoryUsage += sizeof(std::vector<std::valarray<T2>>) 
            + increaseCapacity * (sizeof(std::valarray<T2>) + space.dimensions * sizeof(T2));
        // Fitness
        memoryUsage += 6 * (sizeof(std::vector<T1>) + populationCapacity * sizeof(T1));
        // Utility
        memoryUsage += sizeof(std::vector<double>) + populationCapacity * sizeof(double);
        memoryUsage += sizeof(std::vector<bool>) 
            + static_cast<std::size_t>(std::ceil(populationCapacity / 8.0)); // estimating 1 bit per bool and 8 bits per byte
        memoryUsage += sizeof(std::vector<std::pair<std::size_t, std::size_t>>) 
            + populationCapacity * sizeof(std::pair<std::size_t, std::size_t>);

        const std::size_t memoryUsagePeak = memoryUsage;

        // tracks current population size
        std::size_t populationSize = 0;
        std::size_t increaseSize = mInitialPopulation;

        // stores individuals in mesh space
        std::vector<std::valarray<T2>> 
            populationMesh(populationCapacity),
            increaseMesh = space.rand_n_mesh(mInitialPopulation, mGenerator);
        std::vector<std::valarray<T1>>
            populationReal(populationCapacity),
            increaseReal(increaseCapacity);
        std::vector<T1> 
            // destructible copy of data for std::nth_element
            destructibleDataCopy(populationCapacity),
            temporaryDataCopy(populationCapacity),
            // stores fitness function evaluations
            fitnessValues(populationCapacity),
            fitnessScaled(populationCapacity),
            fitnessWeighted(populationCapacity),
            proximityWeight(populationCapacity);
        auto fitnessValueMin = std::numeric_limits<T1>::max();
        
        std::vector<bool> mask(populationCapacity);
        std::vector<double> probability(populationCapacity);
        std::vector<std::pair<std::size_t, std::size_t>> crossover(mCrossoversPerGeneration);
        
        for (std::size_t k = 0; k < increaseSize; ++k) {
            increaseReal[k] = space(increaseMesh[k]);
        }

        mGeneration = 0;

        while (increaseSize + populationSize < mMaxPopulation) {
            //-- Adding new individuals to population --
            for (std::size_t k = 0; k < increase.size(); ++k) {

                auto fitnessValue = fitness(increaseReal[k]);
                if (fitnessValue < fitnessValueMin) {
                    fitnessValueMin = fitnessValue;
                }
                std::tie(fitnessValues[populationSize + k], populationMesh[populationSize + k], populationReal[populationSize + k]) = std::move(std::tie(fitnessValue, increaseMesh[k], increaseReal[k]));
            }
            populationSize += increaseSize;
            increaseSize = 0;

            onNextGeneration(
                const_cast<const std::size_t>(mGeneration), 
                const_cast<const std::size_t>(populationSize),
                const_cast<const std::vector<std::valarray<T1>>&>(populationReal),
                const_cast<const std::vector<T1>&>(fitnessValues)
            );
            ++mGeneration;

            //-- Finding local optima --
            {
                auto begin = mask.begin();
                auto end = mask.end();
                std::fill(begin, end, true);
            }
            for (std::size_t i = 0; i < populationSize; ++i) {
                if (!mask[i]) {
                    continue;
                }
                for (std::size_t j = 0; j < populationSize; ++j) {
                    temporaryDataCopy[j] = destructibleDataCopy[j] = euclidean_distance(populationReal[i], populationReal[j]);
                }
                {
                    auto begin = destructibleDataCopy.begin();
                    auto end = destructibleDataCopy.end();
                    std::nth_element(begin, std::advance(destructibleDataCopy.begin(), mMin), end);
                }
                {
                    auto locality = destructibleDataCopy[mMin];
                    for (std::size_t j = 0; mask[i] && j < populationSize; ++j) {
                        auto isNonLocal = (i == j) || (temporaryDataCopy[j] > locality);
                        mask[i] &&= !isNonLocal || (fitnessValues[i] > fitnessValues[j]);
                        mask[j] &&= isNonLocal || (fitnessValues[i] < fitnessValues[j]);
                    }
                }
            }

            //-- Stopping criteria --
            auto closestAvg = 0.0;
            {
                std::size_t localOptimaNumber = 0;
                auto accumulator = 0.0;
                for (std::size_t i = 0; i < populationSize; ++i) {
                    if (mask[i]) {
                        ++localOptimaNumber;
                        auto closest = std::numeric_limits<double>::max();
                        for (std::size_t j = 0; j < populationSize; ++j) {
                            auto distance = euclidean_distance(populationMesh[i], populationMesh[j]);
                            if (i != j && closest > distance) {
                                closest = distance;
                            }
                        }
                        accumulator += closest;
                    }
                }
                closestAvg = accumulator / localOptimaNumber;
            }
            if (closestAvg < mTargetClosestAvg) {
                break;
            }

            //-- Calculating adjusted fitness --
            {
                auto begin = proximityWeight.begin();
                auto end = std::advance(proximityWeight.begin(), populationSize);
                std::fill(begin, end, static_cast<T1>(0.0));
            }
            {
                auto begin = fitnessWeighted.begin();
                auto end = std::advance(fitnessWeighted.begin(), populationSize);
                std::fill(begin, end, static_cast<T1>(0.0));
            }
            for (std::size_t j = 0; j < populationSize; ++j) {
                if (!mask[j]) {
                    continue;
                }
                {
                    // Calculating linear-scaled fitness per optima
                    auto begin = fitnessValues.cbegin();
                    auto end = std::advance(fitnessValues.cbegin(), populationSize);
                    auto fitnessOptimusDiff = fitnessValues[j] - fitnessValueMin;
                    std::transform(begin, end, fitnessScaled.begin(), [fitnessValueMin, fitnessOptimusDiff] (const auto f) {
                        return (f - fitnessValueMin) / fitnessOptimusDiff;
                    });
                }
                {
                    // Making a copy of linear-scaled fitness, because nth_element rearranges elements
                    auto begin = fitnessScaled.cbegin();
                    auto end = std::advance(fitnessScaled.cbegin(), populationSize);
                    std::copy(begin, end, destructibleDataCopy.begin());
                }
                {
                    // Finding median of linear-scaled fitness
                    auto begin = destructibleDataCopy.begin();
                    auto end = std::advance(destructibleDataCopy.begin(), populationSize);
                    std::nth_element(begin, std::advance(begin, populationSize / 2), end);
                }
                {
                    // Calculating exponentially-scaled fitness per optima
                    auto power =  std::log(0.5) / std::log(destructibleDataCopy[populationSize / 2]);
                    auto begin = fitnessScaled.cbegin();
                    auto end = std::advance(fitnessScaled.cbegin(), populationSize);
                    std::transform(begin, end, fitnessScaled.begin(), [power] (const auto f) {
                        return std::pow(f, power);
                    });
                }
                // Combining scaled fitness per optima into a scaled fitness overall
                for (std::size_t i = 0; i < populationSize; ++i) {
                    auto proximity = static_cast<T1>(1.0) / euclidean_distance(populationReal[j], populationReal[i]);
                    proximityWeight[i] += proximity;
                    fitnessWeighted[i] += proximity * fitnessScaled[i];
                }
            }
            {
                // Normalizing overall scaled fitness by proximity
                auto begin = proximityWeight.cbegin();
                auto end = std::advance(proximityWeight.cbegin(), populationSize);
                std::transform(begin, end, fitnessWeighted.cbegin(), fitnessWeighted.begin(), [](const auto p, const auto pf) {
                    auto wf = pf / p;
                    return std::isnan(wf) ? static_cast<T1>(1.0) : wf;
                });
            }

            //-- Selection --
            double accumulator = 0.0;
            {
                auto begin = fitnessWeighted.cbegin();
                auto end = std::advance(fitnessWeighted.cbegin(), populationSize);
                accumulator = std::accumulate(begin, end, accumulator);
            }
            {
                auto begin = fitnessWeighted.cbegin();
                auto end = std::advance(fitnessWeighted.cbegin(), populationSize);
                std::transform(begin, end, probability.begin(), [accumulator](const auto w) {
                    return w / accumulator;
                });
            }
            {
                // Fitness-Proportionate selection (FPS) for P1
                auto begin = probability.cbegin();
                auto end = std::advance(probability.cbegin(), populationSize);
                auto p1 = std::discrete_distribution<std::size_t>(begin, end);
                for (std::size_t k = 0; k < mCrossoversPerGeneration; ++k) {
                    crossover[k].first = p1(mGenerator);
                }
            }
            for (std::size_t k = 0; k < mCrossoversPerGeneration; ++k) {
                // Proximity-Proportionate selection (PPS) for P2
                for (std::size_t i = 0; i < populationSize; ++i) {
                    probability[i] = (i == crossover[k].first)
                        ? 0.0
                        : 1.0 / euclidean_distance(populationReal[crossover[k].first], populationReal[i]);
                }
                auto begin = probability.cbegin();
                auto end = std::advance(probability.cbegin(), populationSize);
                auto p2 = std::discrete_distribution<std::size_t>(begin, end);
                auto fitnessMax = std::numeric_limits<T1>::min();
                for (std::size_t crowd = 0; crowd < mCrowd; ++crowd) {
                    auto candidate = p2(mGenerator);
                    if (fitnessMax < fitnessWeighted[candidate]) {
                        std::tie(crossover[k].second, fitnessMax) = std::tie(candidate, fitnessWeighted[candidate]);
                    }
                }
            }
            {
                // Flags for crossover pairs used
                auto begin = mask.begin();
                auto end = std::advance(mask.begin(), mCrossoversPerGeneration);
                std::fill(begin, end, false);
                // Crossover
                std::size_t remainingTries = mTry;
                std::size_t k = 0;
                while (remainingTries > 0 && increaseSize < mCrossoversPerGeneration) {
                    auto pair = k++ % mCrossoversPerGeneration;
                    if (!mask[pair]) {
                        --remainingTries;
                        auto [p1, p2] = crossover[pair];
                        auto offspringMesh = offspring_by_crossover(populationMesh[p1], populationMesh[p2]);
                        auto offspringReal = space(offspringMesh);
                        std::size_t closestIndex = 0;
                        for (std::size_t i = 1; i < populationSize; ++i) {
                            auto distance = euclidean_distance(offspringReal, populationReal[i]);
                            if (minDistance > distance) {
                                std::tie(minDistance, closestIndex) = std::tie(distance, i);
                            }
                        }
                        if (minDistance > r_min(fitnessWeighted[closestIndex], mGeneration)) {
                            std::tie(increaseMesh[increaseSize], increaseReal[increaseSize]) = std::move(std::tie(offspringMesh, offspringReal));
                            mask[pair] = true;
                            ++increaseSize;
                        }
                    }
                }
                //-- Mutation --
                auto d = std::uniform_int_distribution<std::size_t>(0, populationSize);
                while (remainingTries > 0 && increaseSize < mCrossoversPerGeneration + mMutationsPerGeneration) {
                    --remainingTries;
                    auto k = d(mGenerator);
                    auto offspringMesh = offspring_by_mutation(populationMesh[k]);
                    if (space.contains(offspringMesh)) {
                        auto offspringReal = space(offspringMesh);
                        std::size_t closestIndex = 0;
                        for (std::size_t i = 1; i < populationSize; ++i) {
                            auto distance = euclidean_distance(offspringReal, populationReal[i]);
                            if (minDistance > distance) {
                                std::tie(minDistance, closestIndex) = std::tie(distance, i);
                            }
                        }
                        if (minDistance > r_min(fitnessWeighted[closestIndex], mGeneration)) {
                            std::tie(increaseMesh[increaseSize], increaseReal[increaseSize]) = std::move(std::tie(offspringMesh, offspringReal));
                            ++increaseSize;
                        }
                    }
                }
            }
        } // while

        populationReal.resize(populationSize);
        fitnessValues.resize(populationSize);
        return std::make_tuple(populationReal, fitnessValues);
    }

protected:

    template<typename T>
    auto r_min(const T fitness, const std::uint64_t generation) const {
        double factorG = 1.0 - 0.5 * std::pow(0.9, generation);
        return 0.08 * (1.001 - factorG * fitness);
    }

    template<typename T>
    auto offspring_by_crossover(const std::valarray<T>& p1, const std::valarray<T>& p2) const {
        auto dimensions = std::min(p1.size(), p2.size());
        auto offspring = std::valarray<T>(dimensions);
        for (std::size_t dim = 0; dim < dimensions; ++dim) {
            if (p1[dim] == p2[dim]) {
                offspring[dim] = p1[dim];
            }
            else {
                auto [min, max] = std::minmax(p1[dim], p2[dim]);
                auto d = std::uniform_int_distribution<T>(min, max);
                offspring[dim] = d(mGenerator);
            }
        }
        return offspring;
    }

    template<typename T>
    auto offspring_by_mutation(const std::valarray<T>& p1) const {
        auto dimensions = p1.size();
        auto stddev = 0.33 * dimensions;
        auto offspring = std::valarray<T>(dimensions);
        for (std::size_t dim = 0; dim < dimensions; ++dim) {
            auto d = std::normal_distribution<double>(p1[dim], stddev);
            offspring[dim] = static_cast<T>(std::round(d(mGenerator)));
        }
        return offspring;
    }

private:
    const std::size_t mInitialPopulation;
    const std::size_t mMaxPopulation;
	const std::size_t mCrossoversPerGeneration;
	const std::size_t mMutationsPerGeneration;
	const int mMin;
	const int mCrowd;
	const int mTry;
	const double mTargetClosestAvg;

    Generator mGenerator;
    std::uint64_t mGeneration;
};

struct cmn_parameters {
    template<class Generator>
    auto make_optimizer(std::random_device& randomDevice) const {
        return cmn_optimizer<Generator>(
            initialPopulation,
            maxPopulation,
            crossoversPerGeneration,
            mutationsPerGeneration,
            nMin, nCrowd, nTry,
            targetClosestAvg,
            randomDevice
        );
    }

    std::size_t initialPopulation;
    std::size_t maxPopulation;
	std::size_t crossoversPerGeneration;
	std::size_t mutationsPerGeneration;
	int nMin;
	int nCrowd;
	int nTry;
	
	double targetClosestAvg;
};

template<typename CharT>
std::basic_istream<CharT>& operator<<(std::basic_istream<CharT>& is, cmn_parameters &params) {
    is >> params.initialPopulation >> params.crossoversPerGeneration >> params.mutationsPerGeneration;
    is >> params.nMin >> params.nCrowd >> params.nTry;
    is >> params.maxPopulation;
    return is;
}

template<typename CharT>
std::basic_ostream<CharT>& operator>>(std::basic_ostream<CharT>& os, cmn_parameters &params) {
    return os;
}

}   //-- namespace genetic_algorithm --

#endif // GENETIC_ALGORITHM__CMN_GA_H_
