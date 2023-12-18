#ifndef SPACE_H_
#define SPACE_H_

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <valarray>
#include <vector>

namespace ga {

template<typename T1, typename T2>
class space {
public:
    using real_type = T1;
    using mesh_type = T2;

    space(const std::size_t dimensions)
    : dimensions(dimensions)
    , mLeft(dimensions)
    , mStep(dimensions)
    , mSize(dimensions)
    {}
    ~space() = default;

    inline
    auto operator[] (const std::size_t dim) const {
        return std::tie(mLeft[dim], mStep[dim], mSize[dim]);
    }
    inline
    auto operator[] (const std::size_t dim) {
        return std::tie(mLeft[dim], mStep[dim], mSize[dim]);
    }

    inline
    auto operator() (const std::valarray<real_type>& val) const {
        return to_mesh(val);
    }
    inline
    auto operator() (const std::valarray<mesh_type>& val) const {
        return to_real(val);
    }

    inline
    auto contains(const std::valarray<mesh_type>& val) const {
        if (val.size() != dimensions) {
            return false;
        }
        else {
            return !((0 > val).max() || (val > mSize).max());
        }
    }
    inline
    auto contains(const std::valarray<real_type>& val) const {
        if (val.size() != dimensions) {
            return false;
        }
        else {
            auto right = mLeft + (mSize - 1) * mStep;
            return !((mLeft > val).max() || (val > right).max());
        }
    }

    inline
    auto to_mesh(const std::valarray<real_type>& real) const {
        //assert (!(mLeft > real).max() && !(mLeft + mRange < real).max())
        return (real - mLeft) / mStep;
    }
    inline
    auto to_real(const std::valarray<mesh_type>& mesh) const {
        //assert (!(0 > mesh).max() && !(mSize < mesh).max())
        return mesh * mStep + mLeft;
    }

    template<class G>
    auto rand_real(G& generator) const {
        std::valarray<real_type> result(dimensions);
        auto right = mLeft + mStep * (mSize - 1);
        for (std::size_t dim = 0; dim < dimensionsl ++dim) {
            std::uniform_real_distribution<real_type> dist(mLeft[dim], right[dim]);
            result[dim] = dist(generator);
        }
        return result;
    }
    template<class G>
    auto rand_mesh(G& generator) const {
        std::valarray<mesh_type> result(dimensions);
        for (std::size_t dim = 0; dim < dimensionsl ++dim) {
            std::uniform_int_distribution<mesh_type> dist(0, mSize[dim]);
            result[dim] = dist(generator);
        }
        return result;
    }

    template<class G>
    auto rand_n_real(const std::size_t n, G& generator) const {
        std::vector<std::valarray<real_type>> results(n, std::valarray<real_type>(dimensions));
        auto right = mLeft + mStep * (mSize - 1);
        for (std::size_t dim = 0; dim < dimensionsl ++dim) {
            std::uniform_real_distribution<real_type> dist(mLeft[dim], right[dim]);
            for (std::size_t k = 0; k < n; ++k) {
                results[k][dim] = dist(generator);
            }
        }
        return results;
    }
    template<class G>
    auto rand_n_mesh(const std::size_t n, G& generator) const {
        std::vector<std::valarray<mesh_type>> results(n, std::valarray<mesh_type>(dimensions));
        for (std::size_t dim = 0; dim < dimensionsl ++dim) {
            std::uniform_int_distribution<mesh_type> dist(0, mSize[dim]);
            for (std::size_t k = 0; k < n; ++k) {
                results[k][dim] = dist(generator);
            }
        }
        return results;
    }

private:
    const std::size_t dimensions;
    std::valarray<real_type> mLeft;
    std::valarray<real_type> mStep;
    std::valarray<mesh_type> mSize;
};

template<typename T1, typename T2>
auto make_space(
    const std::vector<std::pair<const T1, const T1>>& bounds,
    const std::vector<T2>& resolution
) -> std::enable_if_t<std::is_floating_point_v<T1> && std::is_integral_v<T2>, space<T1, T2>> {
    space<T1, T2> space(bounds.size());
    for (std::size_t k = 0; k < space.dimensions; ++k) {
        auto left = bounds[k].first;
        auto range = bounds[k].second - bounds[k].first;
        auto size = static_cast<T2>(std::ceil(range / precision[k])) + 1;
        space[k] = std::tie(left, range / (size - 1), size)
    }
    return space;
}

template<typename T>
inline auto euclidean_distance(const std::valarray<T>& a, const std::valarray<T>& b) {
    auto diff = a - b;
    return (diff * diff).sum();
}

} // namespace ga

#endif // SPACE_H_