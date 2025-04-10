#ifndef GENETIC_ALGORITHM__SPACE_H_
#define GENETIC_ALGORITHM__SPACE_H_

#include <cstddef>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <valarray>
#include <vector>

namespace genetic_algorithm {

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
    , mTemp(dimensions)
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
    auto to_mesh(const std::valarray<real_type>& real) const -> std::valarray<mesh_type> {
        //assert (!(mLeft > real).max() && !(mLeft + mRange < real).max())
        mTemp = (real - mLeft) / mStep;
        auto mesh = std::valarray<mesh_type>(dimensions);
        std::transform(std::begin(mTemp), std::end(mTemp), std::begin(mesh), [](auto val) {
            return static_cast<mesh_type>(std::round(val));
        });
        return mesh;
    }
    inline
    auto to_real(const std::valarray<mesh_type>& mesh) const -> std::valarray<real_type> {
        //assert (!(0 > mesh).max() && !(mSize < mesh).max())
        std::transform(std::begin(mesh), std::end(mesh), std::begin(mTemp), [](auto val) {
            return static_cast<real_type>(val);
        });

        return std::valarray<real_type>(mTemp * mStep + mLeft);
    }

    template<class G>
    auto rand_real(G& generator) const {
        std::valarray<real_type> result(dimensions);
        auto right = mLeft + mStep * (mSize - 1);
        for (std::size_t dim = 0; dim < dimensions; ++dim) {
            std::uniform_real_distribution<real_type> dist(mLeft[dim], right[dim]);
            result[dim] = dist(generator);
        }
        return result;
    }
    template<class G>
    auto rand_mesh(G& generator) const {
        std::valarray<mesh_type> result(dimensions);
        for (std::size_t dim = 0; dim < dimensions; ++dim) {
            std::uniform_int_distribution<mesh_type> dist(0, mSize[dim]);
            result[dim] = dist(generator);
        }
        return result;
    }

    template<class G>
    auto rand_n_real(G& generator, const std::size_t n) const {
        std::transform(std::begin(mSize), std::end(mSize), std::begin(mTemp), [](auto val) {
            return static_cast<real_type>(val - 1);
        });
        std::vector<std::valarray<real_type>> results(n, std::valarray<real_type>(dimensions));
        auto right = mLeft + mStep * mTemp;
        for (std::size_t dim = 0; dim < dimensions; ++dim) {
            std::uniform_real_distribution<real_type> dist(mLeft[dim], right[dim]);
            for (std::size_t k = 0; k < n; ++k) {
                results[k][dim] = dist(generator);
            }
        }
        return results;
    }
    template<class G>
    auto rand_n_mesh(G& generator, const std::size_t n) const {
        std::vector<std::valarray<mesh_type>> results(n, std::valarray<mesh_type>(dimensions));
        for (std::size_t dim = 0; dim < dimensions; ++dim) {
            std::uniform_int_distribution<mesh_type> dist(0, mSize[dim]);
            for (std::size_t k = 0; k < n; ++k) {
                results[k][dim] = dist(generator);
            }
        }
        return results;
    }
    
    const std::size_t dimensions;

    inline const auto& sizes() const {
        return mSize;
    }

private:
    std::valarray<real_type> mLeft;
    std::valarray<real_type> mStep;
    std::valarray<mesh_type> mSize;
    mutable std::valarray<real_type> mTemp;
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
        auto step = static_cast<T1>(range / (resolution[k] - 1));
        space[k] = std::tie(left, step, resolution[k]);
    }
    return space;
}

template<typename CharT, typename T1, typename T2>
std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os, const space<T1, T2> &space) {
    return os;
}

template<typename T>
inline auto euclidean_distance(const std::valarray<T>& a, const std::valarray<T>& b) {
    auto diff = a - b;
    auto underroot = (diff * diff).sum();
    return std::sqrt(underroot);
}

} // namespace genetic_algorithm

#endif // GENETIC_ALGORITHM__SPACE_H_
