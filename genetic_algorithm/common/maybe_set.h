#ifndef GENETIC_ALGORITHM__COMMON__MAYBE_SET_H_
#define GENETIC_ALGORITHM__COMMON__MAYBE_SET_H_

#include <cstddef>
#include <cstdint>
#include <memory>

#include <json/json.h>

namespace genetic_algorithm::common {

template<typename T>
struct maybe_set {};

template<>
struct maybe_set<int> {
    int& val;

    maybe_set<int> &from_json(const Json::Value &json) {
        if (json.isInt()) {
            val = json.asInt();
        }
        return *this;
    }
};

template<>
struct maybe_set<std::size_t> {
    std::size_t& val;

    maybe_set<std::size_t> &from_json(const Json::Value &json) {
        if (json.isUInt64()) {
            val = json.asUInt64();
        }
        return *this;
    }
};

template<>
struct maybe_set<std::uint64_t> {
    std::uint64_t& val;

    maybe_set<std::uint64_t> &from_json(const Json::Value &json) {
        if (json.isUInt64()) {
            val = json.asUInt64();
        }
        return *this;
    }
};

template<>
struct maybe_set<double> {
    double& val;

    maybe_set<double> &from_json(const Json::Value &json) {
        if (json.isDouble()) {
            val = json.asDouble();
        }
        return *this;
    }
};


}   //-- namespace genetic_algorithm::common --

#endif //GENETIC_ALGORITHM__COMMON__MAYBE_SET_H_
