#ifndef GENETIC_ALGORITHM__PARAMETERS_H_
#define GENETIC_ALGORITHM__PARAMETERS_H_

#include <string>

#include <json/json.h>

namespace genetic_algorithm {

template<class Params>
struct model {};


template<typename T>
class entry {
public:

    entry(const std::string &name)
    : name(name)
    , value()
    {}
    entry(const std::string &name, const T &value)
    : name(name)
    , value(value)
    {}

    bool load_from(const Json::Value &json);

private:
    std::string name;
    T value;
};

template<>
bool entry<int>::load_from(const Json::Value &json) {
    auto jsonVal = json[name];
    if (jsonVal.isInt()) {
        value = jsonVal.asInt();
        return true;
    }
    return false;
}

template<>
bool entry<std::size_t>::load_from(const Json::Value &json) {
    auto jsonVal = json[name];
    if (jsonVal.isUInt64()) {
        value = jsonVal.asUInt64();
        return true;
    }
    return false;
}

template<>
bool entry<std::uint64_t>::load_from(const Json::Value &json) {
    auto jsonVal = json[name];
    if (jsonVal.isUInt64()) {
        value = jsonVal.asUInt64();
        return true;
    }
    return false;
}

template<>
bool entry<float>::load_from(const Json::Value &json) {
    auto jsonVal = json[name];
    if (jsonVal.isNumeric()) {
        value = jsonVal.asFloat();
        return true;
    }
    return false;
}

template<>
bool entry<double>::load_from(const Json::Value &json) {
    auto jsonVal = json[name];
    if (jsonVal.isNumeric()) {
        value = jsonVal.asDouble();
        return true;
    }
    return false;
}

}

#endif //GENETIC_ALGORITHM__PARAMETERS_H_
