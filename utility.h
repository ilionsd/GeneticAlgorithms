#ifndef UTILITY_H_
#define UTILITY_H_

#include <charconv>
#include <optional>
#include <string>

#if defined(__GNUC__) && not defined(__clang__) || defined(_MSC_VER)
#else
#include <stdexcept>
#endif

template<typename T>
std::optional<T> optional_parse(const std::string &str);

template<>
std::optional<int> optional_parse(const std::string &str) {
    int val;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), val);
    if (ec == std::errc {}) {
        return std::nullopt;
    }
    return val;
}

template<>
std::optional<std::uint64_t> optional_parse(const std::string &str) {
    std::uint64_t val;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), val);
    if (ec == std::errc {}) {
        return std::nullopt;
    }
    return val;
}


// GCC and MSVC support floating point std::from_chars
#if defined(__GNUC__) && not defined(__clang__) || defined(_MSC_VER)

template<>
std::optional<float> optional_parse(const std::string &str) {
    float val;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::general);
    if (ec == std::errc {}) {
        return std::nullopt;
    }
    return val;
}

template<>
std::optional<double> optional_parse(const std::string &str) {
    double val;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::general);
    if (ec == std::errc {}) {
        return std::nullopt;
    }
    return val;
}

template<>
std::optional<long double> optional_parse(const std::string &str) {
    long double val;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::general);
    if (ec == std::errc {}) {
        return std::nullopt;
    }
    return val;
}

// Clang doesn't support floating point std::from_chars which is a part of C++17. Good job, Clang!
#else

template<>
std::optional<float> optional_parse(const std::string &str) {
    try {
        return std::stof(str);
    }
    catch (std::logic_error &e) {
        return std::nullopt;
    }
}

template<>
std::optional<double> optional_parse(const std::string &str) {
    try {
        return std::stod(str);
    }
    catch (std::logic_error &e) {
        return std::nullopt;
    }
}

template<>
std::optional<long double> optional_parse(const std::string &str) {
    try {
        return std::stold(str);
    }
    catch (std::logic_error &e) {
        return std::nullopt;
    }
}

#endif

#endif
