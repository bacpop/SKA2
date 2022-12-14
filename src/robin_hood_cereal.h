//
// Created by Johanna Helene von Wachsmann on 12/10/2022.
//

#ifndef SKA_CPP_ROBIN_HOOD_CEREAL_H
#define SKA_CPP_ROBIN_HOOD_CEREAL_H

#include "robin_hood.h"
#include "ska.hpp"
#include <iostream>
#include <cereal/cereal.hpp>

namespace cereal {

//! Saving for std-like pair associative containers
template<class Archive>
inline
void save(Archive &ar, robin_hood::unordered_map<uint64_t, uint8_t> const &map) {
    ar(cereal::make_size_tag(static_cast<cereal::size_type>(map.size())));

    for (const auto &i: map)
        ar(cereal::make_map_item(i.first, i.second));
}

//! Loading for std-like pair associative containers
template<class Archive>
inline
void load(Archive &ar, robin_hood::unordered_map<uint64_t, uint8_t> &map) {
    cereal::size_type size;
    ar(cereal::make_size_tag(size));

    map.clear();

    auto hint = map.begin();
    for (size_t i = 0; i < size; ++i) {
        uint64_t key;
        uint64_t value;

        ar(cereal::make_map_item(key, value));
#ifdef CEREAL_OLDER_GCC
        hint = map.insert( hint, std::make_pair(std::move(key), std::move(value)) );
#else // NOT CEREAL_OLDER_GCC
        hint = map.emplace_hint(hint, std::move(key), std::move(value));
#endif // NOT CEREAL_OLDER_GCC
    }
}

} // namespace cereal
#endif //SKA_CPP_ROBIN_HOOD_CEREAL_H
