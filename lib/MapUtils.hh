#pragma once
#include <map>
namespace labw {
namespace art_modern {

    template <class k, class v>
    std::tuple<std::vector<k>, std::vector<v>> convert_map_to_k_v_list(const std::map<k, v>& map)
    {
        std::vector<k> keys;
        keys.reserve(map.size());
        std::vector<v> values;
        keys.reserve(values.size());
        for (auto& item : map) {
            keys.emplace_back(item.first);
            values.emplace_back(item.second);
        }
        return std::tie(keys, values);
    }

    template <class k, class v>
    std::tuple<std::vector<k>, std::vector<v>> convert_map_to_k_v_list(const std::unordered_map<k, v>& map)
    {
        std::vector<k> keys;
        keys.reserve(map.size());
        std::vector<v> values;
        keys.reserve(values.size());
        for (auto& item : map) {
            keys.emplace_back(item.first);
            values.emplace_back(item.second);
        }
        return std::tie(keys, values);
    }

}
}