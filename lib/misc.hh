#pragma once

#include <iostream>
#include <map>
#include <vector>
namespace labw {
namespace art_modern {
    std::string map_to_str(const std::map<int, char>& m);

    std::vector<int> range(int start, int stop, int step = 1);
} // namespace art_modern
}; // namespace labw