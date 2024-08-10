#include "misc.hh"
#include <sstream>

namespace labw {
namespace art_modern {
    std::string map_to_str(const std::map<int, char>& m)
    {
        std::ostringstream oss;
        oss << "{";

        bool first = true;
        for (const auto& pair : m) {
            if (!first) {
                oss << ", ";
            }
            first = false;
            oss << '"' << pair.first << "\": \"" << pair.second << '"';
        }

        oss << "}";
        return oss.str();
    }

    std::vector<int> range(int start, int stop, int step)
    {
        std::vector<int> result;
        for (int i = start; step > 0 ? i < stop : i > stop; i += step) {
            result.push_back(i);
        }
        return result;
    }
} // namespace art_modern
} // namespace labw
