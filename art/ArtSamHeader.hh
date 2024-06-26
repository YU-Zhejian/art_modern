#pragma once

#include <iterator>
#include <string>
#include <utility>
#include <vector>
namespace labw {
namespace art_modern {

    class ArtSamHeader {
    public:
        ArtSamHeader(std::vector<std::string> SN,
            std::vector<int> LN)
            : SN(std::move(SN))
            , LN(std::move(LN))
        {
        }

        const std::string VN = "1.4";
        const std::string SO = "unsorted";

        std::vector<std::string> SN;
        std::vector<int> LN;

        const std::string ID = "01";
        const std::string PN = "ART_Illumina";
        std::string CL;

        std::string printHeader() const;
    };

}
}
