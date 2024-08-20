#include <atomic>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/asio.hpp>
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>

#include "ArtContig.hh"
#include "art_modern_constants.hh"
#include "main_fn.hh"
#include "seq_utils.hh"

using namespace std;

namespace labw {
namespace art_modern {

    void print_banner()
    {
        BOOST_LOG_TRIVIAL(info) << "YuZJ Modified ART_Illumina";
        BOOST_LOG_TRIVIAL(info) << "Based on: v. 2008-2016, Q Version 2.5.8 (June 6, 2016)";
        BOOST_LOG_TRIVIAL(info) << "Originally written by: Weichun Huang <whduke@gmail.com>";
        BOOST_LOG_TRIVIAL(info) << "Modified by: YU Zhejian <Zhejian.23@intl.zju.edu.cn>";
    }

    void generate_all(const string& contig_name, const string& ref_seq, const ArtParams& art_params,
        const Empdist& qdist, double x_fold, const std::shared_ptr<BaseReadOutput>& output_dispatcher)
    {
        // FIXME
    }
} // namespace art_modern
} // namespace labw
