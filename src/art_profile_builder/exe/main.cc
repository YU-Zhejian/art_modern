// TODO:
// 1. Add argument parsing
// 2. Support PE BAM/SAM

#include "libam_support/CExceptionsProxy.hh"
#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.hh"
#include "libam_support/utils/arithmetic_utils.hh"
#include "libam_support/utils/class_macros_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <absl/base/attributes.h>

#include <boost/log/trivial.hpp>

#ifdef WITH_BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace labw::art_modern;

class IntermediateEmpDistPosition {
public:
    constexpr static std::size_t ALL_IDX = 0;
    constexpr static std::size_t A_IDX = 1;
    constexpr static std::size_t C_IDX = 2;
    constexpr static std::size_t G_IDX = 3;
    constexpr static std::size_t T_IDX = 4;
    constexpr static std::size_t N_IDX = 5;
    constexpr static std::size_t BASE_IDX[] = { ALL_IDX, A_IDX, C_IDX, G_IDX, T_IDX, N_IDX };
    constexpr static std::size_t NUM_BASES = 6; // ACGTN + all
    /** ASCII to index
     *
     *  Generated using Python:
     *  @code
     *  for i in range(0, 256): print((chr(i).upper() if chr(i) in "ACGTacgt" else "N")+"_IDX", end=", ")
     *  @endcode
     */
    constexpr const static std::size_t BASE_ASCII_TO_IDX[] = { N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, A_IDX, N_IDX, C_IDX, N_IDX, N_IDX, N_IDX, G_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, T_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, A_IDX, N_IDX, C_IDX, N_IDX, N_IDX, N_IDX, G_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, T_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX,
        N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX, N_IDX };
    IntermediateEmpDistPosition()
    {
        data_.resize(NUM_BASES * (MAX_QUAL - MIN_QUAL + 1), 0);
        std::fill(data_.begin(), data_.end(), 0);
    }
    ~IntermediateEmpDistPosition() = default;
    DEFAULT_COPY(IntermediateEmpDistPosition)
    DELETE_MOVE(IntermediateEmpDistPosition)

    ABSL_ATTRIBUTE_ALWAYS_INLINE inline void add(const char base, const am_qual_t qual)
    {
        std::size_t const base_idx = BASE_ASCII_TO_IDX[static_cast<unsigned char>(base)];
        if (qual > MAX_QUAL || qual < MIN_QUAL) {
            // Ignored
            return;
        }
        auto const qual_idx = static_cast<unsigned char>(qual);
        data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx]++;
        data_[ALL_IDX * (MAX_QUAL - MIN_QUAL + 1) + qual_idx]++;
    }

    /**
     * Call this before writing.
     */
    void accumulate()
    {
        for (std::size_t base_idx = 0; base_idx < NUM_BASES; ++base_idx) {
            for (std::size_t qual_idx = MIN_QUAL + 1; qual_idx <= MAX_QUAL; ++qual_idx) {
                data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx]
                    += data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx - 1];
            }
        }
    }

    void add(IntermediateEmpDistPosition const& other)
    {
        for (std::size_t i = 0; i < NUM_BASES * (MAX_QUAL - MIN_QUAL + 1); ++i) {
            data_[i] += other.data_[i];
        }
    }

    void write(std::ostream& oss, const std::size_t pos_id, const std::size_t base_idx) const
    {
        char leading_char = 0;
        switch (base_idx) {
        case ALL_IDX:
            leading_char = '.';
            break;
        case A_IDX:
            leading_char = 'A';
            break;
        case C_IDX:
            leading_char = 'C';
            break;
        case G_IDX:
            leading_char = 'G';
            break;
        case T_IDX:
            leading_char = 'T';
            break;
        case N_IDX:
            leading_char = 'N';
            break;
        default:
            throw std::runtime_error("Impossible base index");
        }
        // .	4	2	3	4	5	6	7	8	9	10	11	12	13
        // 14	15	16	17	18	19	20	21	22	23	24	25	26	27
        // 28	29	30	31	32	33	34	35	36	37	38	39	40
        oss << leading_char << '\t' << pos_id << '\t';
        for (std::size_t qual_idx = MIN_QUAL; qual_idx <= MAX_QUAL; ++qual_idx) {
            if (data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx] != 0) {
                oss << qual_idx << '\t';
            }
        }
        oss << '\n';
        //.	4	237014	959096	2455180	4010036	4814541	5640434	6491015	7373885	9329364	10385542	11461730
        // 12674308	13887276	15197241	16550894	17990403	19487422	21239417
        // 22851768	24649672	26592306	28335198	30412329	32599170	35412897
        // 39619438	41699377	72382201	78001945	79233998	80433996	81604748
        // 82738206	83837558	84900848	85931284	86923678	87883854	136948695
        oss << leading_char << '\t' << pos_id << '\t';
        for (std::size_t qual_idx = MIN_QUAL; qual_idx <= MAX_QUAL; ++qual_idx) {
            if (data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx] != 0) {
                oss << data_[base_idx * (MAX_QUAL - MIN_QUAL + 1) + qual_idx] << '\t';
            }
        }
        oss << '\n';
    }

private:
    std::vector<std::size_t> data_;
};

class IntermediateEmpDist {
public:
    explicit IntermediateEmpDist(const std::size_t read_length)
        : read_length_(read_length)
    {
        positions_.resize(read_length);
    }
    ~IntermediateEmpDist() = default;

    IntermediateEmpDist(const IntermediateEmpDist& ied) = default;
    IntermediateEmpDist& operator=(const IntermediateEmpDist&) = delete;

    DELETE_MOVE(IntermediateEmpDist)

    bool parse_read(const bam1_t* b)
    {
        const auto* seq = bam_get_seq(b);
        const auto* qual = bam_get_qual(b);
        if (qual[0] == 0xff) {
            // No quality
            return false;
        }
        for (std::size_t i = 0; i < am_min(static_cast<std::size_t>(b->core.l_qseq), read_length_); ++i) {
            positions_[i].add(seq_nt16_str[bam_seqi(seq, i)], static_cast<am_qual_t>(qual[i]));
        }
        return true;
    }

    void accumulate()
    {
        for (auto& pos : positions_) {
            pos.accumulate();
        }
    }

    void add(const IntermediateEmpDist& other)
    {
        if (read_length_ != other.read_length_) {
            throw std::runtime_error("Incompatible IntermediateEmpDist sizes");
        }
        for (std::size_t i = 0; i < read_length_; ++i) {
            positions_[i].add(other.positions_[i]);
        }
    }

    void write(std::ostream& oss) const
    {
        for (const auto& base_idx : IntermediateEmpDistPosition::BASE_IDX) {
            for (std::size_t pos_id = 0; pos_id < positions_.size(); ++pos_id) {
                positions_[pos_id].write(oss, pos_id, base_idx);
            }
        }
    }

private:
    std::vector<IntermediateEmpDistPosition> positions_;
    const std::size_t read_length_;
};

namespace {
void view_sam(const std::string& file_path, const std::shared_ptr<IntermediateEmpDist>& ied, std::size_t thread_id = 0,
    std::size_t num_threads = 1, std::size_t num_io_threads = 4)
{
    htsThreadPool tpool = {NULL, 0};
    tpool.pool = CExceptionsProxy::assert_not_null(hts_tpool_init(num_io_threads), "HTSLib", "Failed to init HTS thread pool.");

    auto* in
        = CExceptionsProxy::assert_not_null(hts_open(file_path.c_str(), "r"), "HTSLib", "Failed to open HTS file.");
    auto* hdr = CExceptionsProxy::assert_not_null(sam_hdr_read(in), "HTSLib", "Failed to read SAM header.");
    auto* b = CExceptionsProxy::assert_not_null(bam_init1(), "HTSLib", "Failed to init BAM record.");
    am_readnum_t num_valid_reads = 0;
    am_readnum_t num_total_reads = 0;
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &tpool);

    // Skip thread_id reads
    for (std::size_t i = 0; i < thread_id; ++i) {
        if (sam_read1(in, hdr, b) < 0) {
            goto destroy;
        }
    }
    while (true) {
        // read 1 read
        if (sam_read1(in, hdr, b) < 0) {
            goto destroy;
        }
        if (ied->parse_read(b)) {
            num_valid_reads++;
        }
        num_total_reads++;
        if (num_total_reads % 500000 == 0) {
            BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": Processed "
                                    << to_si(static_cast<double>(num_total_reads), 2, 1000) << " reads, "
                                    << to_si(static_cast<double>(num_valid_reads), 2, 1000) << " ("
                                    << static_cast<double>(num_valid_reads) / static_cast<double>(num_total_reads)
                    * 100.0 << ")% valid reads.";
        }
        for (std::size_t i = 1; i < num_threads; ++i) {
            // This is really slow...
            if (sam_read1(in, hdr, b) < 0) {
                goto destroy;
            }
        }
    }
destroy:
    BOOST_LOG_TRIVIAL(info) << "Thread " << thread_id << ": Processed "
                            << to_si(static_cast<double>(num_total_reads), 2, 1000) << " reads, "
                            << to_si(static_cast<double>(num_valid_reads), 2, 1000) << " ("
                            << static_cast<double>(num_valid_reads) / static_cast<double>(num_total_reads) * 100.0
                            << ")% valid reads.";

    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    hts_close(in);
    hts_tpool_destroy(tpool.pool);
}
} // namespace

int main()
{
#ifdef WITH_BOOST_TIMER
    boost::timer::cpu_timer timer;
#endif
    const std::string file_path
        = "/home/yuzj/Documents/pbsim3_modern/explore/benchmark_other_simulators/data/soybean_SRR16074289.bam";
    const std::size_t read_length = 300;
    const std::size_t num_threads
        = 4; // Tell the users more than 4 threads will result in tremendous waste of CPU time.
    const std::size_t num_io_threads = 4;
    // Since all reads needs to be parsed num_threads times.
    htsFile* file
        = CExceptionsProxy::assert_not_null(hts_open(file_path.c_str(), "r"), "HTSLib", "Failed to open HTS file.");

    const auto* format = hts_get_format(file);
    if (format->category != sequence_data) {
        BOOST_LOG_TRIVIAL(error) << "File is not a sequence data file.";
        hts_close(file);
        return EXIT_FAILURE;
    }
    hts_close(file);

    IntermediateEmpDist ied(read_length);
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    std::vector<std::shared_ptr<IntermediateEmpDist>> ieds;
    for (std::size_t i = 0; i < num_threads; ++i) {
        auto this_ied = std::make_shared<IntermediateEmpDist>(read_length);
        ieds.emplace_back(this_ied);
        threads.emplace_back(view_sam, file_path, this_ied, i, num_threads, num_io_threads);
    }
    for (auto& t : threads) {
        t.join();
    }
    for (const auto& other_ied : ieds) {
        BOOST_LOG_TRIVIAL(info) << "Merging one intermediate result..." << std::endl;
        ied.add(*other_ied);
        std::cout << "//" << std::endl;
    }

    ied.accumulate();
    ied.write(std::cout);
    BOOST_LOG_TRIVIAL(info) << "Done.";

#ifdef WITH_BOOST_TIMER
    timer.stop();
    BOOST_LOG_TRIVIAL(info) << "Time elapsed: " << timer.format();
#endif
    // FIXME: Redo this analysis
    // Threads: 1: Time elapsed:  11.390000s wall, 11.170000s user + 0.210000s system = 11.380000s CPU (99.9%)
    // Threads: 2: Time elapsed:  7.330000s wall, 14.270000s user + 0.320000s system = 14.590000s CPU (199.0%)

    return EXIT_SUCCESS;
}
