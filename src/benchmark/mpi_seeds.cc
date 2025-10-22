#include "libam_support/Constants.hh"
#include "libam_support/utils/log_utils.hh"
#include "libam_support/utils/mpi_utils.hh"
#include "libam_support/utils/rand_utils.hh"
#include "libam_support/utils/si_utils.hh"

#include <boost/log/trivial.hpp>

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace labw::art_modern;

class GeneratedSeeds {
public:
    GeneratedSeeds()
        : hostname_(mpi_hostname())
        , thread_id_(std::hash<std::thread::id> { }(std::this_thread::get_id())) { };

    /**
     * Deserialization for MPI communication
     *
     * @param ser MPI buffer
     * @param position Buffer position
     */
    GeneratedSeeds(char* ser, std::size_t& position)
        : thread_id_(0)
    {
        // Deserialize the size of the vector
        std::size_t size = 0;
        std::memcpy(&size, ser + position, sizeof(size));
        position += sizeof(size);
        // Deserialize each element
        seeds_.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            std::uint64_t seed = 0;
            std::memcpy(&seed, ser + position, sizeof(seed));
            position += sizeof(seed);
            seeds_.emplace_back(seed);
        }

        // Deserialize hostname length and content
        std::size_t hostname_length = 0;
        std::memcpy(&hostname_length, ser + position, sizeof(hostname_length));
        position += sizeof(hostname_length);

        hostname_.resize(hostname_length);
        std::memcpy(hostname_.data(), ser + position, hostname_length);
        position += hostname_length;

        // Deserialize thread id
        std::memcpy(&thread_id_, ser + position, sizeof(thread_id_));
        position += sizeof(thread_id_);
    }

    void generate_seeds_thread(const std::size_t& num_seeds)
    {
        for (std::size_t i = 0; i < num_seeds; ++i) {
            seeds_.emplace_back(rand_seed());
        }
    }

    [[nodiscard]] std::size_t ser_size() const
    {
        std::size_t size = sizeof(std::size_t); // for vector size
        size += seeds_.size() * sizeof(std::uint64_t); // for each seed
        size += sizeof(std::size_t); // for hostname length
        size += hostname_.size() * sizeof(char); // for hostname content
        size += sizeof(std::size_t); // for thread id hash
        return size;
    }

    /**
     * Serialization for MPI communication
     */
    void ser(char* buffer, std::size_t& position) const
    {
        // Serialize the size of the vector
        const std::size_t size = seeds_.size();
        std::memcpy(buffer + position, &size, sizeof(size));
        position += sizeof(size);

        // Serialize each element
        for (const auto& seed : seeds_) {
            std::memcpy(buffer + position, &seed, sizeof(seed));
            position += sizeof(seed);
        }

        // Serialize hostname length and content
        const std::size_t hostname_length = hostname_.size();
        std::memcpy(buffer + position, &hostname_length, sizeof(hostname_length));
        position += sizeof(hostname_length);

        std::memcpy(buffer + position, hostname_.data(), hostname_length);
        position += hostname_length;

        // Serialize thread id
        std::memcpy(buffer + position, &thread_id_, sizeof(thread_id_));
        position += sizeof(thread_id_);
    }

    std::vector<std::uint64_t> seeds_;

private:
    std::string hostname_;
    std::uint64_t thread_id_;
};

namespace {
constexpr int NUM_THREADS = 4;
constexpr int NUM_SEEDS = 1024;
void main_worker()
{
    std::vector<std::shared_ptr<GeneratedSeeds>> gss;
    std::vector<std::thread> threads;
    // Generate seeds
    threads.reserve(NUM_THREADS);
    for (std::size_t thread_idx = 0; thread_idx < NUM_THREADS; ++thread_idx) {
        auto gs = std::make_shared<GeneratedSeeds>();
        gss.emplace_back(gs);
        threads.emplace_back([thread_idx, &gss]() {
            BOOST_LOG_TRIVIAL(info) << "Rank " << mpi_rank() << " thread=" << thread_idx << " generating seeds.";
            gss[thread_idx]->generate_seeds_thread(NUM_SEEDS);
            BOOST_LOG_TRIVIAL(info) << "Rank " << mpi_rank() << " thread=" << thread_idx << " generated "
                                    << to_si(gss[thread_idx]->seeds_.size()) << " seeds.";
        });
    }
    for (auto& t : threads) {
        t.join();
    }
    for (std::size_t thread_idx = 0; thread_idx < NUM_THREADS; ++thread_idx) {
        const auto& gs = gss.at(thread_idx);

        // Serialize
        const std::size_t ser_size = gs->ser_size();
        char* buffer = static_cast<char*>(std::calloc(ser_size, sizeof(char)));
        std::size_t position = 0;
        gs->ser(buffer, position);
        if (ser_size != position) {
            BOOST_LOG_TRIVIAL(fatal) << "Rank " << mpi_rank() << " thread=" << thread_idx << " serialized seeds size ("
                                     << ser_size << ") != position (" << to_si(position) << ")!";
            abort_mpi();
        }
        BOOST_LOG_TRIVIAL(info) << "Rank " << mpi_rank() << " thread=" << thread_idx << " serialized "
                                << to_si(ser_size) << " bytes.";

        // Send via MPI
        BOOST_LOG_TRIVIAL(info) << "Rank " << mpi_rank() << " thread=" << thread_idx << " sending via MPI.";
        MPI_Send(
            buffer, static_cast<int>(ser_size), MPI_CHAR, MPI_MAIN_RANK, static_cast<int>(thread_idx), MPI_COMM_WORLD);
        BOOST_LOG_TRIVIAL(info) << "Rank " << mpi_rank() << " thread=" << thread_idx << " SEND SUCCESS.";

        std::free(buffer);
    }
    BOOST_LOG_TRIVIAL(info) << "Rank " << mpi_rank() << " all threads joined.";
}
void main_manager()
{
    std::vector<GeneratedSeeds> all_seeds;
    // Receive and deserialize via MPI
    for (int source_rank = MPI_MAIN_RANK + 1; source_rank < mpi_size(); ++source_rank) {
        for (std::size_t thread_idx = 0; thread_idx < NUM_THREADS; ++thread_idx) {
            BOOST_LOG_TRIVIAL(info) << "Rank " << source_rank << " thread=" << thread_idx << " waiting to RECV.";
            MPI_Status status;
            // Probe for an incoming message from process zero
            MPI_Probe(source_rank, thread_idx, MPI_COMM_WORLD, &status);
            BOOST_LOG_TRIVIAL(info) << "Rank " << source_rank << " thread=" << thread_idx << " PROBE SUCCESS.";

            int ser_size = 0;
            // When probe returns, the status object has the size and other
            // attributes of the incoming message. Get the message size
            MPI_Get_count(&status, MPI_CHAR, &ser_size);
            BOOST_LOG_TRIVIAL(info) << "Rank " << source_rank << " thread=" << thread_idx << " sized "
                                    << to_si(ser_size) << ".";

            // Then receive the actual serialized data
            char* buffer = new char[ser_size];
            MPI_Recv(buffer, ser_size, MPI_CHAR, source_rank, thread_idx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            BOOST_LOG_TRIVIAL(info) << "Rank " << source_rank << " thread=" << thread_idx << " RECV SUCCESS.";

            // Deserialize
            std::size_t position = 0;
            GeneratedSeeds gs(buffer, position);
            if (ser_size != position) {
                BOOST_LOG_TRIVIAL(fatal) << "Rank " << source_rank << " thread=" << thread_idx
                                         << " serialized seeds size (" << ser_size << ") != position ("
                                         << to_si(position) << ")!";
                abort_mpi();
            }
            all_seeds.emplace_back(std::move(gs));
            BOOST_LOG_TRIVIAL(info) << "Rank " << source_rank << " thread=" << thread_idx << " DECODE SUCCESS.";
            delete[] buffer;
        }
    }
    // Find number of collision seeds
    std::set<std::uint64_t> unique_seeds;
    std::size_t collision_count = 0;
    std::size_t total_seeds = 0;
    for (const auto& gs : all_seeds) {
        total_seeds += gs.seeds_.size();
        for (const auto& seed : gs.seeds_) {
            if (!unique_seeds.insert(seed).second) {
                ++collision_count;
            }
        }
    }
    BOOST_LOG_TRIVIAL(info) << "Total seeds generated: " << total_seeds;
    BOOST_LOG_TRIVIAL(info) << "Total collision seeds: " << collision_count;
}
} // namespace

int main(int argc, char** argv)
{
    init_mpi(&argc, &argv);
    init_logger();
    init_file_logger(true);

    if (mpi_rank() != MPI_MAIN_RANK) {
        main_worker();
    } else {
        main_manager();
    }
    exit_mpi();
    return EXIT_SUCCESS;
}
