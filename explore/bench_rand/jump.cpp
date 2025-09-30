#include "jump_matrix.h"

#include <boost/filesystem/directory.hpp>

#include <algorithm>
#include <optional>
#include <sstream>
#include <thread>

using namespace std;

enum GenType : std::uint8_t { mt, sfmt };
const std::string extension = ".bits";

template <GenType> struct GenTraits;

template <> struct GenTraits<mt> {
    using matrix_t = MT19937Matrix;
    static const size_t power2 = 0;
};

template <> struct GenTraits<sfmt> {
    using matrix_t = SFMT19937Matrix;
    static const size_t power2 = 2;
};

namespace {
template <typename Matrix> void square(const Matrix& src, Matrix& dst, std::vector<typename Matrix::buffer_t>& buffers)
{
    const auto start = std::chrono::system_clock::now();
    dst.square(src, buffers);
    const auto end = std::chrono::system_clock::now();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << (end - start).count() << "s" << std::endl;
}

std::string mk_file_name(const GenType gen_type, const std::string& path, size_t n)
{
    std::ostringstream os;
    os << path << (gen_type == mt ? "MT" : "SFMT") << "F" << std::setw(5) << std::setfill('0') << n << extension;
    return os.str();
}

template <GenType gen>
void run(const std::string& filepath, const size_t n_threads, const size_t save_frequency, const size_t stop_index)
{
    typename GenTraits<gen>::matrix_t f[2];

    size_t lastComputed = GenTraits<gen>::power2;

    for (const auto& entry : boost::filesystem::directory_iterator(filepath)) {
        if (boost::filesystem::is_regular_file(entry) && entry.path().has_extension()
            && entry.path().extension().string() == extension) {
            std::string s = entry.path().filename().string();
            s = s.substr(1, s.length() - 1 - extension.length());
            size_t const n = std::atol(s.c_str());
            lastComputed = std::max(n, lastComputed);
        }
    }

    if (lastComputed == GenTraits<gen>::power2) {
        std::cout << "initializing from base matrix: F^(2^" << GenTraits<gen>::power2 << ")\n";
    } else {
        const auto fn = mk_file_name(gen, filepath, lastComputed);
        std::cout << "initializing F^(2 ^ " << lastComputed << ") from file " << fn << "\n";
        std::ifstream is(fn, ios::binary);
        f[lastComputed % 2].fromBin(is);
    }

    std::cout << "  ";
    f[lastComputed % 2].printSparsity();

    std::vector<typename GenTraits<gen>::matrix_t::buffer_t> buffers(n_threads);

    for (size_t i = lastComputed + 1; i <= 19937 && (lastComputed < stop_index); ++i) {
        std::cout << "computing F^(2^" << i << ")\n";
        size_t const in = (i + 1) % 2;
        size_t const out = (i) % 2;
        f[out].resetZero();
        square(f[in], f[out], buffers);
        std::cout << "  ";
        f[out].printSparsity();

        if ((i % save_frequency) == 0 || i > 19930) {
            const auto fn = mk_file_name(gen, filepath, i);
            std::cout << "  saving file: " << fn << " ... ";
            std::ofstream of(fn, ios::binary);
            f[out].toBin(of);
            of.close();
            std::cout << "saved\n";
            lastComputed = i;
        }
    }
}
} // namespace
int main(int /*argc*/, const char** /*argv*/)
{
    // parse command line arguments
    const size_t n_threads = std::thread::hardware_concurrency();
    const std::string filepath = "/home/yuzj/Documents/pbsim3_modern/explore/bench_rand/dat/";
    const size_t save_frequency = 100;
    const size_t stop_index = std::numeric_limits<size_t>::max();
    std::cout << "n_threads = " << n_threads << "\n";
    std::cout << "filepath = " << filepath << "\n";
    std::cout << "save_frequency = " << save_frequency << "\n";
    std::cout << "stop_index = " << stop_index << "\n";
    run<mt>(filepath, n_threads, save_frequency, stop_index);
    run<sfmt>(filepath, n_threads, save_frequency, stop_index);
    return 0;
}
