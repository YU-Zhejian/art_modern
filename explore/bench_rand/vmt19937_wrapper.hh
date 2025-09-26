#pragma once
#include "class_utils.hh"

#include <VRandGen.h>

#include <climits>

namespace {
std::string get_filename()
{
    switch (SIMD_N_BITS) {
    case 128:
        return "F19935.bits";
    case 256:
        return "F19934.bits";
    case 512:
        return "F19933.bits";
    default:
        throw std::runtime_error("Unsupported SIMD_N_BITS");
    }
}
} // namespace

constexpr uint32_t SEED_LENGTH = 4;
constexpr uint32_t SEED_INIT[SEED_LENGTH] = { 0x123, 0x234, 0x345, 0x456 };

class VMT19937RandomDeviceBase {
public:
    using result_type = unsigned int;

    VMT19937RandomDeviceBase() = default;
    ~VMT19937RandomDeviceBase() = default;

    DELETE_COPY_MOVE(VMT19937RandomDeviceBase)
};

class VMT19937RandomDevice : VMT19937RandomDeviceBase {
public:
    DELETE_COPY_MOVE(VMT19937RandomDevice)

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return UINT_MAX; }

    VMT19937RandomDevice()
    {
        MT19937Matrix const* jumpMatrix = new MT19937Matrix(
            "/home/yuzj/Documents/pbsim3_modern/deps/VMT19937-master/dat/mt/dat2/" + get_filename());
        rand_stream_ = new VMT19937<SIMD_N_BITS, QM_Scalar>(SEED_INIT, SEED_LENGTH, 0, nullptr, jumpMatrix);
        delete jumpMatrix;
    }
    ~VMT19937RandomDevice() { delete rand_stream_; };

    result_type operator()() { return rand_stream_->genrand_uint32(); }

private:
    VMT19937<SIMD_N_BITS, QM_Scalar>* rand_stream_;
};

class VSFMT19937RandomDevice : VMT19937RandomDeviceBase {
public:
    DELETE_COPY_MOVE(VSFMT19937RandomDevice)

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return UINT_MAX; }

    VSFMT19937RandomDevice()
    {
        SFMT19937Matrix const* jumpMatrix
            = new SFMT19937Matrix("/home/yuzj/Documents/pbsim3_modern/deps/VMT19937-master/dat/sfmt/" + get_filename());
        rand_stream_ = new VSFMT19937<SIMD_N_BITS, QM_Scalar>(SEED_INIT, SEED_LENGTH, 0, nullptr, jumpMatrix);
        delete jumpMatrix;
    }
    ~VSFMT19937RandomDevice() { delete rand_stream_; };

    result_type operator()() { return rand_stream_->genrand_uint32(); }

private:
    VSFMT19937<SIMD_N_BITS, QM_Scalar>* rand_stream_;
};

class VMT19937BulkRandomDevice : public VMT19937RandomDeviceBase {
public:
    DELETE_COPY_MOVE(VMT19937BulkRandomDevice)

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return UINT_MAX; }

    VMT19937BulkRandomDevice()
    {
        MT19937Matrix const* jumpMatrix = new MT19937Matrix(
            "/home/yuzj/Documents/pbsim3_modern/deps/VMT19937-master/dat/mt/dat2/" + get_filename());
        rand_stream_ = new VMT19937<SIMD_N_BITS, QM_Any>(SEED_INIT, SEED_LENGTH, 0, nullptr, jumpMatrix);
        delete jumpMatrix;
    }

    ~VMT19937BulkRandomDevice() { delete rand_stream_; };

    void gen(std::vector<result_type>& vec)
    {
        result_type* p = vec.data();
        AlignedVector<uint32_t, 64> buffer(vec.size());
        rand_stream_->genrand_uint32_anySize(buffer.data(), vec.size());
        std::memcpy(p, buffer.data(), vec.size() * sizeof(result_type));
    }

private:
    VMT19937<SIMD_N_BITS, QM_Any>* rand_stream_;
};

class VSFMT19937BulkRandomDevice : public VMT19937RandomDeviceBase {
public:
    DELETE_COPY_MOVE(VSFMT19937BulkRandomDevice)

    static constexpr result_type min() { return 0; }

    static constexpr result_type max() { return UINT_MAX; }

    VSFMT19937BulkRandomDevice()
    {
        SFMT19937Matrix const* jumpMatrix
            = new SFMT19937Matrix("/home/yuzj/Documents/pbsim3_modern/deps/VMT19937-master/dat/sfmt/" + get_filename());
        rand_stream_ = new VSFMT19937<SIMD_N_BITS, QM_Any>(SEED_INIT, SEED_LENGTH, 0, nullptr, jumpMatrix);
        delete jumpMatrix;
    }

    ~VSFMT19937BulkRandomDevice() { delete rand_stream_; };

    void gen(std::vector<result_type>& vec)
    {
        result_type* p = vec.data();
        AlignedVector<uint32_t, 64> buffer(vec.size());
        rand_stream_->genrand_uint32_anySize(buffer.data(), vec.size());
        std::memcpy(p, buffer.data(), vec.size() * sizeof(result_type));
    }

private:
    VSFMT19937<SIMD_N_BITS, QM_Any>* rand_stream_;
};
