#pragma once

#include <cstdint>

// nAlign must be a multiple of 2 and no more than 64
template <typename T, size_t nAlign>
inline T* myAlignedNew(size_t n)
{
    // allocate memory
    char *p = new char[n * sizeof(T) + nAlign];

    // align pointer
    auto addr = reinterpret_cast<std::uintptr_t>(p);
    size_t moveFwdBy = (nAlign - (addr % nAlign));
    p += moveFwdBy;
    p[-1] = (char) moveFwdBy;

    return (T*)p;
}

inline void myAlignedDelete(void *p)
{
    if (p) {
        char* p8 = (char*)p;
        p8 -= (size_t)p8[-1];
        delete [] p8;
    }
}

template <typename T, unsigned nAlign>
class AlignedVector
{
    T* m_data;
    size_t m_n;

    void deallocate()
    {
        if (m_data) {
            myAlignedDelete(m_data);
            m_data = nullptr;
        }
    }

public:
    AlignedVector() : m_data(nullptr), m_n(0) {}
    AlignedVector(size_t n) : m_data(nullptr), m_n(0) { init(n); }
    ~AlignedVector() { deallocate(); }

    void init(size_t n)
    {
        if (n != m_n) {
            deallocate();
            m_data = myAlignedNew<T, nAlign>(n);
            m_n = n;
        }
    }

    T& operator[](size_t i) { return m_data[i]; }
    const T& operator[](size_t i) const { return m_data[i]; }
    T* data() { return m_data; }
    const T* data() const { return m_data; }
};

// Given a cube with dimension [outerDim][MidDim][InnerDim]
// Extract values with midDim=midIndex and store in a matrix of dimensions [outerDim][InnerDim]
template <size_t MidDim, size_t InnerDim, typename T>
inline void cubeToMatrix(T* vec, const T* mat, size_t outerDim, size_t midIndex)
{
    mat += InnerDim * midIndex;
    for (size_t w = 0; w < outerDim; ++w)
        for (size_t i = 0; i < InnerDim; ++i)
            vec[w * InnerDim + i] = mat[w * (MidDim * InnerDim) + i];
}

// Given a cube with dimension [outerDim][MidDim][InnerDim]
// Copy a matrix of dimensions [outerDim][InnerDim] to the layer in the cube where midDim=midIndex
template <size_t MidDim, size_t InnerDim, typename T>
inline void matrixToCube(T* mat, const T* vec, size_t outerDim, size_t midIndex)
{
    mat += InnerDim * midIndex;
    for (size_t w = 0; w < outerDim; ++w)
        for (size_t i = 0; i < InnerDim; ++i)
            mat[w * (MidDim * InnerDim) + i] = vec[w * InnerDim + i];
}

