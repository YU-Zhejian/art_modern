/**
 * Modified from <https://github.com/wangeddie67/ringbuffer_opt_demo>
 *
 * TODO: This file have no open-source license.
 */
#include <atomic>
#include <cstdlib>
#include <cstring>

// Entry in ring buffer
template <typename T> struct BufferEntry {
    std::atomic<std::size_t> m_sn { 0 }; // Sequential number.
    std::atomic<T*> mp_ptr { nullptr }; // Pointer to data entry.
    BufferEntry() = default;
} __attribute__((aligned(16)));

// Data structure for ring buffer
template <typename T> class RingBuffer {
public:
    alignas(64) std::atomic<std::size_t> m_size;
    alignas(64) std::atomic<std::size_t> m_head { 0 };
    alignas(64) std::atomic<std::size_t> m_tail { 0 };
    alignas(64) BufferEntry<T>* mp_entries;

    explicit RingBuffer(const std::size_t entry_num)
        : m_size(entry_num)
        , mp_entries(static_cast<BufferEntry<T>*>(std::calloc(entry_num << 2, sizeof(BufferEntry<T>))))
    {
        // No need to memset due to calloc.
        // std::memset(mp_entries, 0, sizeof(BufferEntry) * m_size * 4);

        for (std::size_t i = 0; i < m_size; i++) {
            mp_entries[i << 2].m_sn = i;
        }
    }
    int enqueue_ringbuf(T* entry)
    {
        std::size_t const sn = std::atomic_fetch_add(&m_tail, 1);
        std::size_t const enq_ptr = (sn % m_size) << 2;

        std::size_t entry_sn = 0;
        T* entry_ptr = nullptr;
        do {
            entry_sn = mp_entries[enq_ptr].m_sn.load(std::memory_order_relaxed);
            entry_ptr = mp_entries[enq_ptr].mp_ptr.load(std::memory_order_relaxed);
        } while (entry_sn != sn || entry_ptr != nullptr);

        mp_entries[enq_ptr].mp_ptr.store(entry, std::memory_order_release);

        return 0;
    }

    int dequeue_ringbuf(T** entry)
    {
        std::size_t const sn = std::atomic_fetch_add(&m_head, 1);
        std::size_t const deq_ptr = (sn % m_size) << 2;

        std::size_t entry_sn = 0;
        T* entry_ptr = nullptr;
        do {
            entry_sn = mp_entries[deq_ptr].m_sn.load(std::memory_order_relaxed);
            entry_ptr = mp_entries[deq_ptr].mp_ptr.load(std::memory_order_relaxed);
        } while (entry_ptr == nullptr || entry_sn != sn);

        *entry = entry_ptr;
        mp_entries[deq_ptr].mp_ptr.store(nullptr, std::memory_order_relaxed);
        mp_entries[deq_ptr].m_sn.store(sn + m_size, std::memory_order_release);
        return 0;
    }
};
