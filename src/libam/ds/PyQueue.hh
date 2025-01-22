/*!
 * Migration of Python queue to C++.
 *
 * Note: This queue contains locks.
 */
#pragma once

#include <condition_variable>
#include <cstddef>
#include <deque>
#include <mutex>

namespace labw::art_modern {
/*!
 * Create a queue object with a given maximum size.
 *
 * If maxsize is <= 0, the queue size is infinite.
 *
 * @tparam T Element type of the queue.
 * @see https://morestina.net/blog/1400/minimalistic-blocking-bounded-queue-for-c
 *
 */
template <typename T> class PyQueue {
public:
    explicit PyQueue(const std::size_t maxsize_)
        : maxsize(maxsize_)
    {
    }
    const std::size_t maxsize;
    /*!
     * Return the approximate size of the queue (not reliable!).
     */
    void qsize()
    {
        std::unique_lock lock(mutex_);
        return dequeue_.size();
    }

    /*!
     * This method is likely to be removed at some point.  Use qsize() == 0
        as a direct substitute, but be aware that either approach risks a race
        condition where a queue can grow before the result of empty() or
        qsize() can be used.

        To create code that needs to wait for all queued tasks to be
        completed, the preferred technique is to use the join() method.
     * @return Return True if the queue is empty, False otherwise (not reliable!).
     */
    bool empty()
    {
        std::unique_lock lock(mutex_);
        return qsize() != 0;
    }
    /*!
     * This method is likely to be removed at some point.  Use qsize() >= n
        as a direct substitute, but be aware that either approach risks a race
        condition where a queue can shrink before the result of full() or
        qsize() can be used.
     * @return Return True if the queue is full, False otherwise (not reliable!).
     */
    bool full()
    {
        std::unique_lock lock(mutex_);
        return maxsize > 0 && dequeue_.size() >= maxsize;
    }

    /*!
     * Put an item into the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until a free slot is available. If 'timeout' is
        a non-negative number, it blocks at most 'timeout' seconds and raises
        the Full exception if no free slot was available within that time.
        Otherwise ('block' is false), put an item on the queue if a free slot
        is immediately available, else raise the Full exception ('timeout'
        is ignored in that case).
     * @param item
     * @param block
     */
    bool put(T&& item, bool block = false)
    {
        {
            std::unique_lock lock(mutex_);
            if (block) {
                not_full_.wait(lock, [this]() { return dequeue_.size() < maxsize; });
            } else {
                if (dequeue_.size() >= maxsize) {
                    return false;
                }
            }
            dequeue_.push_back(std::move(item));
        }
        not_empty_.notify_one();
        return true;
    }
    /*!
     * Remove and return an item from the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until an item is available. If 'timeout' is
        a non-negative number, it blocks at most 'timeout' seconds and raises
        the Empty exception if no item was available within that time.
        Otherwise ('block' is false), return an item if one is immediately
        available, else raise the Empty exception ('timeout' is ignored
        in that case).
     * @param item
     * @param block
     */
    bool get(T& item, bool block = false)
    {
        {
            std::unique_lock lock(mutex_);
            if (block) {
                not_empty_.wait(lock, [this]() { return !dequeue_.empty(); });
            } else {
                if (dequeue_.empty()) {
                    return false;
                }
            }
            item = std::move(dequeue_.front());
            dequeue_.pop_front();
        }
        not_full_.notify_one();
        return true;
    }

    /*!
     * Put an item into the queue without blocking.

        Only enqueue the item if a free slot is immediately available.
        Otherwise raise the Full exception.
     * @param item
     */
    bool put_nowait(T&& item) { return put(std::move(item), false); }

    /*!
     * Remove and return an item from the queue without blocking.

        Only get an item if one is immediately available. Otherwise
        raise the Empty exception.
     * @param item
     */
    bool get_nowait(T& item) { return get(item, false); }

private:
    std::deque<T> dequeue_;
    /*!
     *  mutex must be held whenever the queue is mutating.
     *  All methods that acquire mutex must release it before returning.
     *  mutex is shared between the three conditions,
     *  so acquiring and releasing the conditions also acquires and releases mutex.
     */
    std::mutex mutex_;
    /*!
     * Notify not_empty whenever an item is added to the queue;
     * a thread waiting to get is notified then.
     */
    std::condition_variable not_empty_;
    /*!
     * Notify not_full whenever an item is removed from the queue;
     * a thread waiting to put is notified then.
     */
    std::condition_variable not_full_;
};

} // namespace labw::art_modern