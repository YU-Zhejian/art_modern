/*
 * SharedQueue.h
 *
 * A shared producer-consumer queue implementation with blocking put()s and
 * get()s, and size determined at construction.
 */

#ifndef _READ_QUEUE_H
#define _READ_QUEUE_H

#include "Lock.h"
#include <stddef.h>

template <typename T>
class SharedQueue {
private:
	Semaphore  m_filled_slots;
	Semaphore  m_empty_slots;
	Lock       m_lock;
	unsigned   m_front;
	unsigned   m_back;
	T	  *m_array;
	size_t	   m_size;
public:

	SharedQueue(size_t size)
		: m_filled_slots(0),
		  m_empty_slots(size),
		  m_lock(),
		  m_front(0),
		  m_back(size - 1),
		  m_array(new T[size]),
		  m_size(size)
	{ }

	~SharedQueue() {
		int n = (int)m_filled_slots;
		for (int i = 0; i < n; i++)
			delete m_array[(m_front + i) % m_size];
		delete [] m_array;
	}

	/* Retrieve a T from the queue as soon as there is one available. */
	T get() {
		m_filled_slots.wait();
		m_lock.lock();

		T t = m_array[m_front];
		m_front = (m_front + 1) % m_size;

		m_empty_slots.post();
		m_lock.unlock();
		return t;
	}
	/* Place the T in the queue as soon as there is an empty space. */
	void put(T t) {
		m_empty_slots.wait();
		m_lock.lock();

		m_back = (m_back + 1) % m_size;
		m_array[m_back] = t;

		m_filled_slots.post();
		m_lock.unlock();
	}
};

#endif
