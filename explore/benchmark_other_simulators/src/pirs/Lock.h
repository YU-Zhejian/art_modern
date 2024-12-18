/*
 * Lock.h
 *
 * Lock and Semaphore classes to wrap around pthread_mutex_t and sem_t.
 */

#ifndef _LOCK_H
#define _LOCK_H

#include <pthread.h>
#include <semaphore.h>

class Lock {
private:
	pthread_mutex_t _lock;
public:
	Lock()	     { pthread_mutex_init(&_lock, NULL);}
	~Lock()	     { pthread_mutex_destroy(&_lock);	}
	void set()   { pthread_mutex_lock(&_lock);	}
	void lock()  { set();				}
	void unset() { pthread_mutex_unlock(&_lock);	}
	void unlock(){ unset();				}
	bool test()  { return pthread_mutex_trylock(&_lock) == 0; }
};

class Semaphore {
private:
	sem_t sem;
public:
	Semaphore() {
		sem_init(&sem, 0, 1);
	}

	Semaphore(unsigned val) {
		sem_init(&sem, 0, val);
	}

	~Semaphore() {
		sem_destroy(&sem);
	}

	void post() {
		sem_post(&sem);
	}

	void wait() {
		sem_wait(&sem);
	}

	int val() {
		int n;
		sem_getvalue(&sem, &n);
		return n;
	}
	operator int() {
		return val();
	}
};

#endif
