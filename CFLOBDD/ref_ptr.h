#ifndef REF_PTR_H_
#define REF_PTR_H_

#include <cassert>
#include <stdio.h>
#include <limits.h>
#include <iostream>

/* A reference counting pointer class
 *
 * This class is *NOT* thread safe.
 *
 * The templated class must have a member variable named count
 * that can be accessed from ref_ptr and modified by ref_ptr.
 * the count variable should have operator++(),operator--(), and
 * operator==( int ) defined.  As a note, this class was designed
 * with count being an unsigned integer.
 *
 * Count should be initialized to 0 for proper reference
 * couting to work.  If it is desirable for the pointer/object
 * to never be deleted, initialize count >= 1.  WARNING: don't forget
 * to initialize count to 0 in the copy constructor!  (And don't
 * forget to have a copy constructor.)
 *
 * ref_ptr provides a templated function operator<< that allows 
 * the use of std::ostream << with ref_ptr<T> objects.  To take
 * advantage of this, the templated class must provide a method
 * with signature:
 *
 *         std::ostream& print( std::ostream& )
 *
 */

class RefCounter {
    friend std::ostream& operator<<( std::ostream& o, const RefCounter& rc )
    {
        return (o << rc.cnt);
    }
 public:
    RefCounter(unsigned _cnt = 0) : cnt(_cnt) { }
    RefCounter(const RefCounter &that) : cnt(0) { }
    RefCounter& operator++() { ++cnt; return *this; }
    RefCounter& operator--() { --cnt; return *this; }
    RefCounter & operator=(const RefCounter &that) { return *this; }
    bool operator==(unsigned i) const { return (unsigned)i == cnt; }
 private:
    unsigned cnt;
};

inline bool operator==(unsigned i, const RefCounter &rc) { return rc == i; }

template< typename T > class ref_ptr {

    public:
        ref_ptr( T *t = 0 ) {
            acquire( t );
        }

        ref_ptr( const ref_ptr& rp ) {
            acquire( rp.ptr );
        }

        ~ref_ptr() {
            release();
        }

        ref_ptr& operator=( T* t ) {
          if( ptr != t ) {
              T* old_ptr = ptr;
              acquire(t);
              release(old_ptr);
          }
          return *this;
        }

        inline ref_ptr& operator=( const ref_ptr& rp ) {
            return *this = rp.ptr;
        }

        inline bool operator==(const ref_ptr& that) const {
            return ptr == that.ptr;
        }

        T * get_ptr() const {
            return ptr;
        }

        T * operator->() const {
            return ptr;
        }

        T& operator*()   const {
          /*            CND_DBGS( (0 == ptr), 
                        fprintf(stderr, "ref_ptr::operator*: Invalid ptr.\n") ); */
            /* for NDEBUG. Make sure it crashes
               here rather than have some
               undefined behaviour. */
            assert(0 != ptr);
            return *ptr;
        }

    private:
        T *ptr;

        void acquire( T *t ) {
            ptr = t;
            if( t ) {
                ++t->count;
#ifdef DBGREFPTR
                std::cout << "Acquired " << t << " with count = "
                    << t->count << std::endl;
#endif
            }
        }

        static void release(T* old_ptr)
        {
            if( old_ptr ) {
                --old_ptr->count;
#ifdef DBGREFPTR
                std::cout << "Released " << *old_ptr << " with count = "
                    << old_ptr->count << std::endl;
#endif
                if( old_ptr->count == 0 ) {
#ifdef DBGREFPTR
                    /*std::cout << "Deleting ptr: " << *old_ptr << std::endl;*/
					std::cout << "Deleting ptr: " << old_ptr->count << std::endl;
#endif
                    delete old_ptr;
                }
            }
        }

        inline void release() {
          release(ptr);
          ptr = 0;
        }
};

/*
template< typename T > std::ostream& operator<<(
        std::ostream& o,
        const ref_ptr< T >& r )
{
    if( 0 != r.get_ptr() )
        return r->print(o);
    else
        return (o << "NULL");
}
*/

template< typename T > std::ostream& operator<<(
        std::ostream& o,
        const ref_ptr< T > r )
{
    if( 0 != r.get_ptr() )
        return r->print(o);
    else
        return (o << "NULL");
}

#endif    // REF_PTR_H_
/* Yo, Emacs!
;;; Local Variables: ***
;;; tab-width: 4 ***
;;; End: ***
*/
