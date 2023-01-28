#ifndef PAIR_T_GUARD
#define PAIR_T_GUARD

#include <iostream>
#include <fstream>
#include <boost/functional/hash.hpp>
#define MODSIZE 2019

template <typename T, typename T1>
class Pair_T {
 public:
  Pair_T();                              // Default constructor
  Pair_T(const T i1, const T1 i2);    // Constructor
  Pair_T& operator= (const Pair_T& p);  // Overloaded assignment
  Pair_T operator! ();
  bool operator!= (const Pair_T& p) const;     // Overloaded !=
  bool operator== (const Pair_T& p);     // Overloaded ==
  T First() const { return first; }     // Access function
  T1 Second() const { return second; }   // Access function
  struct Pair_T_hash {
	  size_t operator()(const Pair_T& p) const {
          boost::hash<T> hash1;
          boost::hash<T1> hash2;
		  return (117 * (hash1(p.First()) % MODSIZE + 1) + hash2(p.Second()) % MODSIZE) % MODSIZE;
	  }
  };
 private:
  T first;
  T1 second;
};

template <typename T, typename T1>
bool operator==(const Pair_T<T,T1>& lhs, const Pair_T<T,T1>& rhs);
template <typename T, typename T1>
std::ostream& operator<< (std::ostream & out, const Pair_T<T,T1> &p);

#endif
