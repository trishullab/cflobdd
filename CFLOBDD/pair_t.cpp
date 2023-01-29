#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "pair_t.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include "fourier_semiring.h"

template <typename T, typename T1>
Pair_T<T,T1>::Pair_T()
{
}

template <typename T, typename T1>
Pair_T<T,T1>::Pair_T(const T i1, const T1 i2)
  :  first(i1), second(i2)
{
}

template <typename T, typename T1>
std::ostream& operator<< (std::ostream & out, const Pair_T<T,T1> &p)
{
  out << "(" << p.First() << ", " << p.Second() << ")";
  return(out);
}

template <typename T, typename T1>
bool operator==(const Pair_T<T,T1>& lhs, const Pair_T<T,T1>& rhs)
{
	return (lhs.First() == rhs.First()) && (lhs.Second() == rhs.Second());
}

template <typename T, typename T1>
Pair_T<T,T1>& Pair_T<T,T1>::operator= (const Pair_T<T,T1>& i)
{
  if (this != &i)      // don't assign to self!
  {
    first = i.first;
    second = i.second;
  }
  return *this;        
}

template <typename T, typename T1>
Pair_T<T,T1> Pair_T<T,T1>::operator! ()
{
	int newFirst = !this->First();
	int newSecond = !this->Second();
	return Pair_T(newFirst,newSecond);
}

template <>
Pair_T<fourierSemiring, fourierSemiring> Pair_T<fourierSemiring,fourierSemiring>::operator! ()
{
	abort();
}

// Overloaded !=
template <typename T, typename T1>
bool Pair_T<T,T1>::operator!=(const Pair_T<T,T1>& p) const
{
  return (first != p.first) || (second != p.second);
}

// Overloaded ==
template <typename T, typename T1>
bool Pair_T<T,T1>::operator==(const Pair_T<T,T1>& p)
{
  return (first == p.first) && (second == p.second);
}


typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
template class Pair_T<BIG_FLOAT, BIG_FLOAT>;

typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
template class Pair_T<BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT>;
template class Pair_T<fourierSemiring, fourierSemiring>;