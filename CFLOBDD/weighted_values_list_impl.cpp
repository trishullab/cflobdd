#include "weighted_values_list.cpp"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
// typedef double BIG_FLOAT;
template class WeightedValuesListHandle<BIG_FLOAT>;
template std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<BIG_FLOAT> &r);


typedef boost::multiprecision::cpp_complex_double BIG_COMPLEX_FLOAT;
// typedef double BIG_COMPLEX_FLOAT;
template class WeightedValuesListHandle<BIG_COMPLEX_FLOAT>;
template std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<BIG_COMPLEX_FLOAT> &r);