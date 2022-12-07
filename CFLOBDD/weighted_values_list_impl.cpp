#include "weighted_values_list.cpp"
#include <boost/multiprecision/cpp_dec_float.hpp>

typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
template class WeightedValuesListHandle<BIG_FLOAT>;
template std::ostream& operator<< (std::ostream & out, const WeightedValuesListHandle<BIG_FLOAT> &r);