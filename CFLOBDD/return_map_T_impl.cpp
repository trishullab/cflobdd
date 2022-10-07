#include "return_map_T.h"

template class ReturnMapHandle<int>;
template std::ostream& operator<< (std::ostream & out, const ReturnMapHandle<int> &r);

// template class ReturnMapHandle<CFL_OBDD::BIG_COMPLEX_FLOAT>;
// template std::ostream& operator<< (std::ostream & out, const ReturnMapHandle<CFL_OBDD::BIG_COMPLEX_FLOAT> &r);