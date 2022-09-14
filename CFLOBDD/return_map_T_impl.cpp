#include "return_map_T.h"

template class ReturnMapHandle<int>;
template std::ostream& operator<< (std::ostream & out, const ReturnMapHandle<int> &r);