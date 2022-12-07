#include "weighted_rootConnectionT.cpp"
#include "wmatrix1234_fb_mul.h"

namespace CFL_OBDD {

    template class WRootConnection<ReturnMapHandle<BIG_FLOAT>, BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WRootConnection<ReturnMapHandle<BIG_FLOAT>, BIG_FLOAT, std::multiplies<BIG_FLOAT>> &c);
    
}