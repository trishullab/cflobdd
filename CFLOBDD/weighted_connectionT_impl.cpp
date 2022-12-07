#include "weighted_connectionT.cpp"
#include "wmatrix1234_fb_mul.h"

namespace CFL_OBDD {

    template class WConnection<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WConnection<BIG_FLOAT, std::multiplies<BIG_FLOAT>> &c);
    
}