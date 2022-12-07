#include "weighted_cflobdd_top_node_t.cpp"
#include "wmatrix1234_fb_mul.h"

namespace CFL_OBDD {

    typedef WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDTopNode; 
    template class WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNode &d);
}