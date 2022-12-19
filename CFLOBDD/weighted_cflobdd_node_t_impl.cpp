
#include "weighted_cflobdd_node_t.cpp"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_complex_fb_mul.h"

namespace CFL_OBDD {

    template class WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDNodeHandleT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> &d);
    template class WeightedCFLOBDDForkNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template class WeightedCFLOBDDDontCareNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template class WeightedCFLOBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;

    template class WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> &d);
    template class WeightedCFLOBDDForkNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template class WeightedCFLOBDDDontCareNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template class WeightedCFLOBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
}