
#include "weighted_bdd_node_t.cpp"
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"

namespace CFL_OBDD {

    template class WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> &d);
    // template class WeightedBDDNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template class WeightedBDDInternalNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template class WeightedBDDLeafNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;

    template class WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedBDDNodeHandle<BIG_FLOAT, std::multiplies<BIG_FLOAT>> &d);
    // template class WeightedBDDNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template class WeightedBDDInternalNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template class WeightedBDDLeafNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;

    template class WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedBDDNodeHandle<fourierSemiring, std::multiplies<fourierSemiring>> &d);
    // template class WeightedBDDNode<fourierSemiring, std::multiplies<fourierSemiring>>;
    template class WeightedBDDInternalNode<fourierSemiring, std::multiplies<fourierSemiring>>;
    template class WeightedBDDLeafNode<fourierSemiring, std::multiplies<fourierSemiring>>;

}