#include "weighted_cflobdd_top_node_t.cpp"
#include "weighted_cflobdd_t.h"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"

namespace CFL_OBDD {

    typedef WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDTopNode; 
    template class WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNode &d);

    typedef WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDTopNodeComplex; 
    template class WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeComplex &d);

    // template class WEIGHTED_CFLOBDD_T<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    typedef WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>> WeightedCFLOBDDTopNodeFourier; 
    template class WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeFourier &d);
}