#include "weighted_rootConnectionT.cpp"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"

namespace CFL_OBDD {

    template class WRootConnection<ReturnMapHandle<BIG_FLOAT>, BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WRootConnection<ReturnMapHandle<BIG_FLOAT>, BIG_FLOAT, std::multiplies<BIG_FLOAT>> &c);

    template class WRootConnection<ReturnMapHandle<BIG_COMPLEX_FLOAT>, BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WRootConnection<ReturnMapHandle<BIG_COMPLEX_FLOAT>, BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> &c);

    template class WRootConnection<ReturnMapHandle<fourierSemiring>, fourierSemiring, std::multiplies<fourierSemiring>>;
    template std::ostream& operator<< (std::ostream & out, const WRootConnection<ReturnMapHandle<fourierSemiring>, fourierSemiring, std::multiplies<fourierSemiring>> &c);
    
}