#include "weighted_connectionT.cpp"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"

namespace CFL_OBDD {

    template class WConnection<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WConnection<BIG_FLOAT, std::multiplies<BIG_FLOAT>> &c);

    template class WConnection<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WConnection<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> &c);

    template class WConnection<fourierSemiring, std::multiplies<fourierSemiring>>;
    template std::ostream& operator<< (std::ostream & out, const WConnection<fourierSemiring, std::multiplies<fourierSemiring>> &c);
    
}