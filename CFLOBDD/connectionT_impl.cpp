#include "connectionT.cpp"
#include "matrix1234_float_boost.h"
#include "fourier_semiring.h"
#include "matrix1234_complex_float_boost.h"

namespace CFL_OBDD {
    template class  ConnectionT<ReturnMapHandle<int>>;
	template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<int>> &c);

    template class  ConnectionT<ReturnMapHandle<MatMultMapHandle>>;
	template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<MatMultMapHandle>> &c);

    template class  ConnectionT<ReturnMapHandle<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<BIG_FLOAT>> &c);

    template class  ConnectionT<ReturnMapHandle<double>>;
    template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<double>> &c);

    template class  ConnectionT<ReturnMapHandle<fourierSemiring>>;
    template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<fourierSemiring>> &c);

    // template class  ConnectionT<ReturnMapHandle<std::complex<double>>>;
    // template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<std::complex<double>>> &c);

    template class  ConnectionT<ReturnMapHandle<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<BIG_COMPLEX_FLOAT>> &c);
}