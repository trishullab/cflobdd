#include "connectionT.cpp"
#include "matrix1234_float_boost.h"

namespace CFL_OBDD {
    template class  ConnectionT<ReturnMapHandle<int>>;
	template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<int>> &c);

    template class  ConnectionT<ReturnMapHandle<MatMultMapHandle>>;
	template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<MatMultMapHandle>> &c);

    template class  ConnectionT<ReturnMapHandle<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<BIG_FLOAT>> &c);

    template class  ConnectionT<ReturnMapHandle<double>>;
    template std::ostream& operator<< (std::ostream & out, const ConnectionT<ReturnMapHandle<double>> &c);
}