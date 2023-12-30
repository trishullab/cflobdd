#include "cflobdd_top_node_t.cpp"
#include "matmult_map.h"
#include "matrix1234_float_boost.h"
#include "fourier_semiring.h"
#include "matrix1234_complex_float_boost.h"

namespace CFL_OBDD {
    template class CFLOBDDTopNodeT<int>;

    template std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNodeT<int> &d);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr ApplyAndReduce<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr n1,
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr n2,
        BoolOp op);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr ApplyAndReduce<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr n1,
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr n2,
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr n3,
        BoolOp3 op);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr MkPlusTopNode<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr operator+<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f, 
	    CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr MkExorTopNode<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr operator^<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f, 
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<int, int>(
        int c, CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr operator*<int, int>(
        int c, CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr MkRightScalarTimesTopNode<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f, int c);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr operator*<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f, int c);
    template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<int>(
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);
    // template CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr operator*<int>(
    //     CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr f,
    //     CFLOBDDTopNodeT<int>::CFLOBDDTopNodeTRefPtr g);

    template class CFLOBDDTopNodeT<MatMultMapHandle>;
    template CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr operator*<MatMultMapHandle, VAL_TYPE>(
        VAL_TYPE c, CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<MatMultMapHandle, VAL_TYPE>(
        VAL_TYPE c, CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr MkPlusTopNode<MatMultMapHandle>(
        CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr operator+<MatMultMapHandle>(
        CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr f, 
	    CFLOBDDTopNodeT<MatMultMapHandle>::CFLOBDDTopNodeTRefPtr g);

    template class CFLOBDDTopNodeT<BIG_FLOAT>;
    template std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNodeT<BIG_FLOAT> &d);
    template CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr MkPlusTopNode<BIG_FLOAT>(
        CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<BIG_FLOAT>(
        CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<BIG_FLOAT, int>(
        int c, CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<BIG_FLOAT, BIG_FLOAT>(
        BIG_FLOAT c, CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<BIG_FLOAT, double>(
        double c, CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr g);

    template class CFLOBDDTopNodeT<double>;
    template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr MkPlusTopNode<double>(
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<double>(
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<double, int>(
        int c, CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<double, double>(
        double c, CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr g);    
    // template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr ApplyAndReduce<double>(
    //     CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr n1,
    //     CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr n2,
    //     double(*func)(double, double) op);

    template class CFLOBDDTopNodeT<fourierSemiring>;
    template CFLOBDDTopNodeT<fourierSemiring>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<fourierSemiring>(
        CFLOBDDTopNodeT<fourierSemiring>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<fourierSemiring>::CFLOBDDTopNodeTRefPtr g);

    // template class CFLOBDDTopNodeT<std::complex<double>>;
    // template CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<std::complex<double>>(
    //     CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr f,
    //     CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr g);
    // template CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<std::complex<double>, int>(
    //     int c, CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr g);
    // template CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<std::complex<double>, double>(
    //     double c, CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr g);
    // template CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr operator*<std::complex<double>, int>(
    //     int c, CFLOBDDTopNodeT<std::complex<double>>::CFLOBDDTopNodeTRefPtr f);
    
    template class CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>;
    template std::ostream& operator<< (std::ostream & out, const CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT> &d);
    template CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<BIG_COMPLEX_FLOAT>(
        CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<BIG_COMPLEX_FLOAT, int>(
        int c, CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<BIG_COMPLEX_FLOAT, double>(
        double c, CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr MkLeftScalarTimesTopNode<BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT>(
        BIG_COMPLEX_FLOAT c, CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr MkPlusTopNode<BIG_COMPLEX_FLOAT>(
        CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<BIG_COMPLEX_FLOAT>::CFLOBDDTopNodeTRefPtr g);
}