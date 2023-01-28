#include "weighted_cflobdd_top_node_t.cpp"
#include "wmatrix1234_fb_mul.h"
#include "wmatrix1234_complex_fb_mul.h"
#include "wmatrix1234_fourier_mul.h"

namespace CFL_OBDD {

    template class WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> &d);
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkLeftScalarTimesTopNode<BIG_FLOAT, BIG_FLOAT, std::multiplies<BIG_FLOAT>>(BIG_FLOAT c, WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        operator*<BIG_FLOAT, BIG_FLOAT, std::multiplies<BIG_FLOAT>>(BIG_FLOAT c, WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkTimesTopNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
    operator*<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f, WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    // template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
    //     ApplyAndReduce<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n1,
    //         WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n2,
    //         TimesFunc<BIG_FLOAT>
	//     );
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkPlusTopNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        operator+<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f, 
	        WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);


    template class WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> &d);
    template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkLeftScalarTimesTopNode<BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(BIG_COMPLEX_FLOAT c, WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        operator*<BIG_COMPLEX_FLOAT, BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(BIG_COMPLEX_FLOAT c, WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkTimesTopNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
    operator*<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f, WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    // template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
    //     ApplyAndReduce<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n1,
    //         WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n2,
    //         TimesFunc<BIG_COMPLEX_FLOAT>
	//     );
    template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkPlusTopNode<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        operator+<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f, 
	        WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);


    template class WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>;
    template std::ostream& operator<< (std::ostream & out, const WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>> &d);
    template WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr
        MkLeftScalarTimesTopNode<fourierSemiring, fourierSemiring, std::multiplies<fourierSemiring>>(fourierSemiring c, WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr
        operator*<fourierSemiring, fourierSemiring, std::multiplies<fourierSemiring>>(fourierSemiring c, WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr
        MkTimesTopNode<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr
    operator*<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr f, WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr g);
    // template WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
    //     ApplyAndReduce<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n1,
    //         WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n2,
    //         TimesFunc<BIG_COMPLEX_FLOAT>
	//     );
    template WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr
        MkPlusTopNode<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr
        operator+<fourierSemiring, std::multiplies<fourierSemiring>>(WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr f, 
	        WeightedCFLOBDDTopNodeT<fourierSemiring, std::multiplies<fourierSemiring>>::WeightedCFLOBDDTopNodeTRefPtr g);
}