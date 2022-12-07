#include "weighted_cflobdd_top_node_t.cpp"
#include "wmatrix1234_fb_mul.h"

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
    //     ApplyAndReduce(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n1,
    //         WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr n2,
    //         T(*func)(T, T)
	//     );
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        MkPlusTopNode<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f,
                WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);
    template WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr
        operator+<BIG_FLOAT, std::multiplies<BIG_FLOAT>>(WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr f, 
	        WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr g);

}