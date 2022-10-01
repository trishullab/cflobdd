#include "cflobdd_top_node_t.cpp"
#include "matmult_map.h"
#include "matrix1234_float_boost.h"

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

    template class CFLOBDDTopNodeT<double>;
    template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr MkPlusTopNode<double>(
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr g);
    template CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr MkTimesTopNode<double>(
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr f,
        CFLOBDDTopNodeT<double>::CFLOBDDTopNodeTRefPtr g);
}