#ifndef W_MATRIX1234_FB_MUL_TOP_NODE_GUARD
#define W_MATRIX1234_FB_MUL_TOP_NODE_GUARD

#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "wmatrix1234_fb_mul.h"
#include "return_map_T.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
    // typedef double BIG_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef ReturnMapBody<BIG_FLOAT> FloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_FLOAT> FloatBoostReturnMapHandle;
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
}

namespace CFL_OBDD {

	typedef WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDTopNodeFloatBoost;
	typedef WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr WeightedCFLOBDDTopNodeFloatBoostRefPtr;
	

	namespace WeightedMatrix1234FloatBoostMul {

		// Initialization routine
		extern void Matrix1234InitializerTop();

		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i); // Representation of identity relation
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i); // Representation of Walsh matrix
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i);
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i);
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i); // Representation of exchange matrix
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled);
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr KroneckerProduct2VocsTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr m1, WeightedCFLOBDDTopNodeFloatBoostRefPtr m2);
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled);

		extern void MatrixPrintRowMajorTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorInterleavedTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);

		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MatrixMultiplyV4TopNode(WeightedCFLOBDDTopNodeFloatBoostRefPtr c1, WeightedCFLOBDDTopNodeFloatBoostRefPtr c2);
	}
}

#endif

