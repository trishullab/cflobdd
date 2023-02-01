#ifndef W_MATRIX1234_COMPLEX_FB_MUL_TOP_NODE_GUARD
#define W_MATRIX1234_COMPLEX_FB_MUL_TOP_NODE_GUARD

#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_complex.hpp>
#include "wmatrix1234_complex_fb_mul.h"
#include "return_map_T.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
    // typedef double BIG_COMPLEX_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_COMPLEX_FLOAT;
	typedef ReturnMapBody<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapHandle;
	typedef ReturnMapBody<int> CFLOBDDReturnMapBody;
	typedef ReturnMapHandle<int> CFLOBDDReturnMapHandle;
}

namespace CFL_OBDD {

	typedef WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDTopNodeComplexFloatBoost;
	typedef WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr;
	

	namespace WeightedMatrix1234ComplexFloatBoostMul {

		// Initialization routine
		extern void Matrix1234InitializerTop();

		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkIdRelationInterleavedTop(unsigned int i, int cflobdd_kind = 1); // Representation of identity relation
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkWalshInterleavedTop(unsigned int i, int cflobdd_kind = 1); // Representation of Walsh matrix
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkNegationMatrixInterleavedTop(unsigned int i, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkPauliYGateTop(unsigned int i, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkPauliZGateTop(unsigned int i, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkSGateTop(unsigned int i, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkCNOTInterleavedTop(unsigned int i);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkExchangeInterleavedTop(unsigned int i); // Representation of exchange matrix
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkCNOTTopNode(unsigned int level, unsigned int n, long int controller, long int controlled, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr KroneckerProduct2VocsTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr m1, WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr m2);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkCCNOTTop(unsigned int level, long int controller1, long int controller2, long int controlled, int cflobdd_kind = 1);

		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkRestrictTop(unsigned int level, std::string s, int cflobdd_kind = 1);

		extern void MatrixPrintRowMajorTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);
		extern void MatrixPrintRowMajorInterleavedTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);

		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MatrixMultiplyV4TopNode(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr c1, WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr c2);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkSwapGateTop(unsigned int level, long int i, long int j, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkiSwapGateTop(unsigned int level, long int i, long int j, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkCPGateTop(unsigned int level, long int i, long int j, double theta, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkCZGateTop(unsigned int level, long int i, long int j, double theta, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkPhaseShiftGateTop(unsigned int level, double theta, int cflobdd_kind = 1);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkCSwapGateTop(unsigned int level, long int c, long int i, long int j, int cflobdd_kind = 1);
	}
}

#endif

