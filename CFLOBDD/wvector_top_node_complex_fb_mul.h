#ifndef W_VECTOR_FLOAT_BOOST_TOP_NODE_GUARD
#define W_VECTOR_FLOAT_BOOST_TOP_NODE_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include <boost/multiprecision/cpp_complex.hpp>
#include "wvector_complex_fb_mul.h"
#include "return_map_T.h"
#include "weighted_cflobdd_top_node_t.h"
// namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_complex_100 BIG_COMPLEX_FLOAT;
	// typedef double BIG_COMPLEX_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_COMPLEX_FLOAT;
	typedef ReturnMapBody<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_COMPLEX_FLOAT> ComplexFloatBoostReturnMapHandle;
}
namespace CFL_OBDD {

    typedef WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>> WeightedCFLOBDDTopNodeComplexFloatBoost;
	typedef WeightedCFLOBDDTopNodeT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr;

	namespace WeightedVectorComplexFloatBoostMul {

		// Initialization routine
		extern void VectorInitializerTop();

		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkBasisVectorTop(unsigned int level, unsigned int index, int cflobdd_kind = 1); // Representation of basis vector at index
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkBasisVectorTop(unsigned int level, std::string s, int cflobdd_kind = 1); // Representation of basis vector at index
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr VectorToMatrixInterleavedTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n); // Convert vector to matrix with variables interleaved
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr MkColumn1MatrixTop(unsigned int level);
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr NoDistinctionNodeTop(unsigned int level, BIG_COMPLEX_FLOAT val);
		extern std::string SamplingTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n, std::mt19937 mt, std::uniform_real_distribution<double> dis, bool voctwo = false, std::string = "");
		extern WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr VectorWithAmplitudeTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n);
		extern void VectorPrintColumnMajorTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);
		extern void VectorPrintColumnMajorInterleavedTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n, std::ostream & out);
		extern long double getNonZeroProbabilityTop(WeightedCFLOBDDTopNodeComplexFloatBoostRefPtr n);
	}
}

#endif

