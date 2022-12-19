#ifndef W_VECTOR_FLOAT_BOOST_TOP_NODE_GUARD
#define W_VECTOR_FLOAT_BOOST_TOP_NODE_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "wvector_fb_mul.h"
#include "return_map_T.h"
#include "weighted_cflobdd_top_node_t.h"
// namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	// typedef double BIG_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef ReturnMapBody<BIG_FLOAT> FloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_FLOAT> FloatBoostReturnMapHandle;
}
namespace CFL_OBDD {

    typedef WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>> WeightedCFLOBDDTopNodeFloatBoost;
	typedef WeightedCFLOBDDTopNodeT<BIG_FLOAT, std::multiplies<BIG_FLOAT>>::WeightedCFLOBDDTopNodeTRefPtr WeightedCFLOBDDTopNodeFloatBoostRefPtr;

	namespace WeightedVectorFloatBoostMul {

		// Initialization routine
		extern void VectorInitializerTop();

		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkBasisVectorTop(unsigned int level, unsigned int index); // Representation of basis vector at index
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkBasisVectorTop(unsigned int level, std::string s); // Representation of basis vector at index
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr VectorToMatrixInterleavedTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n); // Convert vector to matrix with variables interleaved
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr MkColumn1MatrixTop(unsigned int level);
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr NoDistinctionNodeTop(unsigned int level, BIG_FLOAT val);
		extern std::string SamplingTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, bool voctwo = false, std::string = "");
		extern WeightedCFLOBDDTopNodeFloatBoostRefPtr VectorWithAmplitudeTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n);
		extern void VectorPrintColumnMajorTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);
		extern void VectorPrintColumnMajorInterleavedTop(WeightedCFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);
	}
}

#endif

