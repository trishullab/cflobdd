#ifndef VECTOR_FLOAT_BOOST_TOP_NODE_GUARD
#define VECTOR_FLOAT_BOOST_TOP_NODE_GUARD

#include <iostream>
#include <string>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "vector_float_boost.h"
#include "return_map_T.h"
// namespace mp = boost::multiprecision;

namespace CFL_OBDD {
	typedef boost::multiprecision::cpp_dec_float_100 BIG_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef ReturnMapBody<BIG_FLOAT> FloatBoostReturnMapBody;
	typedef ReturnMapHandle<BIG_FLOAT> FloatBoostReturnMapHandle;
}
#include "connectionT.h"
namespace CFL_OBDD {
	typedef ConnectionT<CFLOBDDReturnMapHandle> Connection;
}

namespace CFL_OBDD {

	typedef CFLOBDDTopNodeT<BIG_FLOAT> CFLOBDDTopNodeFloatBoost;
	typedef CFLOBDDTopNodeT<BIG_FLOAT>::CFLOBDDTopNodeTRefPtr CFLOBDDTopNodeFloatBoostRefPtr;

	namespace VectorFloatBoost {

		// Initialization routine
		extern void VectorInitializerTop();

		extern CFLOBDDTopNodeFloatBoostRefPtr MkBasisVectorTop(unsigned int level, unsigned int index); // Representation of basis vector at index
		extern CFLOBDDTopNodeFloatBoostRefPtr MkBasisVectorTop(unsigned int level, std::string s); // Representation of basis vector at index
		extern CFLOBDDTopNodeFloatBoostRefPtr VectorToMatrixInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n); // Convert vector to matrix with variables interleaved
		extern CFLOBDDTopNodeFloatBoostRefPtr MkColumn1MatrixTop(unsigned int level);
		extern CFLOBDDTopNodeFloatBoostRefPtr MkVectorWithVoc12Top(CFLOBDDTopNodeFloatBoostRefPtr n);
		extern CFLOBDDTopNodeFloatBoostRefPtr VectorShiftVocs1To2Top(CFLOBDDTopNodeFloatBoostRefPtr n);
		extern CFLOBDDTopNodeFloatBoostRefPtr MatrixToVectorTop(CFLOBDDTopNodeFloatBoostRefPtr c);
		extern CFLOBDDTopNodeFloatBoostRefPtr NoDistinctionNodeTop(unsigned int level, int val);
		extern CFLOBDDTopNodeFloatBoostRefPtr ConvertToDoubleTop(CFLOBDDTopNodeFloatBoostRefPtr c);
//#ifdef PATH_COUNTING_ENABLED
		extern std::string SamplingTop(CFLOBDDTopNodeFloatBoostRefPtr n, bool voctwo = false, std::string func = "");
		extern std::string SamplingV2Top(CFLOBDDTopNodeFloatBoostRefPtr n);
//#endif
		extern CFLOBDDTopNodeFloatBoostRefPtr VectorWithAmplitudeTop(CFLOBDDTopNodeFloatBoostRefPtr n);
		extern void VectorPrintColumnMajorTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);
		extern void VectorPrintColumnMajorInterleavedTop(CFLOBDDTopNodeFloatBoostRefPtr n, std::ostream & out);
	}
}

#endif

